"""Central configuration for the KIRC-Hetionet project.

Everything that the notebook, the library code and the scripts need to agree on
lives here: filesystem layout, download URLs, the Google Drive target and -
most importantly - the single ``GRAPH_CONTEXT`` switch that selects between the
Pathway and the Biological Process view of the graph.
"""

from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Graph context switch
# ---------------------------------------------------------------------------
# Change this single value to run the whole pipeline in the other context.
# It is also the default used by ``run_context_experiment()`` when no explicit
# context is passed.
GRAPH_CONTEXT = "pathway"

CONTEXT_CONFIG = {
    "pathway": {
        "node_kind": "Pathway",
        "metaedge": "GpPW",
    },
    "biological_process": {
        "node_kind": "Biological Process",
        "metaedge": "GpBP",
    },
}

# Every context the dual-run driver executes, in order.
ALL_GRAPH_CONTEXTS = ["pathway", "biological_process"]

# ---------------------------------------------------------------------------
# 1b. Subgraph extraction (carried over from the original notebook)
# ---------------------------------------------------------------------------
# The original notebook kept 4 node types and 7 edge types at once. The only
# change here is that the Gene->Context layer is no longer hard-wired to both
# GpPW and GpBP: the selected context contributes exactly one of them, while
# everything below stays as it was.
BASE_NODE_KINDS = ["Gene", "Disease"]

BASE_EDGE_TYPES = [
    "GiG",    # Gene - interacts - Gene
    "Gr>G",   # Gene - regulates -> Gene
    "DaG",    # Disease - associates - Gene
    "DuG",    # Disease - upregulates - Gene
    "DdG",    # Disease - downregulates - Gene
]


def keep_node_kinds(graph_context: str) -> set:
    """Node kinds kept in the subgraph for a context (original: 4 kinds)."""
    return set(BASE_NODE_KINDS) | {CONTEXT_CONFIG[graph_context]["node_kind"]}


def keep_edge_types(graph_context: str) -> set:
    """Edge types kept in the subgraph for a context (original: 7 types)."""
    return set(BASE_EDGE_TYPES) | {CONTEXT_CONFIG[graph_context]["metaedge"]}


# Full metaedge glossary from the original notebook, kept for the edge-type
# report in step 4.
METAEDGE_NAMES = {
    "GpBP": "Gene -> participates -> Biological Process",
    "AeG": "Anatomy -> expresses -> Gene",
    "Gr>G": "Gene -> regulates -> Gene",
    "GiG": "Gene -> interacts -> Gene",
    "CcSE": "Compound -> causes -> Side Effect",
    "AdG": "Anatomy -> downregulates -> Gene",
    "AuG": "Anatomy -> upregulates -> Gene",
    "GpMF": "Gene -> participates -> Molecular Function",
    "GpPW": "Gene -> participates -> Pathway",
    "GpCC": "Gene -> participates -> Cellular Component",
    "GcG": "Gene -> covaries -> Gene",
    "CdG": "Compound -> downregulates -> Gene",
    "CuG": "Compound -> upregulates -> Gene",
    "DaG": "Disease -> associates -> Gene",
    "CbG": "Compound -> binds -> Gene",
    "DuG": "Disease -> upregulates -> Gene",
    "DdG": "Disease -> downregulates -> Gene",
    "CrC": "Compound -> resembles -> Compound",
    "DlA": "Disease -> localizes -> Anatomy",
    "DpS": "Disease -> presents -> Symptom",
    "PCiC": "Pharmacologic Class -> includes -> Compound",
    "CtD": "Compound -> treats -> Disease",
    "DrD": "Disease -> resembles -> Disease",
    "CpD": "Compound -> palliates -> Disease",
}

# KIRC in Hetionet's Disease vocabulary (found by the original notebook).
KIRC_DISEASE_ID = "Disease::DOID:263"
KIDNEY_DISEASE_KEYWORDS = ["kidney", "renal", "clear cell"]


# ---------------------------------------------------------------------------
# 2. Project layout
# ---------------------------------------------------------------------------
def _detect_project_dir() -> Path:
    """Project root: env override > this file's directory."""
    override = os.environ.get("KIRC_HETIONET_PROJECT_DIR")
    if override:
        return Path(override).expanduser().resolve()
    return Path(__file__).resolve().parent


PROJECT_DIR = _detect_project_dir()

NOTEBOOK_DIR = PROJECT_DIR / "notebooks"
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"
RESULTS_DIR = PROJECT_DIR / "results"
COMPARISON_DIR = RESULTS_DIR / "comparison"
README_PATH = PROJECT_DIR / "README.md"
NOTEBOOK_NAME = "Download_and_subgraph.ipynb"
NOTEBOOK_PATH = NOTEBOOK_DIR / NOTEBOOK_NAME


def results_dir(graph_context: str) -> Path:
    """results/<graph_context>/ - one directory per context so the two runs
    can never overwrite each other."""
    if graph_context not in CONTEXT_CONFIG:
        raise KeyError(
            f"Unknown graph_context {graph_context!r}. "
            f"Expected one of {sorted(CONTEXT_CONFIG)}."
        )
    return RESULTS_DIR / graph_context


def ensure_dirs() -> None:
    """Create the project directory tree if it does not exist yet."""
    dirs = [NOTEBOOK_DIR, RAW_DIR, PROCESSED_DIR, EXTERNAL_DIR, COMPARISON_DIR]
    dirs += [results_dir(ctx) for ctx in ALL_GRAPH_CONTEXTS]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 3. Downloads
# ---------------------------------------------------------------------------
# Hetionet v1.0.  The edge file is stored with git-lfs, so it must be pulled
# from the media endpoint - the plain raw URL returns a 133-byte LFS pointer.
HETIONET_NODES_URL = (
    "https://raw.githubusercontent.com/hetio/hetionet/main/"
    "hetnet/tsv/hetionet-v1.0-nodes.tsv"
)
HETIONET_EDGES_URL = (
    "https://media.githubusercontent.com/media/hetio/hetionet/main/"
    "hetnet/tsv/hetionet-v1.0-edges.sif.gz"
)
HETIONET_EDGES_SHA256 = (
    "611f0411ac666be4e0270e54bbca67407e14720bfc3099311f0f5eec25c5a947"
)

HETIONET_NODES_FILE = RAW_DIR / "hetionet-v1.0-nodes.tsv"
HETIONET_EDGES_FILE = RAW_DIR / "hetionet-v1.0-edges.sif.gz"

# The original notebook wrote data/nodes.tsv and data/edges.sif (uncompressed).
# Those are still accepted so an existing checkout is not re-downloaded.
HETIONET_NODES_ALT = [DATA_DIR / "nodes.tsv", RAW_DIR / "nodes.tsv"]
HETIONET_EDGES_ALT = [DATA_DIR / "edges.sif", RAW_DIR / "edges.sif"]

# HGNC complete set: the Ensembl -> HGNC -> Entrez bridge.
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)
HGNC_FILE = RAW_DIR / "hgnc_complete_set.txt"


# TCGA-KIRC from the UCSC Xena GDC hub.  Two endpoints are tried in order:
# the official host first, then the S3 bucket that serves the same data on
# networks where the official host is blocked.
XENA_GDC_HOSTS = [
    "https://gdc.xenahubs.net/download",
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download",
]
KIRC_EXPRESSION_FILE_NAME = "TCGA-KIRC.star_tpm.tsv.gz"  # log2(TPM + 1)
KIRC_CLINICAL_FILE_NAME = "TCGA-KIRC.clinical.tsv.gz"
KIRC_EXPRESSION_FILE = RAW_DIR / KIRC_EXPRESSION_FILE_NAME
KIRC_CLINICAL_FILE = RAW_DIR / KIRC_CLINICAL_FILE_NAME

# TCGA barcode field 4 = sample type.
TCGA_SAMPLE_TYPES = {
    "01": "Primary Tumor",
    "05": "Additional - New Primary",
    "06": "Metastatic",
    "11": "Solid Tissue Normal",
}


# ---------------------------------------------------------------------------
# 4. Gene ID standardization
# ---------------------------------------------------------------------------
# The graph join key.  Gene symbols are annotation only - never a join key.
GENE_JOIN_KEY = "hetionet_gene_id"  # e.g. "Gene::7157"
HETIONET_GENE_PREFIX = "Gene::"

GENE_MAPPING_COLUMNS = [
    "ensembl_id",
    "gene_symbol",
    "entrez_id",
    "hetionet_gene_id",
    "hetionet_gene_symbol",
    "mapping_status",
]

# Re-used before anything is regenerated (checked in this order, in
# PROCESSED_DIR, EXTERNAL_DIR, PROJECT_DIR and the Drive project dir).
EXISTING_MAPPING_FILES = {
    "gene_mapping_all": "kirc_gene_mapping_all.tsv",
    "hetionet_gene_nodes": "kirc_hetionet_gene_nodes.tsv",
    "expression_standardized": "kirc_expression_standardized.tsv.gz",
}

# Prior-run artifact shared by a collaborator: Hetionet genes that carry no
# KIRC feature.  Used only as a documented fallback when the KIRC expression
# matrix itself is not reachable - see README section 8.
SHARED_ABSENT_GENES_FILE = EXTERNAL_DIR / "hetionet_genes_absent_from_kirc.tsv"


# ---------------------------------------------------------------------------
# 5. Google Drive
# ---------------------------------------------------------------------------
DRIVE_PROJECT_DIR = "/content/drive/MyDrive/KIRC_Hetionet_Project"
DRIVE_MOUNT_POINT = "/content/drive"

# Relative paths that ``save_project_to_drive()`` mirrors into Drive.
DRIVE_SYNC_PATHS = [
    "README.md",
    "config.py",
    "requirements.txt",
    "notebooks",
    "src",
    "scripts",
    "results",
    "data/processed",
    "data/external",
]

# Files whose presence in Drive is asserted after a backup.
DRIVE_VERIFY_PATHS = [
    ("README.md", "README.md"),
    (f"notebooks/{NOTEBOOK_NAME}", "Notebook"),
    ("results/pathway/context_nodes.tsv", "Pathway results"),
    ("results/biological_process/context_nodes.tsv", "Biological Process results"),
    ("results/comparison/context_comparison.tsv", "Comparison"),
]


def in_colab() -> bool:
    try:
        import google.colab  # noqa: F401
        return True
    except Exception:
        return False
