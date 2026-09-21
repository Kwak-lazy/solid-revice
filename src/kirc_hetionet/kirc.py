"""KIRC expression loading.

The pipeline only needs two things from KIRC in this stage:
  * the gene id universe (Ensembl IDs) that feeds gene-ID standardization, and
  * optionally the expression matrix itself, for later stages.

Loading is deliberately tolerant: whichever of the known KIRC artifacts is
present gets used, and the resolved source is reported so the README can record
exactly what the numbers came from.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import config

from . import gene_mapping

_CACHE: dict[str, object] = {}

# Filenames searched for a KIRC expression matrix, most-preferred first.
KIRC_EXPRESSION_CANDIDATES = [
    "kirc_expression_standardized.tsv.gz",
    "kirc_expression_standardized.tsv",
    "kirc_expression.tsv.gz",
    "TCGA-KIRC.htseq_fpkm.tsv.gz",
]


def find_kirc_expression(extra_dirs=()) -> Path | None:
    for name in KIRC_EXPRESSION_CANDIDATES:
        path = gene_mapping.search_existing_file(name, extra_dirs)
        if path is not None:
            return path
    # Last resort: an expression-looking file in data/.  The pattern is
    # deliberately narrow - a loose "*kirc*" glob happily matches this
    # pipeline's own mapping outputs and would feed them back in as if they
    # were an expression matrix.
    for d in (config.RAW_DIR, config.PROCESSED_DIR):
        for path in sorted(d.glob("*kirc*expression*")):
            if path.suffix in {".tsv", ".gz", ".csv"}:
                return path
    return None


def load_kirc_expression(path: Path | None = None, force_reload: bool = False):
    """Load the KIRC expression matrix (genes x samples).

    Returns ``(dataframe, path)``, or ``(None, None)`` when no matrix is
    reachable - the caller then falls back to a gene-list-only source.
    """
    if "expression" in _CACHE and not force_reload and path is None:
        return _CACHE["expression"]  # type: ignore[return-value]

    path = find_kirc_expression() if path is None else Path(path)
    if path is None:
        print("[miss] no KIRC expression matrix found in data/ or Drive")
        return None, None

    print(f"[load] KIRC expression from {path}")
    expr = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    expr.index = expr.index.astype(str)
    print(f"[ok  ] KIRC expression {expr.shape[0]:,} genes x {expr.shape[1]:,} samples")
    _CACHE["expression"] = (expr, path)
    return expr, path


def kirc_ensembl_ids(expr: pd.DataFrame | None = None):
    """The Ensembl gene id universe of the KIRC matrix (versions stripped)."""
    if expr is None:
        expr, _ = load_kirc_expression()
    if expr is None:
        return None
    ids = gene_mapping.strip_ensembl_version(pd.Series(expr.index, dtype=str))
    ids = ids[ids.notna() & ~ids.isin({"", "nan", "None"})]
    return sorted({str(i) for i in ids})


# ---------------------------------------------------------------------------
# Fallback: prior-run gene coverage shared by a collaborator
# ---------------------------------------------------------------------------
def kirc_hetionet_genes_from_shared(path: Path | None = None):
    """Hetionet genes that KIRC *does* cover, derived from a prior run.

    ``hetionet_genes_absent_from_kirc.tsv`` lists the Hetionet Gene nodes that
    have no KIRC feature.  Its complement against the Hetionet Gene node set is
    therefore exactly the KIRC n Hetionet intersection established by that run.

    This is a documented fallback for environments where the KIRC expression
    matrix itself is not reachable.  It does NOT exercise the
    Ensembl -> Entrez leg, so the HGNC/Entrez counters are not meaningful on
    this path and are reported as such.
    """
    from . import hetionet

    path = config.SHARED_ABSENT_GENES_FILE if path is None else Path(path)
    if not path.exists():
        return None, None

    absent = pd.read_csv(path, sep="\t", dtype=str)
    het = hetionet.gene_nodes()
    present = het[~het["hetionet_gene_id"].isin(set(absent["hetionet_gene_id"]))].copy()
    print(
        f"[load] KIRC gene coverage derived from prior run {path.name}: "
        f"{len(het):,} Hetionet genes - {absent['hetionet_gene_id'].nunique():,} absent "
        f"= {len(present):,} covered by KIRC"
    )
    return present, path


def resolve_kirc_genes():
    """Resolve the KIRC gene universe from the best available source.

    Returns ``{"ensembl_ids", "hetionet_gene_ids", "source", "path"}``.
    ``ensembl_ids`` is ``None`` when only the Hetionet-side fallback is
    available.
    """
    expr, path = load_kirc_expression()
    if expr is not None:
        return {
            "ensembl_ids": kirc_ensembl_ids(expr),
            "hetionet_gene_ids": None,
            "source": "kirc_expression_matrix",
            "path": path,
        }

    present, shared_path = kirc_hetionet_genes_from_shared()
    if present is not None:
        return {
            "ensembl_ids": None,
            "hetionet_gene_ids": sorted(set(present["hetionet_gene_id"])),
            "source": "prior_run:hetionet_genes_absent_from_kirc.tsv",
            "path": shared_path,
        }

    print("[miss] no KIRC gene source available at all")
    return {"ensembl_ids": None, "hetionet_gene_ids": None, "source": None, "path": None}
