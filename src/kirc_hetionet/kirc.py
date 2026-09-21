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
# A standardized matrix from a prior run wins over the raw Xena download.
KIRC_EXPRESSION_CANDIDATES = [
    "kirc_expression_standardized.tsv.gz",
    "kirc_expression_standardized.tsv",
    "kirc_expression.tsv.gz",
    config.KIRC_EXPRESSION_FILE_NAME,
]


# ---------------------------------------------------------------------------
# Download (carried over from the original notebook's fetch_xena)
# ---------------------------------------------------------------------------
def fetch_xena(fname: str, dest: Path | None = None, min_bytes: int = 100_000) -> Path:
    """Download ``fname`` from the UCSC Xena GDC hub into data/raw/.

    Endpoints are tried in order: the official host first, then the S3 bucket
    serving the same data (used on networks where the official host is
    blocked). An already-downloaded file is kept.
    """
    import urllib.error
    import urllib.request

    dest = (config.RAW_DIR / fname) if dest is None else Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > min_bytes:
        print(f"  [skip] {fname} already present ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest

    errors = []
    for base in config.XENA_GDC_HOSTS:
        url = f"{base}/{fname}"
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"  [ok  ] {fname} ({dest.stat().st_size / 1e6:.1f} MB)")
            print(f"         <- {base}")
            return dest
        except (urllib.error.URLError, OSError) as exc:
            print(f"  [fail] {type(exc).__name__}: {url}")
            errors.append(f"{url}: {exc}")

    raise RuntimeError(
        f"Every Xena endpoint failed for {fname}:\n  " + "\n  ".join(errors)
    )


def download_kirc(include_clinical: bool = True) -> dict[str, Path]:
    """Download the TCGA-KIRC expression matrix (and clinical table)."""
    out = {"expression": fetch_xena(config.KIRC_EXPRESSION_FILE_NAME)}
    if include_clinical:
        out["clinical"] = fetch_xena(
            config.KIRC_CLINICAL_FILE_NAME, min_bytes=10_000
        )
    return out


def load_kirc_clinical(path: Path | None = None):
    """Load the KIRC clinical table, if it has been downloaded."""
    path = config.KIRC_CLINICAL_FILE if path is None else Path(path)
    if not path.exists():
        return None
    clin = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    print(f"[load] KIRC clinical {clin.shape[0]:,} samples x {clin.shape[1]:,} fields")
    return clin


def sample_type_counts(expr: pd.DataFrame) -> pd.DataFrame:
    """TCGA barcode field 4 -> sample type distribution (original step 10)."""
    codes = pd.Series(expr.columns, dtype=str).str.split("-").str[3].str[:2]
    counts = codes.value_counts().rename_axis("code").reset_index(name="n_samples")
    counts["sample_type"] = counts["code"].map(config.TCGA_SAMPLE_TYPES).fillna("?")
    return counts[["code", "sample_type", "n_samples"]]


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


def load_kirc_expression(
    path: Path | None = None,
    force_reload: bool = False,
    download_if_missing: bool = True,
):
    """Load the KIRC expression matrix (genes x samples).

    Returns ``(dataframe, path)``, or ``(None, None)`` when no matrix is
    reachable - the caller then falls back to a gene-list-only source.
    """
    if "expression" in _CACHE and not force_reload and path is None:
        return _CACHE["expression"]  # type: ignore[return-value]

    path = find_kirc_expression() if path is None else Path(path)
    if path is None and download_if_missing:
        print("[get ] no local KIRC matrix - downloading from the Xena GDC hub")
        try:
            path = download_kirc()["expression"]
        except (RuntimeError, OSError) as exc:
            print(f"[FAIL] KIRC download failed: {exc}")
            path = None
    if path is None:
        print("[miss] no KIRC expression matrix found in data/, Drive or Xena")
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


def resolve_kirc_genes(download_if_missing: bool = True):
    """Resolve the KIRC gene universe from the best available source.

    Returns ``{"ensembl_ids", "hetionet_gene_ids", "source", "path"}``.
    ``ensembl_ids`` is ``None`` when only the Hetionet-side fallback is
    available.
    """
    expr, path = load_kirc_expression(download_if_missing=download_if_missing)
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
