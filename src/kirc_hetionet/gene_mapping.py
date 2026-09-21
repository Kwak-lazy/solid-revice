"""KIRC <-> Hetionet gene ID standardization.

    KIRC Ensembl gene ID
            |  HGNC complete set (annotation)
            v
        Entrez gene ID
            |  "Gene::" + entrez
            v
      Hetionet Gene::Entrez

The final graph join key is ``hetionet_gene_id`` (``Gene::<entrez>``).
``gene_symbol`` is annotation only and is never used to join.

Existing mapping artifacts (``kirc_gene_mapping_all.tsv``,
``kirc_hetionet_gene_nodes.tsv``, ``kirc_expression_standardized.tsv.gz``) are
re-used when they are found and usable; they are never overwritten in place.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import config

from . import hetionet

_CACHE: dict[str, object] = {}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def strip_ensembl_version(series: pd.Series) -> pd.Series:
    """``ENSG00000134086.7`` -> ``ENSG00000134086``."""
    return series.astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)


def _clean_entrez(series: pd.Series) -> pd.Series:
    """Entrez ids arrive as floats/str/NaN depending on the source; normalise
    everything to a bare integer string or ``pd.NA``."""
    s = pd.to_numeric(series, errors="coerce")
    out = s.astype("Int64").astype(str)
    return out.where(s.notna(), pd.NA)


def search_existing_file(filename: str, extra_dirs=()) -> Path | None:
    """Look for an already-produced artifact before regenerating anything."""
    candidates = [
        config.PROCESSED_DIR,
        config.EXTERNAL_DIR,
        config.DATA_DIR,
        config.PROJECT_DIR,
        Path(config.DRIVE_PROJECT_DIR) / "data" / "processed",
        Path(config.DRIVE_PROJECT_DIR) / "data",
        Path(config.DRIVE_PROJECT_DIR),
        *[Path(d) for d in extra_dirs],
    ]
    for d in candidates:
        try:
            p = d / filename
            if p.exists() and p.stat().st_size > 0:
                return p
        except OSError:
            continue
    return None


def find_existing_mapping_artifacts(extra_dirs=()) -> dict[str, Path | None]:
    """Report which of the three known prior artifacts are available."""
    found = {
        key: search_existing_file(name, extra_dirs)
        for key, name in config.EXISTING_MAPPING_FILES.items()
    }
    print("--- existing gene mapping artifacts ---")
    for key, name in config.EXISTING_MAPPING_FILES.items():
        path = found[key]
        print(f"  [{'FOUND ' if path else 'ABSENT'}] {name}"
              + (f"  ->  {path}" if path else ""))
    return found


# ---------------------------------------------------------------------------
# HGNC
# ---------------------------------------------------------------------------
def load_hgnc(force_reload: bool = False) -> pd.DataFrame:
    """HGNC complete set reduced to ``hgnc_id, gene_symbol, entrez_id,
    ensembl_id`` with one row per (ensembl, entrez) pair."""
    if "hgnc" in _CACHE and not force_reload:
        return _CACHE["hgnc"]  # type: ignore[return-value]

    if not config.HGNC_FILE.exists():
        hetionet.download_hgnc()

    hgnc = pd.read_csv(
        config.HGNC_FILE,
        sep="\t",
        dtype=str,
        low_memory=False,
        usecols=["hgnc_id", "symbol", "entrez_id", "ensembl_gene_id", "status"],
    )
    hgnc = hgnc[hgnc["status"] == "Approved"].copy()
    hgnc = hgnc.rename(
        columns={"symbol": "gene_symbol", "ensembl_gene_id": "ensembl_id"}
    )
    hgnc["entrez_id"] = _clean_entrez(hgnc["entrez_id"])
    hgnc["ensembl_id"] = strip_ensembl_version(hgnc["ensembl_id"]).replace(
        {"nan": pd.NA, "": pd.NA}
    )
    hgnc = hgnc.dropna(subset=["ensembl_id"])
    hgnc = hgnc[["hgnc_id", "gene_symbol", "entrez_id", "ensembl_id"]].drop_duplicates()
    _CACHE["hgnc"] = hgnc
    return hgnc


# ---------------------------------------------------------------------------
# Mapping construction
# ---------------------------------------------------------------------------
def build_gene_mapping(
    kirc_ensembl_ids,
    hgnc: pd.DataFrame | None = None,
    het_genes: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the full Ensembl -> HGNC -> Entrez -> Hetionet mapping table.

    Every input Ensembl id produces at least one row, so unmapped genes stay
    visible instead of silently disappearing.  Columns are exactly
    ``config.GENE_MAPPING_COLUMNS``.

    ``mapping_status`` is one of
      ``mapped``            - reached a Hetionet ``Gene::<entrez>`` node
      ``no_hgnc``           - Ensembl id absent from the HGNC complete set
      ``no_entrez``         - HGNC record carries no Entrez id
      ``not_in_hetionet``   - Entrez id has no Hetionet Gene node
    """
    hgnc = load_hgnc() if hgnc is None else hgnc
    het_genes = hetionet.gene_nodes() if het_genes is None else het_genes

    kirc = pd.DataFrame({"ensembl_id_raw": pd.Series(list(kirc_ensembl_ids), dtype=str)})
    kirc["ensembl_id"] = strip_ensembl_version(kirc["ensembl_id_raw"])
    kirc = kirc.drop_duplicates(subset=["ensembl_id"])

    # Ensembl -> HGNC (symbol + Entrez)
    m = kirc.merge(
        hgnc[["ensembl_id", "gene_symbol", "entrez_id"]], on="ensembl_id", how="left"
    )

    # Entrez -> Hetionet Gene::Entrez
    m = m.merge(
        het_genes[["entrez_id", "hetionet_gene_id", "hetionet_gene_symbol"]],
        on="entrez_id",
        how="left",
    )

    status = pd.Series("mapped", index=m.index, dtype="object")
    status[m["gene_symbol"].isna()] = "no_hgnc"
    status[m["gene_symbol"].notna() & m["entrez_id"].isna()] = "no_entrez"
    status[m["entrez_id"].notna() & m["hetionet_gene_id"].isna()] = "not_in_hetionet"
    m["mapping_status"] = status

    m = m[["ensembl_id"] + config.GENE_MAPPING_COLUMNS[1:]]
    return m.sort_values(["mapping_status", "ensembl_id"]).reset_index(drop=True)


def mapping_from_hetionet_gene_nodes(het_gene_ids) -> pd.DataFrame:
    """Build a mapping table when only the Hetionet side is known.

    Used for the documented fallback path where the KIRC expression matrix is
    unavailable but a prior run already established which Hetionet genes are
    covered by KIRC.  ``ensembl_id`` is left empty on purpose - those rows did
    not travel the Ensembl -> Entrez leg in this run.
    """
    het_genes = hetionet.gene_nodes()
    wanted = pd.DataFrame({"hetionet_gene_id": pd.Series(list(het_gene_ids), dtype=str)})
    out = wanted.merge(het_genes, on="hetionet_gene_id", how="left")
    hgnc = load_hgnc()
    sym = (
        hgnc.dropna(subset=["entrez_id"])
        .drop_duplicates(subset=["entrez_id"])[["entrez_id", "gene_symbol", "ensembl_id"]]
    )
    out = out.merge(sym, on="entrez_id", how="left")
    out["mapping_status"] = "mapped"
    out.loc[out["hetionet_gene_symbol"].isna(), "mapping_status"] = "not_in_hetionet"
    for col in config.GENE_MAPPING_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    return out[config.GENE_MAPPING_COLUMNS].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_gene_mapping(
    mapping: pd.DataFrame,
    kirc_ensembl_ids=None,
    het_genes: pd.DataFrame | None = None,
    n_kirc_genes: int | None = None,
) -> dict:
    """Print the required counts and return them, plus a list of issues.

    Nothing is silently repaired here: every anomaly is reported so it can be
    written into the README.
    """
    het_genes = hetionet.gene_nodes() if het_genes is None else het_genes

    if kirc_ensembl_ids is not None:
        n_kirc = len(
            set(strip_ensembl_version(pd.Series(list(kirc_ensembl_ids), dtype=str)))
        )
    elif n_kirc_genes is not None:
        # Fallback path: the KIRC gene universe was established on the Hetionet
        # side, so it must be passed in rather than counted from ensembl_id
        # (which is a reverse HGNC lookup there, not the KIRC input).
        n_kirc = int(n_kirc_genes)
    else:
        n_kirc = int(mapping["ensembl_id"].notna().sum())
    mapped = mapping[mapping["mapping_status"] == "mapped"]

    counts = {
        "KIRC genes": n_kirc,
        "HGNC mapped": int(mapping["gene_symbol"].notna().sum()),
        "Entrez mapped": int(mapping["entrez_id"].notna().sum()),
        "Hetionet mapped": int(mapping["hetionet_gene_id"].notna().sum()),
        "Unmapped": int((mapping["mapping_status"] != "mapped").sum()),
        "Duplicated": int(mapping["ensembl_id"].dropna().duplicated().sum()),
        "Final usable genes": int(mapped["hetionet_gene_id"].nunique()),
    }

    print("--- gene mapping summary ---")
    width = max(len(k) for k in counts)
    for key, value in counts.items():
        print(f"  {key:<{width}} : {value:,}")

    print("\n--- mapping_status breakdown ---")
    print(mapping["mapping_status"].value_counts().to_string())

    # ---- structural checks -------------------------------------------------
    issues: list[str] = []

    dup_ens = mapping["ensembl_id"].dropna()
    dup_ens = dup_ens[dup_ens.duplicated(keep=False)]
    if len(dup_ens):
        issues.append(
            f"{dup_ens.nunique():,} Ensembl IDs appear on more than one row "
            "(one Ensembl -> several Entrez)."
        )

    ens2entrez = (
        mapped.dropna(subset=["ensembl_id", "entrez_id"])
        .groupby("ensembl_id")["entrez_id"]
        .nunique()
    )
    multi_entrez = ens2entrez[ens2entrez > 1]
    if len(multi_entrez):
        issues.append(
            f"{len(multi_entrez):,} Ensembl IDs map to multiple Entrez IDs "
            f"(e.g. {', '.join(multi_entrez.index[:5])})."
        )

    entrez2ens = (
        mapped.dropna(subset=["ensembl_id", "entrez_id"])
        .groupby("entrez_id")["ensembl_id"]
        .nunique()
    )
    multi_ens = entrez2ens[entrez2ens > 1]
    if len(multi_ens):
        issues.append(
            f"{len(multi_ens):,} Entrez IDs map to multiple Ensembl IDs "
            f"(e.g. {', '.join(multi_ens.index[:5])})."
        )

    n_not_in_het = int((mapping["mapping_status"] == "not_in_hetionet").sum())
    if n_not_in_het:
        issues.append(
            f"{n_not_in_het:,} Entrez IDs have no Hetionet Gene node."
        )

    het_ids = set(het_genes["hetionet_gene_id"])
    usable = set(mapped["hetionet_gene_id"].dropna())
    counts["Hetionet genes"] = len(het_ids)
    counts["KIRC n Hetionet intersection"] = len(usable & het_ids)
    counts["Hetionet genes without KIRC"] = len(het_ids - usable)

    print("\n--- intersection with Hetionet Gene nodes ---")
    print(f"  Hetionet Gene nodes          : {len(het_ids):,}")
    print(f"  KIRC-mapped Hetionet genes   : {len(usable & het_ids):,}")
    print(f"  Hetionet genes without KIRC  : {len(het_ids - usable):,}")

    print("\n--- issues ---")
    if issues:
        for issue in issues:
            print(f"  [!] {issue}")
        print("  (not auto-corrected; record these in README section 8)")
    else:
        print("  none")

    return {"counts": counts, "issues": issues}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def standardize_gene_ids(
    kirc_ensembl_ids=None,
    reuse_existing: bool = True,
    write: bool = True,
) -> dict:
    """Produce the KIRC -> Hetionet gene mapping, re-using prior output first.

    Returns ``{"mapping", "usable", "validation", "source", "path"}`` where
    ``usable`` holds only the rows with ``mapping_status == "mapped"``.
    """
    config.ensure_dirs()
    artifacts = find_existing_mapping_artifacts()

    mapping = None
    source = None
    path = None

    if reuse_existing and artifacts.get("gene_mapping_all") is not None:
        path = artifacts["gene_mapping_all"]
        candidate = pd.read_csv(path, sep="\t", dtype=str)
        missing = [c for c in config.GENE_MAPPING_COLUMNS if c not in candidate.columns]
        if missing:
            print(
                f"[warn] {path.name} is missing {missing}; it is left untouched "
                "and the mapping is rebuilt instead (record the reason in README)."
            )
        else:
            mapping = candidate[config.GENE_MAPPING_COLUMNS].copy()
            source = f"reused:{path}"
            print(f"[reuse] gene mapping loaded from {path}")

    if mapping is None:
        if kirc_ensembl_ids is None:
            raise ValueError(
                "No reusable kirc_gene_mapping_all.tsv was found and no "
                "kirc_ensembl_ids were supplied - cannot build the mapping."
            )
        mapping = build_gene_mapping(kirc_ensembl_ids)
        source = "built:ensembl->hgnc->entrez->hetionet"

    validation = validate_gene_mapping(mapping, kirc_ensembl_ids=kirc_ensembl_ids)
    usable = mapping[mapping["mapping_status"] == "mapped"].copy()

    if write:
        out_map = config.PROCESSED_DIR / "kirc_gene_mapping_all.tsv"
        out_nodes = config.PROCESSED_DIR / "kirc_hetionet_gene_nodes.tsv"
        if source.startswith("reused") and out_map.resolve() == Path(path).resolve():
            print(f"[keep ] existing {out_map.name} left untouched")
        else:
            mapping.to_csv(out_map, sep="\t", index=False)
            print(f"[write] {out_map}")
        usable[["hetionet_gene_id", "hetionet_gene_symbol", "entrez_id",
                "ensembl_id", "gene_symbol"]].drop_duplicates().to_csv(
            out_nodes, sep="\t", index=False
        )
        print(f"[write] {out_nodes}")

    return {
        "mapping": mapping,
        "usable": usable,
        "validation": validation,
        "source": source,
        "path": path,
    }
