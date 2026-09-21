#!/usr/bin/env python3
"""Headless driver for the KIRC-Hetionet pipeline.

Runs the same 12 steps as notebooks/Download_and_subgraph.ipynb:

    1. Environment            7. Mapping validation
    2. Drive mount            8. Context configuration
    3. Project paths          9. Pathway / Biological Process pipeline
    4. Hetionet loading      10. Result comparison
    5. KIRC loading          11. README update
    6. Gene ID standardization  12. Google Drive backup

Usage
    python scripts/run_pipeline.py                      # both contexts
    python scripts/run_pipeline.py --context pathway    # a single context
    python scripts/run_pipeline.py --no-drive           # skip the backup step
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import config  # noqa: E402
from kirc_hetionet import context, drive, gene_mapping, hetionet, kirc, report  # noqa: E402


def banner(step: int, title: str) -> None:
    print(f"\n{'#' * 70}\n# [{step:02d}] {title}\n{'#' * 70}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--context",
        choices=config.ALL_GRAPH_CONTEXTS,
        help="run a single context instead of both",
    )
    parser.add_argument("--no-drive", action="store_true", help="skip the Drive backup")
    parser.add_argument("--no-readme", action="store_true", help="skip the README update")
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()

    issues: list[dict] = []

    # ---- 1-3 environment / drive / paths ---------------------------------
    banner(1, "Environment and project paths")
    config.ensure_dirs()
    print(f"project dir : {config.PROJECT_DIR}")
    print(f"results dir : {config.RESULTS_DIR}")
    print(f"colab       : {config.in_colab()}")

    banner(2, "Google Drive mount")
    drive_ok = drive.mount_drive()

    # ---- 4 Hetionet -------------------------------------------------------
    banner(4, "Hetionet download and loading")
    hetionet.download_hetionet(force=args.force_download)
    nodes = hetionet.load_nodes()
    edges = hetionet.load_edges()
    hetionet.describe_hetionet(nodes, edges)
    het_genes = hetionet.gene_nodes(nodes)

    # ---- 5 KIRC -----------------------------------------------------------
    banner(5, "KIRC loading")
    kirc_source = kirc.resolve_kirc_genes()
    print(f"KIRC gene source: {kirc_source['source']}")

    # ---- 6-7 gene ID standardization + validation -------------------------
    banner(6, "Gene ID standardization (Ensembl -> HGNC -> Entrez -> Gene::Entrez)")
    if kirc_source["ensembl_ids"] is not None:
        std = gene_mapping.standardize_gene_ids(kirc_source["ensembl_ids"])
        kirc_gene_ids = sorted(set(std["usable"]["hetionet_gene_id"].dropna()))
        mapping_path = "ensembl -> hgnc -> entrez -> hetionet"
    elif kirc_source["hetionet_gene_ids"] is not None:
        print(
            "[warn] KIRC expression matrix not reachable in this environment.\n"
            "       Falling back to the Hetionet-side gene coverage established\n"
            "       by a prior run; the Ensembl -> Entrez leg is NOT exercised\n"
            "       on this path, so its counters are reported as n/a."
        )
        issues.append(
            {
                "what": "KIRC expression matrix not reachable; gene coverage taken "
                        "from a prior run instead",
                "where": "scripts/run_pipeline.py step 5-6, kirc.resolve_kirc_genes()",
                "cause": "TCGA/GDC hosts are blocked by this environment's egress "
                         "policy and the Drive copies are shortcuts that the Drive "
                         "API refuses to download",
                "status": "open - re-run steps 6-7 in Colab with the mounted "
                          "kirc_expression_standardized.tsv.gz to get the real "
                          "Ensembl/HGNC/Entrez counters",
            }
        )
        mapping = gene_mapping.mapping_from_hetionet_gene_nodes(
            kirc_source["hetionet_gene_ids"]
        )
        banner(7, "Gene mapping validation")
        validation = gene_mapping.validate_gene_mapping(
            mapping,
            het_genes=het_genes,
            n_kirc_genes=len(kirc_source["hetionet_gene_ids"]),
        )
        out = config.PROCESSED_DIR / "kirc_hetionet_gene_nodes_from_prior_run.tsv"
        mapping.to_csv(out, sep="\t", index=False)
        print(f"[write] {out}")
        std = {"mapping": mapping, "usable": mapping, "validation": validation,
               "source": kirc_source["source"], "path": kirc_source["path"]}
        kirc_gene_ids = sorted(set(mapping["hetionet_gene_id"].dropna()))
        mapping_path = "prior run (Hetionet side only)"
    else:
        print("[FAIL] no KIRC gene source; the subgraph cannot be restricted to KIRC")
        std, kirc_gene_ids, mapping_path = None, None, "unavailable"
        issues.append(
            {
                "what": "No KIRC gene source available",
                "where": "step 5",
                "cause": "neither an expression matrix nor a prior-run artifact "
                         "was found",
                "status": "open",
            }
        )

    if std is not None and kirc_source["ensembl_ids"] is not None:
        banner(7, "Gene mapping validation")
        for issue in std["validation"]["issues"]:
            issues.append(
                {
                    "what": issue,
                    "where": "gene_mapping.validate_gene_mapping()",
                    "cause": "Ensembl/Entrez identifier systems are not 1:1",
                    "status": "recorded, not auto-corrected",
                }
            )

    # ---- 8-10 context pipeline + comparison -------------------------------
    banner(8, "Context configuration")
    for name, cfg in config.CONTEXT_CONFIG.items():
        print(f"  {name:<20} node_kind={cfg['node_kind']!r:<22} metaedge={cfg['metaedge']!r}")
    print(f"  default GRAPH_CONTEXT = {config.GRAPH_CONTEXT!r}")

    banner(9, "Gene -> Context pipeline")
    if args.context:
        single = context.run_context_experiment(
            args.context, kirc_gene_ids=kirc_gene_ids, nodes=nodes, edges=edges
        )
        results = {args.context: single}
        comparison = context.build_comparison(results, write=False)
    else:
        results = context.run_all_context_experiments(
            kirc_gene_ids=kirc_gene_ids, nodes=nodes, edges=edges
        )
        comparison = results["comparison"]

    banner(10, "Result comparison")
    print(context.format_comparison(comparison))

    # ---- 11 README --------------------------------------------------------
    if not args.no_readme:
        banner(11, "README update")
        report.update_readme(results, comparison, std, kirc_source, mapping_path, issues)
        print(f"[ok  ] {config.README_PATH}")

    # ---- 12 Drive ---------------------------------------------------------
    if not args.no_drive:
        banner(12, "Google Drive backup")
        if drive_ok or drive.drive_available():
            drive.save_project_to_drive()
            print()
            drive.verify_drive_backup()
        else:
            print("[skip] Drive not mounted; nothing was backed up")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
