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
    parser.add_argument(
        "--gene-universe",
        choices=["graph", "expression"],
        default=config.GENE_UNIVERSE,
        help="graph = 19,425 Hetionet-connected nodes; "
             "expression = 19,297 of those with a non-zero expression row",
    )
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
    std = gene_mapping.standardize_gene_ids(
        kirc_source["ensembl_ids"], gene_universe=args.gene_universe
    )
    kirc_gene_ids = std["gene_ids"]
    mapping_path = std["source"]

    banner(7, "Gene mapping validation")
    if std["source"].startswith("collaborator"):
        print("  Verified on load against the expected counts in")
        print(f"  {config.COLLAB_SUMMARY_FILE.name}: features / graph nodes /")
        print("  expression nodes / absent nodes, plus the disjointness and")
        print("  union checks against Hetionet's 20,945 Gene nodes. All passed.")
        print(f"\n  gene_universe = {args.gene_universe!r} -> "
              f"{len(kirc_gene_ids):,} genes")
        if args.gene_universe == "graph":
            print("  (128 of these have no non-zero expression row; use")
            print("   --gene-universe expression for anything needing values)")
    else:
        for issue in std["validation"]["issues"]:
            issues.append({
                "what": issue,
                "where": "gene_mapping.validate_gene_mapping()",
                "cause": "Ensembl/Entrez identifier systems are not 1:1",
                "status": "recorded, not auto-corrected",
            })
        issues.append({
            "what": "In-house gene mapping used instead of the collaborator's "
                    "standardization",
            "where": "gene_mapping.standardize_gene_ids() fallback branch",
            "cause": f"{config.COLLAB_MAPPING_FILE} not present",
            "status": "open - results are provisional until the authoritative "
                      "files are restored",
        })

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
