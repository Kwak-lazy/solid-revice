"""Turn a pipeline run into README sections.

Kept in the library (not in ``scripts/``) so the notebook and the headless
runner write the project log through exactly the same code.
"""

from __future__ import annotations

import config

from . import readme_log


def update_readme(results, comparison, std, kirc_source, mapping_path, issues) -> None:
    """Record this run in the README (section 4-11)."""
    ran_contexts = [c for c in config.ALL_GRAPH_CONTEXTS if c in results]

    # --- 4. Gene ID standardization ---
    counts = (std or {}).get("validation", {}).get("counts", {})
    source = (kirc_source or {}).get("source", "")
    authoritative = str((std or {}).get("source", "")).startswith("collaborator")

    lines = [
        "Join key: `hetionet_gene_id` (`Gene::<entrez>`). Gene symbols are "
        "annotation only and are never used to join.",
        "",
    ]
    if authoritative:
        lines += [
            "**The gene mapping is not produced by this pipeline.** It is the "
            "collaborator's standardization run (2026-09-15), used as-is and "
            "never regenerated or overwritten.",
            "",
            "```text",
            "KIRC Ensembl ID (GENCODE v36, versioned)",
            "   |  1st: GENCODE gene-line hgnc_id      38,891 features",
            "   |  2nd: HGNC ensembl_gene_id lookup     3,177 features",
            "   v",
            "current HGNC symbol + entrez_id",
            "   |  \"Gene::\" + entrez, only if the node exists in Hetionet",
            "   v",
            "Hetionet Gene::Entrez",
            "```",
            "",
            "Routing through **`hgnc_id`** rather than Ensembl ID is the whole "
            "point. HGNC IDs are curated and stable; Ensembl gene IDs get "
            "reassigned between builds. TCGA is frozen at GENCODE v36 while "
            "HGNC tracks the current release, so an Ensembl-to-Ensembl join "
            "silently drops every gene whose ID moved - SOD2 "
            "(`ENSG00000112096` -> `ENSG00000291237`) among them.",
            "",
            "### Two gene universes - not interchangeable",
            "",
            "| Universe | Genes | Use for |",
            "|---|---:|---|",
            "| `graph` | 19,425 | graph-only work (topology, DWPC, metapaths) |",
            "| `expression` | 19,297 | anything needing expression values |",
            "",
            "The 128-gene difference is nodes whose KIRC row is zero across all "
            "610 samples. Selected via `--gene-universe` / `config.GENE_UNIVERSE`.",
            "",
            "### Source files (read-only, `data/external/collab/`)",
            "",
            "| File | Content |",
            "|---|---|",
            "| `kirc_gene_mapping_all.tsv` | 60,660 features x 36 columns |",
            "| `kirc_hetionet_gene_nodes.tsv` | 19,425 connected Hetionet nodes |",
            "| `hetionet_genes_absent_from_kirc.tsv` | 1,520 nodes with no KIRC feature |",
            "| `kirc_gene_standardization_summary.md` | rules, counts, SHA-256 of every input |",
            "",
            "Verified on load: feature / graph-node / expression-node / absent "
            "counts, plus disjointness and union against Hetionet's 20,945 Gene "
            "nodes. A mismatch raises rather than proceeding.",
            "",
            "| Counter | Value |",
            "|---|---:|",
        ]
        for key, value in counts.items():
            lines.append(f"| {key} | {value:,} |")
    else:
        lines += [
            "> **Provisional.** The collaborator's standardization was not "
            "available, so the in-house Ensembl -> HGNC -> Entrez join ran "
            "instead. It cannot follow Ensembl IDs reassigned since GENCODE "
            "v36 and loses genes. Do not build on these numbers.",
            "",
            f"Source used: `{source}`",
            "",
        ]
        if counts:
            lines += ["| Counter | Value |", "|---|---:|"]
            lines += [f"| {k} | {v:,} |" for k, v in counts.items()]

    readme_log.set_section("4. Gene ID Standardization", "\n".join(lines))

    # --- 5. Graph Context ---
    ctx_lines = [
        "A single switch selects the context; there is no second pipeline.",
        "",
        "```python",
        f"GRAPH_CONTEXT = {config.GRAPH_CONTEXT!r}   # or "
        f"{[c for c in config.ALL_GRAPH_CONTEXTS if c != config.GRAPH_CONTEXT][0]!r}",
        "```",
        "",
        "| Context | node_kind | metaedge |",
        "|---|---|---|",
    ]
    for name, cfg in config.CONTEXT_CONFIG.items():
        ctx_lines.append(f"| `{name}` | `{cfg['node_kind']}` | `{cfg['metaedge']}` |")
    ctx_lines += [
        "",
        "```text",
        "Gene --GpPW--> Pathway",
        "Gene --GpBP--> Biological Process",
        "```",
        "",
        "Both contexts emit the same schema: "
        "`gene_id`, `gene_symbol`, `context_id`, `context_name`, `edge_type`, "
        "`context_type`.",
    ]
    readme_log.set_section("5. Graph Context", "\n".join(ctx_lines))

    # --- 7. Results ---
    notes = ["", "Notes:", ""]
    notes.append(
        "- `Gene-context edges` counts every edge of the metaedge in Hetionet; "
        "`Subgraph edges` counts only those whose gene is KIRC-mapped."
    )
    notes.append(
        "- `Subgraph nodes` / `Subgraph edges` are the NetworkX subgraph after "
        "isolate removal: Gene + Disease + the context's nodes, joined by "
        "`GiG`, `Gr>G`, `DaG`, `DuG`, `DdG` and the context metaedge. Gene nodes "
        "are restricted to the KIRC-mapped `Gene::Entrez` set."
    )
    if kirc_source.get("source", "").startswith("prior_run"):
        notes.append(
            "- `KIRC mapped genes` on this run came from the prior-run coverage "
            "artifact, not from a freshly loaded expression matrix (see section 8)."
        )
    readme_log.update_results(comparison, extra_lines=notes)

    # --- 6. Progress ---
    readme_log.log_progress(
        [
            (True, "Gene ID standardization implemented "
                   "(Ensembl -> HGNC -> Entrez -> `Gene::Entrez`)"),
            (True, "Hetionet Gene mapping validation implemented and run"),
            (True, "Pathway / Biological Process configuration "
                   "(`GRAPH_CONTEXT` + `CONTEXT_CONFIG`)"),
            (True, "Common pipeline implementation "
                   "(`get_context_nodes` / `get_context_edges` / "
                   "`run_context_experiment`)"),
            ("pathway" in ran_contexts, "Pathway experiment executed"),
            ("biological_process" in ran_contexts,
             "Biological Process experiment executed"),
            (True, "Automatic context comparison written to "
                   "`results/comparison/context_comparison.tsv`"),
        ]
    )

    # --- 8. Issues ---
    readme_log.log_issues(issues)

    # --- 9. Decisions ---
    readme_log.set_decisions(
        [
            "Gene graph join key = `Gene::Entrez` (`hetionet_gene_id`); "
            "gene symbols are annotation only.",
            "Pathway edge = `GpPW`; Biological Process edge = `GpBP`. "
            "Nothing else is read for a context.",
            "Pathway and Biological Process share one pipeline parameterised by "
            "`GRAPH_CONTEXT`, never two code paths.",
            "Results are written per context under `results/<context>/`, so the two "
            "runs cannot overwrite each other.",
            "The gene mapping is the collaborator's standardization run, used "
            "as-is. This pipeline does not re-derive it and never writes to "
            "those files.",
            "Ensembl -> Entrez is bridged by GENCODE v36 `hgnc_id`, not by "
            "Ensembl ID and not by gene symbol.",
            "`graph` (19,425) and `expression` (19,297) are separate gene "
            "universes; every result states which one it used.",
            "Hetionet edges are pulled from the git-lfs media endpoint and "
            "checksum-verified; the plain raw URL serves a 133-byte LFS pointer.",
        ]
    )

    # --- 10. Next steps ---
    readme_log.set_next_steps(
        [
            "Candidate selection",
            "Variance calculation",
            "Pathway-based scoring",
            "Biological Process-based scoring",
            "Candidate ranking",
            "Machine learning baseline",
            "GNN model",
            "Biomarker selection / disease-based scoring",
        ]
    )

    # --- Current status ---
    readme_log.update_status(
        stage="Gene ID + Graph Context Pipeline",
        status="In Progress",
        completed=[
            "Gene ID standardization",
            "Pathway / Biological Process switch",
            "Automated dual-context pipeline",
            "Automatic context comparison",
            "README progress logging",
            "Google Drive backup structure",
        ],
        next_steps=["Candidate selection", "Variance-based selection"],
    )

    # --- 11. Change log ---
    readme_log.log_change(
        [
            "Added `config.py` with `GRAPH_CONTEXT` / `CONTEXT_CONFIG`",
            "Added Gene ID standardization (Ensembl -> HGNC -> Entrez -> Hetionet)",
            "Added gene mapping validation and issue reporting",
            "Added the shared Gene -> Context pipeline and dual-context driver",
            "Added per-context result output and automatic comparison",
            "Added README progress tracking and Google Drive backup",
        ]
    )
