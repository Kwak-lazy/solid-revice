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
    lines = [
        "Join key: `hetionet_gene_id` (`Gene::<entrez>`). Gene symbols are "
        "annotation only and are never used to join.",
        "",
        "```text",
        "KIRC Ensembl Gene ID",
        "        v  HGNC complete set",
        "Entrez Gene ID",
        "        v  \"Gene::\" + entrez",
        "Hetionet Gene::Entrez",
        "```",
        "",
        f"Source used in the last run: `{kirc_source.get('source')}`  "
        f"(path: `{kirc_source.get('path')}`)",
        f"Mapping path exercised: {mapping_path}",
        "",
    ]
    if counts:
        lines += ["| Counter | Value |", "|---|---:|"]
        for key, value in counts.items():
            lines.append(f"| {key} | {value:,} |")
    lines += [
        "",
        "Mapping table columns: "
        "`ensembl_id`, `gene_symbol`, `entrez_id`, `hetionet_gene_id`, "
        "`hetionet_gene_symbol`, `mapping_status`.",
    ]
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
            "Ensembl -> Entrez uses the HGNC complete set "
            "(approved records only), not gene-symbol matching.",
            "Existing mapping artifacts are re-used when present and are never "
            "overwritten in place; corrections are written to a new filename and "
            "explained here.",
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
