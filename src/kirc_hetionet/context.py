"""Gene -> Context subgraph extraction, driven by a single config switch.

Pathway and Biological Process are NOT separate pipelines.  They are the same
code path parameterised by ``graph_context``:

    pathway             Gene --GpPW--> Pathway
    biological_process  Gene --GpBP--> Biological Process

Both produce the identical result schema, so downstream stages never need to
know which context they are looking at.
"""

from __future__ import annotations

import pandas as pd

import config

from . import hetionet, subgraph as subgraph_mod

# The shared result schema for every context (config section 8 of the spec).
GENE_CONTEXT_COLUMNS = [
    "gene_id",
    "gene_symbol",
    "context_id",
    "context_name",
    "edge_type",
    "context_type",
]


# ---------------------------------------------------------------------------
# Common context accessors
# ---------------------------------------------------------------------------
def get_context_nodes(nodes: pd.DataFrame, graph_context: str) -> pd.DataFrame:
    config_ = config.CONTEXT_CONFIG[graph_context]
    return nodes[nodes["kind"] == config_["node_kind"]].copy()


def get_context_edges(edges: pd.DataFrame, graph_context: str) -> pd.DataFrame:
    config_ = config.CONTEXT_CONFIG[graph_context]
    return edges[edges["metaedge"] == config_["metaedge"]].copy()


# ---------------------------------------------------------------------------
# Subgraph assembly
# ---------------------------------------------------------------------------
def build_gene_context_frame(
    context_edges: pd.DataFrame,
    context_nodes: pd.DataFrame,
    gene_nodes: pd.DataFrame,
    graph_context: str,
) -> pd.DataFrame:
    """Join edges to node names and return the shared 6-column schema."""
    cfg = config.CONTEXT_CONFIG[graph_context]

    frame = context_edges.rename(columns={"source": "gene_id", "target": "context_id"})
    frame = frame.merge(
        gene_nodes.rename(
            columns={
                "hetionet_gene_id": "gene_id",
                "hetionet_gene_symbol": "gene_symbol",
            }
        )[["gene_id", "gene_symbol"]],
        on="gene_id",
        how="left",
    )
    frame = frame.merge(
        context_nodes.rename(columns={"id": "context_id", "name": "context_name"})[
            ["context_id", "context_name"]
        ],
        on="context_id",
        how="left",
    )
    frame["edge_type"] = cfg["metaedge"]
    frame["context_type"] = cfg["node_kind"]
    return frame[GENE_CONTEXT_COLUMNS].reset_index(drop=True)


# ---------------------------------------------------------------------------
# One experiment
# ---------------------------------------------------------------------------
def run_context_experiment(
    graph_context: str = None,
    kirc_gene_ids=None,
    nodes: pd.DataFrame = None,
    edges: pd.DataFrame = None,
    write: bool = True,
    verbose: bool = True,
    build_subgraph: bool = True,
) -> dict:
    """Run the full Gene -> Context extraction for one context.

    ``kirc_gene_ids`` is an iterable of Hetionet ``Gene::<entrez>`` ids - the
    graph join key.  When it is ``None`` the subgraph is the unrestricted
    Gene -> Context graph and the KIRC counters are reported as ``None``.

    Returns a dict with the three result frames and a ``summary`` of counts.
    """
    graph_context = config.GRAPH_CONTEXT if graph_context is None else graph_context
    if graph_context not in config.CONTEXT_CONFIG:
        raise KeyError(
            f"Unknown graph_context {graph_context!r}; "
            f"expected one of {sorted(config.CONTEXT_CONFIG)}"
        )
    cfg = config.CONTEXT_CONFIG[graph_context]

    nodes = hetionet.load_nodes() if nodes is None else nodes
    edges = hetionet.load_edges() if edges is None else edges
    gene_nodes = hetionet.gene_nodes(nodes)

    if verbose:
        print(f"\n{'=' * 64}")
        print(f"context: {graph_context}  "
              f"(node kind = {cfg['node_kind']!r}, metaedge = {cfg['metaedge']!r})")
        print("=" * 64)

    context_nodes = get_context_nodes(nodes, graph_context)
    context_edges = get_context_edges(edges, graph_context)

    # Only the selected metaedge may appear - guards against a config slip.
    stray = set(context_edges["metaedge"]) - {cfg["metaedge"]}
    assert not stray, f"unexpected metaedges for {graph_context}: {stray}"

    full = build_gene_context_frame(context_edges, context_nodes, gene_nodes, graph_context)

    # ---- restrict to KIRC-mapped genes (the Gene::Entrez join) -------------
    if kirc_gene_ids is not None:
        kirc_set = set(map(str, kirc_gene_ids))
        sub = full[full["gene_id"].isin(kirc_set)].copy()
        n_kirc_mapped = sub["gene_id"].nunique()
    else:
        kirc_set = None
        sub = full.copy()
        n_kirc_mapped = None

    sub_context_ids = set(sub["context_id"])
    sub_gene_ids = set(sub["gene_id"])

    # ---- context node table with subgraph membership ----------------------
    gene_counts_all = full.groupby("context_id")["gene_id"].nunique()
    gene_counts_kirc = sub.groupby("context_id")["gene_id"].nunique()
    context_out = context_nodes.rename(
        columns={"id": "context_id", "name": "context_name", "kind": "context_type"}
    )[["context_id", "context_name", "context_type"]].copy()
    context_out["n_genes_total"] = (
        context_out["context_id"].map(gene_counts_all).fillna(0).astype(int)
    )
    context_out["n_genes_kirc"] = (
        context_out["context_id"].map(gene_counts_kirc).fillna(0).astype(int)
    )
    context_out["in_subgraph"] = context_out["context_id"].isin(sub_context_ids)
    context_out = context_out.sort_values(
        ["n_genes_kirc", "context_id"], ascending=[False, True]
    ).reset_index(drop=True)

    # ---- NetworkX subgraph (the original notebook's extraction) -----------
    # Node kinds: Gene, Disease + this context's kind.
    # Edge types: GiG, Gr>G, DaG, DuG, DdG + this context's metaedge.
    graph = None
    subgraph_edges_df = None
    if build_subgraph:
        if verbose:
            print("  --- subgraph extraction ---")
        extracted = subgraph_mod.extract(
            nodes, edges, graph_context, kirc_gene_ids=kirc_gene_ids, verbose=verbose
        )
        graph = extracted["graph"]
        subgraph_nodes = extracted["nodes"].rename(
            columns={"id": "node_id", "name": "node_name", "kind": "node_kind"}
        )
        subgraph_edges_df = extracted["edges"]
    else:
        gene_part = (
            gene_nodes[gene_nodes["hetionet_gene_id"].isin(sub_gene_ids)]
            .rename(
                columns={
                    "hetionet_gene_id": "node_id",
                    "hetionet_gene_symbol": "node_name",
                }
            )[["node_id", "node_name"]]
            .assign(node_kind="Gene")
        )
        context_part = context_nodes.rename(
            columns={"id": "node_id", "name": "node_name", "kind": "node_kind"}
        )
        context_part = context_part[context_part["node_id"].isin(sub_context_ids)][
            ["node_id", "node_name", "node_kind"]
        ]
        subgraph_nodes = pd.concat([gene_part, context_part], ignore_index=True)

    summary = {
        "graph_context": graph_context,
        "node_kind": cfg["node_kind"],
        "metaedge": cfg["metaedge"],
        "context_nodes": int(len(context_nodes)),
        "gene_context_edges": int(len(full)),
        "kirc_mapped_genes": None if n_kirc_mapped is None else int(n_kirc_mapped),
        "subgraph_nodes": int(len(subgraph_nodes)),
        "subgraph_edges": int(
            len(subgraph_edges_df) if subgraph_edges_df is not None else len(sub)
        ),
        "subgraph_context_nodes": int(len(sub_context_ids)),
        "subgraph_gene_nodes": int(len(sub_gene_ids)),
        "context_nodes_without_kirc_gene": int(
            len(context_nodes) - len(sub_context_ids)
        ),
    }
    if graph is not None:
        kinds = subgraph_mod.node_kind_counts(graph).set_index("kind")["n_nodes"]
        for kind, value in kinds.items():
            summary[f"subgraph_{kind.lower().replace(' ', '_')}_nodes"] = int(value)

    if verbose:
        width = max(len(k) for k in summary)
        for key, value in summary.items():
            shown = "n/a" if value is None else (
                f"{value:,}" if isinstance(value, int) else value
            )
            print(f"  {key:<{width}} : {shown}")

    outputs = {}
    if write:
        out_dir = config.results_dir(graph_context)
        out_dir.mkdir(parents=True, exist_ok=True)
        targets = {
            "context_nodes": (out_dir / "context_nodes.tsv", context_out),
            "gene_context_edges": (out_dir / "gene_context_edges.tsv", sub),
            "subgraph_nodes": (out_dir / "subgraph_nodes.tsv", subgraph_nodes),
        }
        if subgraph_edges_df is not None:
            targets["subgraph_edges"] = (
                out_dir / "subgraph_edges.tsv",
                subgraph_edges_df,
            )
        for key, (path, frame) in targets.items():
            frame.to_csv(path, sep="\t", index=False)
            outputs[key] = path
            if verbose:
                print(f"  [write] {path}  ({len(frame):,} rows)")

    return {
        "graph_context": graph_context,
        "config": cfg,
        "context_nodes": context_out,
        "gene_context_edges": sub,
        "gene_context_edges_all": full,
        "subgraph_nodes": subgraph_nodes,
        "subgraph_edges": subgraph_edges_df,
        "graph": graph,
        "summary": summary,
        "outputs": outputs,
    }


# ---------------------------------------------------------------------------
# Both experiments + comparison
# ---------------------------------------------------------------------------
def run_all_context_experiments(
    kirc_gene_ids=None,
    nodes: pd.DataFrame = None,
    edges: pd.DataFrame = None,
    write: bool = True,
    verbose: bool = True,
    build_subgraph: bool = True,
) -> dict:
    """Run every context in ``config.ALL_GRAPH_CONTEXTS`` and compare them."""
    nodes = hetionet.load_nodes() if nodes is None else nodes
    edges = hetionet.load_edges() if edges is None else edges

    results = {}
    for graph_context in config.ALL_GRAPH_CONTEXTS:
        results[graph_context] = run_context_experiment(
            graph_context,
            kirc_gene_ids=kirc_gene_ids,
            nodes=nodes,
            edges=edges,
            write=write,
            verbose=verbose,
            build_subgraph=build_subgraph,
        )

    comparison = build_comparison(results, write=write)
    if verbose:
        print(f"\n{'=' * 64}\ncontext comparison\n{'=' * 64}")
        print(format_comparison(comparison))

    results["comparison"] = comparison
    return results


COMPARISON_ROWS = [
    ("Context nodes", "context_nodes"),
    ("Gene-context edges", "gene_context_edges"),
    ("KIRC mapped genes", "kirc_mapped_genes"),
    ("Subgraph nodes", "subgraph_nodes"),
    ("Subgraph edges", "subgraph_edges"),
]


def build_comparison(results: dict, write: bool = True) -> pd.DataFrame:
    contexts = [c for c in config.ALL_GRAPH_CONTEXTS if c in results]
    rows = []
    for label, key in COMPARISON_ROWS:
        row = {"metric": label}
        for ctx in contexts:
            row[ctx] = results[ctx]["summary"][key]
        rows.append(row)
    comparison = pd.DataFrame(rows)

    if write:
        config.COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
        path = config.COMPARISON_DIR / "context_comparison.tsv"
        comparison.to_csv(path, sep="\t", index=False)
        print(f"  [write] {path}")

        detail = (
            pd.DataFrame([results[ctx]["summary"] for ctx in contexts])
            .set_index("graph_context")
            .T.reset_index(names="metric")
        )
        detail_path = config.COMPARISON_DIR / "context_comparison_detail.tsv"
        detail.to_csv(detail_path, sep="\t", index=False)
        print(f"  [write] {detail_path}")

    return comparison


def format_comparison(comparison: pd.DataFrame) -> str:
    """Render the comparison the way the spec shows it."""
    headers = {
        "pathway": "Pathway",
        "biological_process": "Biological Process",
    }
    cols = [c for c in comparison.columns if c != "metric"]
    label_w = max(len(str(v)) for v in comparison["metric"])
    label_w = max(label_w, len("Metric"))
    col_w = {c: max(len(headers.get(c, c)), 12) for c in cols}

    def fmt(value):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return "n/a"
        return f"{int(value):,}"

    lines = [
        " " * label_w + "  " + "  ".join(headers.get(c, c).rjust(col_w[c]) for c in cols)
    ]
    lines.append("-" * len(lines[0]))
    for _, row in comparison.iterrows():
        cells = "  ".join(fmt(row[c]).rjust(col_w[c]) for c in cols)
        lines.append(f"{str(row['metric']).ljust(label_w)}  {cells}")
    return "\n".join(lines)
