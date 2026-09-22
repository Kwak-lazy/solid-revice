"""NetworkX subgraph extraction, carried over from the original notebook.

The original kept 4 node types (Gene, Pathway, Disease, Biological Process)
and 7 edge types (GpPW, GiG, Gr>G, DaG, DuG, DdG, GpBP) in one MultiDiGraph.

The only change: the Gene->Context layer is now selected by ``graph_context``,
so a run keeps ``GpPW`` **or** ``GpBP``, never both. Everything else - the
pandas pre-filter, the MultiDiGraph, isolate removal, the statistics - is the
original logic.
"""

from __future__ import annotations

import time
from collections import Counter

import pandas as pd

import config

from . import hetionet


def _require_networkx():
    try:
        import networkx as nx
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "networkx is required for subgraph extraction: pip install networkx"
        ) from exc
    return nx


# ---------------------------------------------------------------------------
# pandas pre-filter (original step 6-1)
# ---------------------------------------------------------------------------
def filter_edges(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    graph_context: str,
    kirc_gene_ids=None,
    verbose: bool = True,
):
    """Filter 2.25M edges down to the ones the subgraph needs.

    Filtering in pandas before touching NetworkX is what keeps this at a few
    seconds instead of minutes - the original notebook's key optimisation.

    When ``kirc_gene_ids`` is given, Gene nodes are additionally restricted to
    the KIRC-mapped ones (the ``Gene::Entrez`` join).
    """
    t0 = time.time()
    keep_kinds = config.keep_node_kinds(graph_context)
    keep_edges = config.keep_edge_types(graph_context)

    keep_nodes_df = nodes[nodes["kind"].isin(keep_kinds)]
    keep_node_ids = set(keep_nodes_df["id"])

    if kirc_gene_ids is not None:
        kirc_set = set(map(str, kirc_gene_ids))
        gene_ids = set(keep_nodes_df.loc[keep_nodes_df["kind"] == "Gene", "id"])
        keep_node_ids = (keep_node_ids - gene_ids) | (gene_ids & kirc_set)

    mask = (
        edges["metaedge"].isin(keep_edges)
        & edges["source"].isin(keep_node_ids)
        & edges["target"].isin(keep_node_ids)
    )
    filtered = edges[mask]

    if verbose:
        print(f"  keep node kinds ({len(keep_kinds)}): {sorted(keep_kinds)}")
        print(f"  keep edge types ({len(keep_edges)}): {sorted(keep_edges)}")
        print(f"  candidate nodes : {len(keep_node_ids):,}")
        print(f"  edges {len(edges):,} -> {len(filtered):,}  ({time.time() - t0:.1f}s)")

    return filtered, keep_node_ids


# ---------------------------------------------------------------------------
# MultiDiGraph (original step 6-2 / 6-3)
# ---------------------------------------------------------------------------
def build_graph(
    nodes: pd.DataFrame,
    filtered_edges: pd.DataFrame,
    keep_node_ids: set,
    drop_isolates: bool = True,
    verbose: bool = True,
):
    """Build the MultiDiGraph and drop isolated nodes.

    MultiDiGraph because the same node pair can carry several edge types
    (Gene A -> Gene B can be both ``GiG`` and ``Gr>G``) and because several
    edge types are directed (``DaG`` runs Disease -> Gene).
    """
    nx = _require_networkx()
    t0 = time.time()

    node_info = {
        row.id: {"name": row.name, "kind": row.kind}
        for row in nodes.itertuples(index=False)
    }

    G = nx.MultiDiGraph()
    # sorted(): keep_node_ids is a set, so iterating it directly makes node
    # insertion order - and therefore every exported row order - vary between
    # runs on identical input.
    for nid in sorted(keep_node_ids):
        info = node_info[nid]
        G.add_node(nid, name=info["name"], kind=info["kind"])
    for row in filtered_edges.itertuples(index=False):
        G.add_edge(row.source, row.target, metaedge=row.metaedge)

    if verbose:
        print(f"  graph built     : {G.number_of_nodes():,} nodes / "
              f"{G.number_of_edges():,} edges  ({time.time() - t0:.1f}s)")

    if drop_isolates:
        isolates = list(nx.isolates(G))
        if isolates and verbose:
            iso_types = Counter(G.nodes[n]["kind"] for n in isolates)
            print(f"  isolated nodes  : {len(isolates):,}  "
                  f"({', '.join(f'{k} {v:,}' for k, v in iso_types.most_common())})")
        G.remove_nodes_from(isolates)

    if verbose:
        print(f"  final subgraph  : {G.number_of_nodes():,} nodes / "
              f"{G.number_of_edges():,} edges")
    return G


# ---------------------------------------------------------------------------
# Statistics (original steps 7-9)
# ---------------------------------------------------------------------------
def node_kind_counts(G) -> pd.DataFrame:
    counts = Counter(d["kind"] for _, d in G.nodes(data=True))
    return pd.DataFrame(
        sorted(counts.items(), key=lambda kv: -kv[1]), columns=["kind", "n_nodes"]
    )


def edge_type_counts(G) -> pd.DataFrame:
    counts = Counter(d["metaedge"] for _, _, d in G.edges(data=True))
    out = pd.DataFrame(
        sorted(counts.items(), key=lambda kv: -kv[1]), columns=["metaedge", "n_edges"]
    )
    out["meaning"] = out["metaedge"].map(config.METAEDGE_NAMES).fillna("")
    return out


def context_gene_counts(G, graph_context: str) -> pd.DataFrame:
    """Genes per context node (original step 8, generalised over the context).

    The context metaedge runs Gene -> Context, so the genes of a context node
    are the sources of its incoming edges of that metaedge.
    """
    cfg = config.CONTEXT_CONFIG[graph_context]
    metaedge = cfg["metaedge"]
    context_nodes = [n for n, d in G.nodes(data=True) if d["kind"] == cfg["node_kind"]]

    rows = []
    for node in context_nodes:
        n_genes = sum(
            1 for _, _, d in G.in_edges(node, data=True)
            if d.get("metaedge") == metaedge
        )
        rows.append(
            {
                "context_id": node,
                "context_name": G.nodes[node]["name"],
                "context_type": cfg["node_kind"],
                "n_genes": n_genes,
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values("n_genes", ascending=False).reset_index(drop=True)


def describe_context_sizes(counts: pd.DataFrame, graph_context: str, top: int = 15):
    sizes = counts["n_genes"]
    print(f"  {graph_context} nodes: {len(counts):,}")
    print(f"    genes per node  min {sizes.min()} / max {sizes.max()} / "
          f"median {int(sizes.median())} / mean {sizes.mean():.1f}")
    print(f"\n  top {top} by gene count:")
    for i, row in enumerate(counts.head(top).itertuples(index=False), 1):
        print(f"    {i:<3d} {row.n_genes:>6,}  {row.context_id:<35s} {row.context_name}")


def kidney_disease_report(G) -> pd.DataFrame:
    """Kidney-related Disease nodes and their gene links (original step 9)."""
    rows = []
    for node, data in G.nodes(data=True):
        if data["kind"] != "Disease":
            continue
        name = data["name"].lower()
        if not any(kw in name for kw in config.KIDNEY_DISEASE_KEYWORDS):
            continue
        by_type = Counter(d["metaedge"] for _, _, d in G.out_edges(node, data=True))
        rows.append(
            {
                "disease_id": node,
                "disease_name": data["name"],
                "n_genes_total": sum(by_type.values()),
                **{f"n_{k}": v for k, v in by_type.items()},
            }
        )
    out = pd.DataFrame(rows).fillna(0)
    if len(out):
        out = out.sort_values("n_genes_total", ascending=False).reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Export (original step 11)
# ---------------------------------------------------------------------------
def to_frames(G):
    """Graph -> (nodes, edges) dataframes, in the original's column layout.

    Both frames are sorted so a re-run on identical input produces byte-identical
    files.
    """
    nodes_df = pd.DataFrame(
        [{"id": n, "name": d["name"], "kind": d["kind"]} for n, d in G.nodes(data=True)]
    ).sort_values(["kind", "id"]).reset_index(drop=True)
    edges_df = pd.DataFrame(
        [
            {"source": u, "metaedge": d["metaedge"], "target": v}
            for u, v, d in G.edges(data=True)
        ]
    ).sort_values(["metaedge", "source", "target"]).reset_index(drop=True)
    return nodes_df, edges_df


def extract(
    nodes: pd.DataFrame = None,
    edges: pd.DataFrame = None,
    graph_context: str = None,
    kirc_gene_ids=None,
    verbose: bool = True,
):
    """Full extraction for one context: filter -> graph -> isolates -> frames."""
    graph_context = config.GRAPH_CONTEXT if graph_context is None else graph_context
    nodes = hetionet.load_nodes() if nodes is None else nodes
    edges = hetionet.load_edges() if edges is None else edges

    filtered, keep_ids = filter_edges(
        nodes, edges, graph_context, kirc_gene_ids=kirc_gene_ids, verbose=verbose
    )
    G = build_graph(nodes, filtered, keep_ids, verbose=verbose)
    nodes_df, edges_df = to_frames(G)
    return {"graph": G, "nodes": nodes_df, "edges": edges_df}
