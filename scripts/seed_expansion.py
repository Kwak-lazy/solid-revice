#!/usr/bin/env python3
"""How far does the kidney-cancer seed reach, and how much must be filtered?

Starts at ``Disease::DOID:263`` and counts what each hop adds, then how the
candidate pool shrinks under the filters an experiment could apply. Everything
is read from the output of the collaborator's ``apply_kirc_to_hetionet_subgraph.py``
(the Hetionet base subgraph restricted to KIRC-expressed genes), so these are the
numbers the experiments will actually see.

    python scripts/seed_expansion.py --graph-dir <handoff>/output

Writes results/seed_expansion/*.tsv.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
KIDNEY_CANCER = "Disease::DOID:263"
SEED_METAEDGES = ("DaG", "DuG", "DdG")
SIZE_CAPS = (None, 500, 200, 100, 50)
MIN_OVERLAPS = (1, 2, 3, 5, 10)
# ccRCC genes from the TCGA-KIRC landmark study plus the hypoxia axis.
KNOWN_GENES = ("VHL", "PBRM1", "SETD2", "BAP1", "KDM5C", "MTOR", "PTEN", "TP53",
               "PIK3CA", "TSC1", "TSC2", "HIF1A", "EPAS1", "VEGFA", "CA9", "FLT1", "KDR")
READ = dict(sep="\t", dtype=str, keep_default_na=False)


def context_membership(edges, kind, metaedge, ctx_kind):
    """context_id -> set of Gene ids, for one Gene->Context metaedge."""
    e = edges[edges["metaedge"] == metaedge]
    gene_is_source = e["source"].str.startswith("Gene::")
    gene = e["source"].where(gene_is_source, e["target"])
    ctx = e["target"].where(gene_is_source, e["source"])
    pairs = pd.DataFrame({"g": gene, "c": ctx})
    pairs = pairs[pairs["c"].map(kind) == ctx_kind]
    return pairs.groupby("c")["g"].apply(set).to_dict()


def benjamini_hochberg(p):
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    running = 1.0
    for rank, i in zip(range(m, 0, -1), order[::-1]):
        running = min(running, p[i] * m / rank)
        q[i] = running
    return q


def enriched(members, seed, alpha=0.05):
    """Contexts over-represented among seed genes (hypergeometric, BH-FDR)."""
    universe = set().union(*members.values())
    N, n = len(universe), len(seed & universe)
    ids = [c for c in members if members[c] & seed]
    p = np.array([hypergeom.sf(len(members[c] & seed) - 1, N, len(members[c]), n) for c in ids])
    q = benjamini_hochberg(p)
    return {ids[i] for i in range(len(ids)) if q[i] < alpha}


def candidates(members, contexts, seed):
    return (set().union(*[members[c] for c in contexts]) - seed) if contexts else set()


def funnel(members, reached, seed):
    rows = []
    sizes = {c: len(members[c]) for c in reached}
    overlap = {c: len(members[c] & seed) for c in reached}
    for cap in SIZE_CAPS:
        for m in MIN_OVERLAPS:
            keep = [c for c in reached if overlap[c] >= m and (cap is None or sizes[c] <= cap)]
            rows.append({"size_cap": "none" if cap is None else cap, "min_seed_overlap": m,
                         "n_context": len(keep),
                         "n_candidate_genes": len(candidates(members, keep, seed))})
    sig = enriched(members, seed)
    rows.append({"size_cap": "none", "min_seed_overlap": "FDR<0.05", "n_context": len(sig),
                 "n_candidate_genes": len(candidates(members, sig, seed))})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph-dir", type=Path, required=True,
                    help="folder holding nodes.tsv and edges.tsv from the handoff script")
    ap.add_argument("--hetionet-nodes", type=Path, default=ROOT / "data/raw/hetionet-v1.0-nodes.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "results/seed_expansion")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    nodes = pd.read_csv(args.graph_dir / "nodes.tsv", **READ)
    edges = pd.read_csv(args.graph_dir / "edges.tsv", **READ)
    kind = dict(zip(nodes["id"], nodes["kind"]))
    n_graph_genes = int((nodes["kind"] == "Gene").sum())

    # hop 1: kidney cancer -> Gene
    d = edges[((edges["source"] == KIDNEY_CANCER) | (edges["target"] == KIDNEY_CANCER))
              & edges["metaedge"].isin(SEED_METAEDGES)]
    gene = d["source"].where(d["source"] != KIDNEY_CANCER, d["target"])
    by_type = {t: set(gene[d["metaedge"] == t]) for t in SEED_METAEDGES}
    seed = set(gene)
    daG_only_types = by_type["DaG"]
    expr_only = (by_type["DuG"] | by_type["DdG"]) - by_type["DaG"]

    pw = context_membership(edges, kind, "GpPW", "Pathway")
    bp = context_membership(edges, kind, "GpBP", "Biological Process")
    g2pw = defaultdict(set)
    g2bp = defaultdict(set)
    for c, gs in pw.items():
        for g in gs:
            g2pw[g].add(c)
    for c, gs in bp.items():
        for g in gs:
            g2bp[g].add(c)

    # hop 2: seed -> Pathway / BP / neighbouring Gene
    reached_pw = set().union(*[g2pw.get(g, set()) for g in seed])
    reached_bp = set().union(*[g2bp.get(g, set()) for g in seed])
    gg = edges[edges["metaedge"].isin(["GiG", "Gr>G"])]
    adj = defaultdict(set)
    for s, t in zip(gg["source"], gg["target"]):
        adj[s].add(t)
        adj[t].add(s)
    nb1 = set().union(*[adj[g] for g in seed]) - seed
    nb2 = set().union(*[adj[g] for g in nb1]) - seed - nb1

    # hop 3: candidate genes
    cand_pw, cand_bp = candidates(pw, reached_pw, seed), candidates(bp, reached_bp, seed)
    annotated_pw, annotated_bp = set(g2pw), set(g2bp)

    hop = [
        ("hop 0", "kidney cancer node", 1),
        ("hop 1", "seed genes (DaG | DuG | DdG)", len(seed)),
        ("hop 1", "  via DaG", len(by_type["DaG"])),
        ("hop 1", "  via DuG", len(by_type["DuG"])),
        ("hop 1", "  via DdG", len(by_type["DdG"])),
        ("hop 1", "  DaG only (no expression-derived edge)", len(by_type["DaG"] - by_type["DuG"] - by_type["DdG"])),
        ("hop 1", "  DuG/DdG only (no DaG)", len(expr_only)),
        ("hop 2", "pathways reached via GpPW", len(reached_pw)),
        ("hop 2", "  seed genes with a pathway", sum(g in g2pw for g in seed)),
        ("hop 2", "biological processes reached via GpBP", len(reached_bp)),
        ("hop 2", "  seed genes with a BP", sum(g in g2bp for g in seed)),
        ("hop 2", "neighbouring genes (GiG | Gr>G)", len(nb1)),
        ("hop 3", "genes one further step out (neighbours of neighbours)", len(nb2)),
        ("hop 3", "candidate genes via pathway", len(cand_pw)),
        ("hop 3", "candidate genes via BP", len(cand_bp)),
        ("hop 3", "candidate genes via any of the three routes", len(cand_pw | cand_bp | nb1)),
        ("ref", "genes in graph", n_graph_genes),
        ("ref", "genes with any pathway annotation", len(annotated_pw)),
        ("ref", "genes with any BP annotation", len(annotated_bp)),
    ]
    pd.DataFrame(hop, columns=["stage", "item", "count"]).to_csv(args.out / "hop_expansion.tsv", sep="\t", index=False)
    funnel(pw, reached_pw, seed).to_csv(args.out / "funnel_pathway.tsv", sep="\t", index=False)
    funnel(bp, reached_bp, seed).to_csv(args.out / "funnel_biological_process.tsv", sep="\t", index=False)

    # how much does the seed definition matter?
    rows = []
    for label, s in (("all (DaG+DuG+DdG)", seed), ("DaG only", by_type["DaG"]),
                     ("DuG/DdG only", expr_only)):
        sp, sb = enriched(pw, s), enriched(bp, s)
        rows.append({"seed_definition": label, "n_seed": len(s),
                     "enriched_pathways": len(sp), "candidate_genes_pathway": len(candidates(pw, sp, s)),
                     "enriched_bp": len(sb), "candidate_genes_bp": len(candidates(bp, sb, s))})
    pd.DataFrame(rows).to_csv(args.out / "seed_definition_comparison.tsv", sep="\t", index=False)

    # where do known ccRCC genes sit?
    het = pd.read_csv(args.hetionet_nodes, **READ)
    sym2id = dict(zip(het.loc[het["kind"] == "Gene", "name"], het.loc[het["kind"] == "Gene", "id"]))
    kinds = defaultdict(set)
    for g, m in zip(gene, d["metaedge"]):
        kinds[g].add(m)
    rows = []
    for sym in KNOWN_GENES:
        gid = sym2id.get(sym)
        rows.append({"gene": sym, "hetionet_gene_id": gid, "in_seed": gid in seed,
                     "seed_edge_types": ",".join(sorted(kinds.get(gid, ()))),
                     "n_pathways": len(g2pw.get(gid, ())), "n_bp": len(g2bp.get(gid, ()))})
    pd.DataFrame(rows).to_csv(args.out / "known_gene_check.tsv", sep="\t", index=False)

    print(pd.DataFrame(hop, columns=["stage", "item", "count"]).to_string(index=False))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
