#!/usr/bin/env python3
"""D5: context filter + gene-level support filter, with the thresholds fixed in one place.

Context stage (Pathway only):  q < Q_MAX  and  fold >= FOLD_MIN  and  SIZE_MIN <= size <= SIZE_MAX
                            and  seed_overlap >= MIN_OVERLAP
Gene stage:                 candidate gene (not a seed) kept if it belongs to >= MIN_SUPPORT kept contexts;
                            ranked by that support count.

Seed = DaG genes only (212), the decided primary seed (DuG/DdG edges come from expression meta-analysis).
The all-seed set (DaG|DuG|DdG, 698) is rerun with the same thresholds as a sensitivity check.
Biological Process is out of scope (decided 2026-10-06). q is Benjamini-Hochberg over every Pathway
(same convention as seed_expansion.py).

    python scripts/candidate_filter.py --graph-dir <handoff>/output

Writes results/seed_expansion/d5_filter_summary.tsv and d5_candidates_pathway.tsv.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

import seed_expansion as se

ROOT = Path(__file__).resolve().parent.parent
Q_MAX = 0.05
FOLD_MIN = 2.0
SIZE_MIN = 10
SIZE_MAX = {"pathway": 500}
MIN_OVERLAP = 3
MIN_SUPPORT = 2
PRIMARY_SEED = "DaG_only"
ROUTES = (("pathway", "GpPW", "Pathway"),)


def context_table(members, seed):
    universe = set().union(*members.values())
    N, n = len(universe), len(seed & universe)
    ids = list(members)
    ov = np.array([len(members[c] & seed) for c in ids])
    size = np.array([len(members[c]) for c in ids])
    p = np.where(ov > 0, hypergeom.sf(ov - 1, N, size, n), 1.0)
    return pd.DataFrame({"id": ids, "size": size, "seed_overlap": ov, "p": p,
                         "q": se.benjamini_hochberg(p), "fold": ov / (size * n / N)})


def kept_contexts(table, route):
    return table[(table["q"] < Q_MAX) & (table["fold"] >= FOLD_MIN) & (table["size"] >= SIZE_MIN)
                 & (table["size"] <= SIZE_MAX[route]) & (table["seed_overlap"] >= MIN_OVERLAP)]


def support(members, ids, seed):
    return Counter(g for c in ids for g in members[c] - seed)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=ROOT / "results/seed_expansion")
    args = ap.parse_args()

    nodes = pd.read_csv(args.graph_dir / "nodes.tsv", **se.READ)
    edges = pd.read_csv(args.graph_dir / "edges.tsv", **se.READ)
    kind, name = dict(zip(nodes["id"], nodes["kind"])), dict(zip(nodes["id"], nodes["name"]))
    k = se.KIDNEY_CANCER
    d = edges[((edges["source"] == k) | (edges["target"] == k)) & edges["metaedge"].isin(se.SEED_METAEDGES)]
    gene = d["source"].where(d["source"] != k, d["target"])
    seeds = {"DaG_only": set(gene[d["metaedge"] == "DaG"]), "all": set(gene)}

    rows = []
    for route, meta, ctx_kind in ROUTES:
        members = se.context_membership(edges, kind, meta, ctx_kind)
        kept_ids = {}
        for label, seed in seeds.items():
            kept = kept_contexts(context_table(members, seed), route)
            ids = list(kept["id"])
            kept_ids[label] = set(ids)
            sup = support(members, ids, seed)
            rows.append({"route": route, "seed": label, "n_seed": len(seed), "n_contexts": len(ids),
                         "role": "primary" if label == PRIMARY_SEED else "sensitivity",
                         "candidates_support>=1": len(sup),
                         f"candidates_support>={MIN_SUPPORT}": sum(v >= MIN_SUPPORT for v in sup.values()),
                         "candidates_support>=3": sum(v >= 3 for v in sup.values()),
                         "candidates_support>=5": sum(v >= 5 for v in sup.values())})
            if label == PRIMARY_SEED:
                out = pd.DataFrame([(g, name.get(g, ""), v) for g, v in sup.items() if v >= MIN_SUPPORT],
                                   columns=["gene_id", "gene_name", "n_supporting_contexts"])
                out.sort_values(["n_supporting_contexts", "gene_id"], ascending=[False, True]).to_csv(
                    args.out / f"d5_candidates_{route}.tsv", sep="\t", index=False)
        a, b = kept_ids[PRIMARY_SEED], kept_ids["all"]
        rows[-2]["jaccard_contexts_primary_vs_all_seed"] = round(len(a & b) / len(a | b), 3) if a | b else 0.0
    pd.DataFrame(rows).to_csv(args.out / "d5_filter_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
