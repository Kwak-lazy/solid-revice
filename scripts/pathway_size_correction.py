#!/usr/bin/env python3
"""Rank pathways by how over-represented kidney-cancer seed genes are, corrected for size.

A pathway's raw seed count mostly reflects how big the pathway is. This script
turns the count into an enrichment: how many seed genes the pathway holds versus
how many a pathway of that size would hold if seeds were scattered at random.

    python scripts/pathway_size_correction.py --graph-dir <handoff>/output

Writes results/seed_expansion/pathway_enrichment.tsv (one row per pathway) and
results/seed_expansion/size_correction_summary.tsv.

Per pathway (K genes, k of them seeds), with N = genes that have any pathway
annotation and n = seeds among them:
  expected       K * n / N
  fold           k / expected
  p              P(X >= k), hypergeometric
  q              Benjamini-Hochberg over all tested pathways
  matched_fold   k / mean overlap of random gene sets that copy the seeds'
                 pathway-count profile (controls for annotation bias, not only size)
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import hypergeom, spearmanr

ROOT = Path(__file__).resolve().parent.parent
KIDNEY_CANCER = "Disease::DOID:263"
READ = dict(sep="\t", dtype=str, keep_default_na=False)
N_PERM = 2000
SEED_RNG = 20261006
COUNT_BINS = [0, 1, 2, 3, 5, 8, 12, 20, 40, 10**9]   # right-closed pathway-count bins
KNOWN = {"HIF-1": "HIF-1-alpha transcription factor network", "VEGF": "Signaling by VEGF",
         "mTOR": "mTOR signaling pathway"}


def bh(p):
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    running = 1.0
    for rank, i in zip(range(m, 0, -1), order[::-1]):
        running = min(running, p[i] * m / rank)
        q[i] = running
    return q


def seeds_from_edges(edges, metaedges):
    d = edges[((edges["source"] == KIDNEY_CANCER) | (edges["target"] == KIDNEY_CANCER))
              & edges["metaedge"].isin(metaedges)]
    return set(d["source"].where(d["source"] != KIDNEY_CANCER, d["target"]))


def enrichment(members, seed, universe, total_tests):
    N, n = len(universe), len(seed & universe)
    rows = []
    for pid, genes in members.items():
        K, k = len(genes), len(genes & seed)
        exp = K * n / N
        rows.append((pid, K, k, exp, k / exp if exp else np.nan, hypergeom.sf(k - 1, N, K, n)))
    df = pd.DataFrame(rows, columns=["pathway_id", "size", "seed_overlap", "expected", "fold", "p"])
    # pathways with k = 0 have p = 1 and sit at the tail, so BH over all of them
    # only changes the m/rank factor; total_tests makes that explicit.
    p = df["p"].to_numpy()
    padded = np.concatenate([p, np.ones(max(total_tests - len(p), 0))])
    df["q"] = bh(padded)[: len(p)]
    return df


def matched_null(genes, gene_index, seed_idx, pw_count, matrix, rng):
    """Random gene sets with the same pathway-count profile as the seeds."""
    bins = np.digitize(pw_count, COUNT_BINS[1:-1], right=True)
    seed_mask = np.zeros(len(genes), bool)
    seed_mask[seed_idx] = True
    pools, need = {}, {}
    for b in np.unique(bins[seed_idx]):
        pools[b] = np.where((bins == b) & ~seed_mask)[0]
        need[b] = int((bins[seed_idx] == b).sum())
        assert len(pools[b]) >= need[b], f"bin {b}: pool {len(pools[b])} < seeds {need[b]}"
    rows, cols = [], []
    for i in range(N_PERM):
        for b, pool in pools.items():
            pick = rng.choice(pool, size=need[b], replace=False)
            rows.extend([i] * len(pick))
            cols.extend(pick.tolist())
    P = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N_PERM, len(genes)))
    overlap = (P @ matrix).toarray()
    return overlap.mean(0), overlap.std(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=ROOT / "results/seed_expansion")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    nodes = pd.read_csv(args.graph_dir / "nodes.tsv", **READ)
    edges = pd.read_csv(args.graph_dir / "edges.tsv", **READ)
    kind, name = dict(zip(nodes["id"], nodes["kind"])), dict(zip(nodes["id"], nodes["name"]))

    e = edges[edges["metaedge"] == "GpPW"]
    gene_is_src = e["source"].str.startswith("Gene::")
    pairs = pd.DataFrame({"g": e["source"].where(gene_is_src, e["target"]),
                          "c": e["target"].where(gene_is_src, e["source"])})
    pairs = pairs[pairs["c"].map(kind) == "Pathway"]
    members = pairs.groupby("c")["g"].apply(set).to_dict()
    universe = set(pairs["g"])
    n_pathways = len(members)

    seed_all = seeds_from_edges(edges, ("DaG", "DuG", "DdG"))
    seed_da = seeds_from_edges(edges, ("DaG",))

    df = enrichment(members, seed_all, universe, n_pathways)
    df_1230 = enrichment({k: v for k, v in members.items() if v & seed_all}, seed_all, universe,
                         sum(1 for v in members.values() if v & seed_all))
    df["pathway_name"] = df["pathway_id"].map(name)

    # annotation-matched null (controls for seeds being better annotated than average)
    genes = sorted(universe)
    gidx = {g: i for i, g in enumerate(genes)}
    pids = list(members)
    rows_, cols_ = zip(*[(gidx[g], j) for j, p in enumerate(pids) for g in members[p]])
    M = sparse.csr_matrix((np.ones(len(rows_)), (rows_, cols_)), shape=(len(genes), len(pids)))
    pw_count = np.asarray(M.sum(1)).ravel()
    seed_idx = np.array([gidx[g] for g in seed_all & universe])
    mean_null, sd_null = matched_null(genes, gidx, seed_idx, pw_count, M, np.random.default_rng(SEED_RNG))
    pos = {p: j for j, p in enumerate(pids)}
    df["matched_expected"] = df["pathway_id"].map(lambda p: mean_null[pos[p]])
    df["matched_fold"] = df["seed_overlap"] / df["matched_expected"].replace(0, np.nan)
    df["matched_z"] = (df["seed_overlap"] - df["matched_expected"]) / \
        df["pathway_id"].map(lambda p: sd_null[pos[p]]).replace(0, np.nan)

    df["rank_raw"] = df["seed_overlap"].rank(ascending=False, method="min").astype(int)
    df = df.sort_values(["q", "p", "fold"], ascending=[True, True, False]).reset_index(drop=True)
    df["rank_corrected"] = np.arange(1, len(df) + 1)
    cols = ["rank_corrected", "rank_raw", "pathway_id", "pathway_name", "size", "seed_overlap",
            "expected", "fold", "p", "q", "matched_expected", "matched_fold", "matched_z"]
    df[cols].to_csv(args.out / "pathway_enrichment.tsv", sep="\t", index=False, float_format="%.6g")

    # ---- summary ----
    sig = df["q"] < 0.05
    top_raw = df.nsmallest(50, "rank_raw")
    da_df = enrichment(members, seed_da, universe, n_pathways).sort_values(["q", "p"])
    top20, top20_da = set(df.head(20)["pathway_id"]), set(da_df.head(20)["pathway_id"])
    # Annotation-matched comparison is only fair among pathways big enough for a z-score to
    # mean something: with K = 3 and 2 seeds the null sd is tiny and z explodes.
    fair = df[df["size"].between(10, 500) & (df["seed_overlap"] >= 3)]
    fair_q20 = set(fair.head(20)["pathway_id"])          # df is already sorted by q
    fair_m20 = set(fair.nlargest(20, "matched_z")["pathway_id"])
    crit = sig & (df["fold"] >= 2) & df["size"].between(5, 500)

    def pool(ids):
        return len(set().union(*[members[x] for x in ids]) - seed_all)
    ids = list(members)
    sizes = np.array([len(members[p]) for p in ids])
    raw_k = np.array([len(members[p] & seed_all) for p in ids])
    S = [
        ("pathways_tested", n_pathways),
        ("universe_genes_N", len(universe)), ("seeds_in_universe_n", len(seed_all & universe)),
        ("spearman_size_vs_raw_seed_count", round(spearmanr(sizes, raw_k)[0], 3)),
        ("spearman_size_vs_fold", round(spearmanr(df["size"], df["fold"], nan_policy="omit")[0], 3)),
        ("spearman_size_vs_neglog10p", round(spearmanr(df["size"], -np.log10(df["p"]))[0], 3)),
        ("raw_top20_median_size", float(df.nsmallest(20, "rank_raw")["size"].median())),
        ("corrected_top20_median_size", float(df.head(20)["size"].median())),
        ("raw_top50_also_q<0.05", int((top_raw["q"] < 0.05).sum())),
        ("q<0.05_all_1822_tests", int(sig.sum())),
        ("q<0.05_only_pathways_with_overlap", int((df_1230["q"] < 0.05).sum())),
        ("q<0.05_and_fold>=2", int((sig & (df["fold"] >= 2)).sum())),
        ("q<0.05_and_fold>=2_and_size_5-500", int((sig & (df["fold"] >= 2) & df["size"].between(5, 500)).sum())),
        ("top20_overlap_with_DaG_only_seed", len(top20 & top20_da)),
        ("fair_set_size_10-500_overlap>=3", len(fair)),
        ("fair_top20_overlap_q_vs_annotation_matched_z", len(fair_q20 & fair_m20)),
        ("fair_spearman_neglog10p_vs_matched_z",
         round(spearmanr(-np.log10(fair["p"]), fair["matched_z"], nan_policy="omit")[0], 3)),
        ("top20_median_matched_expected_over_expected",
         round(float((df.head(20)["matched_expected"] / df.head(20)["expected"]).median()), 3)),
        ("candidate_genes_from_raw_top10_pathways", pool(df.nsmallest(10, "rank_raw")["pathway_id"])),
        ("candidate_genes_from_corrected_top10_pathways", pool(df.head(10)["pathway_id"])),
        ("candidate_genes_from_criteria_set", pool(df[crit]["pathway_id"])),
    ]
    for label, pat in KNOWN.items():
        r = df[df["pathway_name"] == pat]
        if len(r):
            S.append((f"rank_corrected_{label}", int(r["rank_corrected"].iloc[0])))
            S.append((f"rank_raw_{label}", int(r["rank_raw"].iloc[0])))
    pd.DataFrame(S, columns=["metric", "value"]).to_csv(args.out / "size_correction_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(S, columns=["metric", "value"]).to_string(index=False))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
