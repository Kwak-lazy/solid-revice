# KIRC-Hetionet Project

## Current Status
Stage: Gene ID + Graph Context Pipeline

Status: In Progress

Last updated: 2026-09-22

Completed:
- Gene ID standardization
- Pathway / Biological Process switch
- Automated dual-context pipeline
- Automatic context comparison
- README progress logging
- Google Drive backup structure

Next:
- Candidate selection
- Variance-based selection

## 1. Project Goal
Build a reproducible KIRC (kidney renal clear cell carcinoma) x Hetionet
pipeline in which every gene is identified by a single standardized key, and in
which the biological context layer - **Pathway** or **Biological Process** - is
a configuration value rather than a separate code path.

This stage delivers only the identifier and graph-context foundation.
Variance scoring, candidate ranking, ML and GNN stages come later
(see section 10).

## 2. Current Pipeline
```text
 0. Environment + Google Drive mount + project paths
 1. Hetionet download           (nodes.tsv / edges.sif.gz, checksum-verified)
 2. Load + structure            (columns, row counts)
 3. Node type analysis          (11 kinds)
 4. Identifier analysis         (Gene::Entrez, Pathway PC7/WP, Disease DOID)
 5. Edge type analysis          (24 metaedges)
 6. TCGA-KIRC download          (star_tpm + clinical, Xena -> S3 fallback)
 7. Gene ID standardization     (Ensembl -> HGNC -> Entrez -> Gene::Entrez)
 8. Mapping validation
 9. Context configuration       (GRAPH_CONTEXT / CONTEXT_CONFIG)
10. Single-context run          (pandas filter -> MultiDiGraph -> drop isolates)
11. Subgraph analysis           (node/edge counts, gene-per-context, kidney Disease)
12. Dual-context run + comparison
13. Result files
14. README update
15. Google Drive backup
16. Summary
```

This notebook is a revision of the original `Download_and_subgraph.ipynb`, not a
replacement. Preserved unchanged: the download strategy, the pandas pre-filter,
the `nx.MultiDiGraph` construction, isolate removal, the 24-metaedge glossary,
the Pathway PC7/WP prefix analysis, the Disease DOID analysis, the kidney-disease
report, and the TCGA barcode sample-type breakdown.

Changed, and only where the new requirements demanded it:

| Item | Original | Now |
|---|---|---|
| Gene ID matching | deferred to a separate `step4` | done in this notebook (step 7) |
| Context layer | `GpPW` **and** `GpBP` kept together | one selected by `GRAPH_CONTEXT` |
| Kept node kinds | 4 (Gene, Pathway, Disease, BP) | 3 (Gene, Disease + context kind) |
| Kept edge types | 7 (fixed) | 6 (5 fixed + the context metaedge) |
| Gene nodes | all Hetionet genes | restricted to KIRC-mapped `Gene::Entrez` |
| Output | one `data/subgraph_*.tsv` | `results/<context>/` per context |
| Comparison | none | automatic, into `results/comparison/` |
| Record keeping | none | README auto-updated + Drive backup |

Layout:

```text
project/
├── notebooks/Download_and_subgraph.ipynb   # the 16 steps, interactively
├── config.py                               # paths, GRAPH_CONTEXT, CONTEXT_CONFIG,
│                                           #   KEEP node/edge sets, METAEDGE_NAMES
├── src/kirc_hetionet/
│   ├── hetionet.py     # download / load / node kinds / metaedges
│   ├── kirc.py         # Xena download, expression, sample types, clinical
│   ├── gene_mapping.py # Ensembl -> HGNC -> Entrez -> Gene::Entrez + validation
│   ├── subgraph.py     # pandas filter -> MultiDiGraph -> isolates -> stats
│   ├── context.py      # get_context_nodes / get_context_edges / experiments
│   ├── readme_log.py   # this file, kept up to date automatically
│   ├── report.py       # run -> README sections
│   └── drive.py        # mount / sync / verify
├── scripts/run_pipeline.py                 # the same steps, headless
├── data/{raw,processed,external}/
└── results/{pathway,biological_process,comparison}/
```

Run it:

```python
results = run_all_context_experiments()   # both contexts + comparison
```

```bash
python scripts/run_pipeline.py            # the whole flow, headless
```

Repository / branch: `Kwak-lazy/solid-revice`, branch
`claude/determined-lovelace-mj90py`.

## 3. Data
### Division of labour

Gene-ID standardization belongs to the collaborator's pipeline; graph
construction belongs to this one. The handoff is a fixed set of files and a
fixed key.

```text
collaborator                                    this repo
─────────────────────────────────────────────   ──────────────────────────────
1  build_kirc_feature_guide.py
     all_features.tsv (60,660)
2  build_kirc_gene_mapping.py
     kirc_gene_mapping_all.tsv  (hetionet cols blank)
3-4  build_kirc_hetionet_mapping.py
     kirc_hetionet_gene_nodes.tsv  ───────────▶  gene_mapping.load_collaborator_mapping()
     hetionet_genes_absent_from_kirc.tsv ─────▶  (graph-only nodes)
4-2  verify_auxiliary_candidates_ncbi.py
     ..._verified.tsv
5  build_kirc_standardized_outputs.py
     kirc_expression_standardized.tsv.gz
     kirc_sample_metadata.tsv
     kirc_mapping_exceptions.tsv
     kirc_gene_standardization_summary.md ────▶  (expected counts, verified on load)
                                                 subgraph.py -> context.py
```

Join key across the boundary: **`hetionet_gene_id`** (`Gene::<entrez>`), never a
symbol.

### Handoff files (`data/external/collab/`, read-only)

| File | Rows | Used for |
|---|---:|---|
| `kirc_hetionet_gene_nodes.tsv` | 19,425 | the gene universe; `all_samples_zero` splits it |
| `hetionet_genes_absent_from_kirc.tsv` | 1,520 | nodes kept in the graph but carrying no expression |
| `kirc_gene_mapping_all.tsv` | 60,660 | traceability from any Hetionet node back to its KIRC feature |
| `kirc_gene_standardization_summary.md` | — | rules, counts and source SHA-256s, asserted on load |

Produced by the collaborator but outside this handoff, so not present here:
`kirc_expression_standardized.tsv.gz`, `kirc_sample_metadata.tsv`,
`kirc_mapping_exceptions.tsv`, `kirc_duplicate_gene_review.tsv`,
`kirc_auxiliary_name_candidates*.tsv`, `all_features.tsv`,
`TCGA-KIRC.sample_labels.tsv`. The first two are needed for anything that reads
expression values; the rest are their internal audit trail.

### Shared sources - byte-identical where it matters

Both pipelines download Hetionet and TCGA-KIRC independently, so they were
checked against the SHA-256s recorded in the collaborator's summary:

| Source | Match |
|---|---|
| `hetionet-v1.0-nodes.tsv` | identical |
| `hetionet-v1.0-edges.sif.gz` | identical |
| `TCGA-KIRC.star_tpm.tsv.gz` | identical |
| `hgnc_complete_set.txt` | **differs** - theirs 2026-09-14, ours fetched later |

The HGNC difference is inert: HGNC feeds only the retired in-house mapping path,
which no longer runs. It is not an input to anything this pipeline now produces.

### Reading convention

The collaborator's scripts read every TSV with `keep_default_na=False`; this
pipeline now does the same. Without it pandas turns empty cells - and any
literal `NA`, `None` or `null` - into NaN. Hetionet's own files contain no such
strings, but `kirc_hetionet_gene_nodes.tsv` has one row with an empty
`current_hgnc_symbol` (`RSC1A1`, linked via NCBI cross-reference rather than
HGNC), and keeping both sides on the same convention removes the class of bug
rather than the one instance.

### Downloaded sources

| Dataset | Source | Local path |
|---|---|---|
| Hetionet v1.0 nodes (47,031) | `raw.githubusercontent.com/hetio/hetionet` | `data/raw/hetionet-v1.0-nodes.tsv` |
| Hetionet v1.0 edges (2,250,197) | git-lfs **media** endpoint, sha256-verified | `data/raw/hetionet-v1.0-edges.sif.gz` |
| TCGA-KIRC expression, `log2(TPM+1)` | UCSC Xena GDC hub (S3 fallback) | `data/raw/TCGA-KIRC.star_tpm.tsv.gz` |
| TCGA-KIRC clinical / survival | UCSC Xena GDC hub (S3 fallback) | `data/raw/TCGA-KIRC.{clinical,survival}.tsv.gz` |
| HGNC complete set | `storage.googleapis.com/public-download-files/hgnc` | `data/raw/hgnc_complete_set.txt` |

Measured: Hetionet Gene 20,945 · Biological Process 11,381 · Pathway 1,822 ·
Disease 137 (+7 further kinds). TCGA-KIRC 60,660 features x 610 samples;
537 Primary Tumor, 72 Solid Tissue Normal, 1 Additional New Primary.

`data/raw/` is git-ignored and re-downloaded on demand.
`results/*/gene_context_edges.tsv` and `results/*/subgraph_edges.tsv` are
git-ignored for size.

## 4. Gene ID Standardization
Join key: `hetionet_gene_id` (`Gene::<entrez>`). Gene symbols are annotation only and are never used to join.

**The gene mapping is not produced by this pipeline.** It is the collaborator's standardization run (2026-09-15), used as-is and never regenerated or overwritten.

```text
KIRC Ensembl ID (GENCODE v36, versioned)
   |  1st: GENCODE gene-line hgnc_id      38,891 features
   |  2nd: HGNC ensembl_gene_id lookup     3,177 features
   v
current HGNC symbol + entrez_id
   |  "Gene::" + entrez, only if the node exists in Hetionet
   v
Hetionet Gene::Entrez
```

Routing through **`hgnc_id`** rather than Ensembl ID is the whole point. HGNC IDs are curated and stable; Ensembl gene IDs get reassigned between builds. TCGA is frozen at GENCODE v36 while HGNC tracks the current release, so an Ensembl-to-Ensembl join silently drops every gene whose ID moved - SOD2 (`ENSG00000112096` -> `ENSG00000291237`) among them.

### Two gene universes - not interchangeable

| Universe | Genes | Use for |
|---|---:|---|
| `graph` | 19,425 | graph-only work (topology, DWPC, metapaths) |
| `expression` | 19,297 | anything needing expression values |

The 128-gene difference is nodes whose KIRC row is zero across all 610 samples. Selected via `--gene-universe` / `config.GENE_UNIVERSE`.

### Source files (read-only, `data/external/collab/`)

| File | Content |
|---|---|
| `kirc_gene_mapping_all.tsv` | 60,660 features x 36 columns |
| `kirc_hetionet_gene_nodes.tsv` | 19,425 connected Hetionet nodes |
| `hetionet_genes_absent_from_kirc.tsv` | 1,520 nodes with no KIRC feature |
| `kirc_gene_standardization_summary.md` | rules, counts, SHA-256 of every input |

Verified on load: feature / graph-node / expression-node / absent counts, plus disjointness and union against Hetionet's 20,945 Gene nodes. A mismatch raises rather than proceeding.

| Counter | Value |
|---|---:|
| features | 60,660 |
| graph_nodes | 19,425 |
| expression_nodes | 19,297 |
| absent_nodes | 1,520 |
| all_samples_zero | 128 |

## 5. Graph Context
A single switch selects the context; there is no second pipeline.

```python
GRAPH_CONTEXT = 'pathway'   # or 'biological_process'
```

| Context | node_kind | metaedge |
|---|---|---|
| `pathway` | `Pathway` | `GpPW` |
| `biological_process` | `Biological Process` | `GpBP` |

```text
Gene --GpPW--> Pathway
Gene --GpBP--> Biological Process
```

Both contexts emit the same schema: `gene_id`, `gene_symbol`, `context_id`, `context_name`, `edge_type`, `context_type`.

## 6. Experiment Progress
- **Note (2026-10-06): D3 and D5 decided - DaG-only seed, Pathway route only.**
  - Where: `scripts/candidate_filter.py`, `results/seed_expansion/d5_filter_summary.tsv`,
    `d5_candidates_pathway.tsv`.
  - Decision (user): seed = `DaG` genes (212); Biological Process is dropped from the experiments.
    The all-seed set (698) stays only as a sensitivity run. `GpBP` code paths in the pipeline are
    untouched; only the candidate filter no longer produces a BP list (the earlier BP candidate
    file was removed).
  - Thresholds: q < 0.05 (BH over all 1,822 Pathways), fold >= 2, size 10-500, seed overlap >= 3;
    gene stage keeps genes supported by >= 2 kept Pathways, ranked by that support count.
  - Measured (DaG seed): 221 Pathways -> 3,066 candidate genes -> 1,651 with support >= 2
    (>= 3: 1,044; >= 5: 609). All-seed run: 163 Pathways -> 2,741 -> 1,332. Kept-Pathway sets of the
    two seeds have Jaccard 0.488.
  - Check: 100 random seeds of 698 annotated genes gave on average 0.01 Pathways with q < 0.05
    (uniform random seeds only; not annotation-matched).
  - Caveat: top-supported genes are signalling hubs (MAP2K1 106, PIK3R1 101, HRAS 97, AKT1 96,
    GRB2 91 Pathways), so support count alone mixes cancer relevance with hub-ness.
  - Still open: D7 (evaluation) and a gene-level score; BP rows in earlier tables are historical.

- **Note (2026-10-06): size-correction run repeated for the DaG seed (Pathway).**
  - Where: `scripts/pathway_size_correction.py --seed DaG` -> `results/seed_expansion/pathway_enrichment_DaG.tsv`,
    `size_correction_summary_DaG.tsv`. The all-seed files keep their old names; the all-seed
    enrichment table was re-run and is byte-identical, the summary gained two rows
    (D5 filter set) and one renamed row (`top20_overlap_with_DaG_seed`).
  - Measured (DaG, 165 of 212 seeds annotated, 1,822 Pathways): Spearman(size, raw seed count) 0.629;
    raw top-20 median size 275 vs corrected 70; q<0.05 for 236; D5 filter set 221 (3,066 candidate genes);
    annotation-matched expected over uniform expected 1.87 (top 20); Spearman(-log10 p, matched z) 0.879.
  - Ranks (corrected / raw): VEGF 14 / 19, mTOR 24 / 60, HIF-1-alpha TF network 71 / 81,
    HIF-2-alpha TF network 13 / 69.
  - Resolved: yes (documented on the Notion decisions page, section 8).

### 2026-09-21

- [x] Gene ID standardization implemented (Ensembl -> HGNC -> Entrez -> `Gene::Entrez`)
- [x] Hetionet Gene mapping validation implemented and run
- [x] Pathway / Biological Process configuration (`GRAPH_CONTEXT` + `CONTEXT_CONFIG`)
- [x] Common pipeline implementation (`get_context_nodes` / `get_context_edges` / `run_context_experiment`)
- [x] Pathway experiment executed
- [x] Biological Process experiment executed
- [x] Automatic context comparison written to `results/comparison/context_comparison.tsv`

### 2026-09-22

- [x] Gene ID standardization implemented (Ensembl -> HGNC -> Entrez -> `Gene::Entrez`)
- [x] Hetionet Gene mapping validation implemented and run
- [x] Pathway / Biological Process configuration (`GRAPH_CONTEXT` + `CONTEXT_CONFIG`)
- [x] Common pipeline implementation (`get_context_nodes` / `get_context_edges` / `run_context_experiment`)
- [x] Pathway experiment executed
- [x] Biological Process experiment executed
- [x] Automatic context comparison written to `results/comparison/context_comparison.tsv`

### 2026-10-05

- [x] Filed the handoff and seed-analysis results in Drive `BML/KIRC_Hetionet_Project/2026-10-05_handoff_and_seed_analysis` (small files only; see its MANIFEST.md)
- [x] Wrote the experiment-decision summary to Notion and to `analysis/decisions.md` in that folder

### 2026-10-06

- [x] Size-corrected pathway ranking (`scripts/pathway_size_correction.py` -> `results/seed_expansion/pathway_enrichment.tsv`): fold enrichment, hypergeometric p, BH q, plus an annotation-matched null

## 7. Results
Measured on 2026-09-22 from an actual pipeline run.

| Metric | Pathway | Biological Process |
|---|---:|---:|
| Context nodes | 1,822 | 11,381 |
| Gene-context edges | 84,372 | 559,504 |
| KIRC mapped genes | 8,947 | 14,742 |
| Subgraph nodes | 19,456 | 29,631 |
| Subgraph edges | 524,711 | 999,523 |

Result files:

```text
results/
├── pathway/
│   ├── context_nodes.tsv
│   ├── gene_context_edges.tsv
│   └── subgraph_nodes.tsv
├── biological_process/
│   ├── context_nodes.tsv
│   ├── gene_context_edges.tsv
│   └── subgraph_nodes.tsv
└── comparison/
    ├── context_comparison.tsv
    └── context_comparison_detail.tsv
```

Notes:

- `Gene-context edges` counts every edge of the metaedge in Hetionet; `Subgraph edges` counts only those whose gene is KIRC-mapped.
- `Subgraph nodes` / `Subgraph edges` are the NetworkX subgraph after isolate removal: Gene + Disease + the context's nodes, joined by `GiG`, `Gr>G`, `DaG`, `DuG`, `DdG` and the context metaedge. Gene nodes are restricted to the KIRC-mapped `Gene::Entrez` set.

## 8. Problems / Issues
**Standing issues** (re-checked every run; dated entries below are per-run findings.)

- **RESOLVED (2026-09-22): the gene mapping is now the collaborator's standardization.**
  - Where: `data/external/collab/`, loaded by `gene_mapping.load_collaborator_mapping()`.
  - Cause: this pipeline had been deriving its own Ensembl -> Entrez mapping by
    joining HGNC's `ensembl_gene_id` column. That join is wrong for TCGA data:
    TCGA is quantified against a frozen GENCODE v36 build while HGNC tracks the
    current Ensembl release, so every gene whose Ensembl ID was reassigned in
    between was dropped without warning - 33 genes, **SOD2 among them**.
  - Resolved: yes. The collaborator's files are authoritative and read-only.
    Their chain goes Ensembl -> GENCODE v36 `hgnc_id` -> current HGNC -> Entrez,
    which is immune to Ensembl ID drift. All 22 counts in their summary were
    re-verified against the delivered files before adoption.

- **RESOLVED (2026-09-22): `_PAR_Y` identifiers were silently dropped.**
  - Where: `gene_mapping.strip_ensembl_version()`.
  - Cause: the version-stripping regex anchored on `\.\d+$`, but GENCODE writes
    pseudoautosomal duplicates as `ENSG00000002586.20_PAR_Y` - the version is not
    at the end of the string. All 44 such ids survived unstripped and then failed
    every join, landing in `no_hgnc` with no error.
  - Resolved: yes, the suffix is stripped first. The bug only ever affected the
    fallback path, which is no longer used, but it would have recurred.

- **RESOLVED (2026-09-21): the KIRC expression matrix is downloaded directly.**
  - Cause: an earlier run could not reach TCGA/GDC.
  - Resolved: the original notebook's `fetch_xena()` endpoint pair was adopted;
    `gdc.xenahubs.net` is still unreachable here but the S3 fallback serves the
    same files.

- **RESOLVED (2026-09-21): the baseline `Download_and_subgraph.ipynb` is in hand.**
  - Resolved: the original was supplied and the notebook was rebuilt from it;
    14 of its cells are reused verbatim.

- **RESOLVED (2026-09-22): the survival endpoint's encoding is now established.**
  - Where: `TCGA-KIRC.survival.tsv.gz`, loaded by `kirc.load_kirc_survival()`.
  - Cause: the collaborator's own field guide flagged both `OS` and `OS.time` as
    unverified - "Xena's usual convention is 1=death, 0=censored, but confirm this
    build's ETL encoding" and "no unit stated in the metadata". Experiment 4
    (survival validation) inverts entirely if the flag is backwards.
  - Resolved: yes, cross-checked against the clinical table rather than assumed.
    `OS=1` is `vital_status == Dead` for 336/336 rows and `OS=0` is `Alive` for
    608/608 - an exact split. `OS.time` is in days: it equals `days_to_death` for
    all 336 deaths and `days_to_last_follow_up` for all 608 censored rows. The
    conventional reading holds. Recorded in the loader's docstring.

- **Note: `duplicate_target` rows keep their `hetionet_gene_id`.**
  - Where: `kirc_gene_mapping_all.tsv`.
  - Cause: by design - the 26 dropped duplicate features stay in the table for
    traceability and retain the node id they resolved to. Filtering the table on
    `hetionet_gene_id != ""` therefore yields 19,451 rows for 19,425 nodes.
  - Resolved: guarded. `gene_mapping.representative_rows()` selects on
    `mapping_status` instead and asserts uniqueness. This pipeline reads the
    delivered node file, which is already de-duplicated, so it was never exposed.

- **Open: `GSTT1` (`Gene::2952`) has no KIRC expression row.**
  - Where: the 699 genes linked to `Disease::DOID:263` (kidney cancer).
  - Cause: not a defect. GSTT1 is annotated only on a GRCh38 alternate locus, so
    it is absent from the primary-assembly GENCODE quantification TCGA uses. It
    is also a well-known copy-number-variable gene.
  - Resolved: open by design. **Graph-only experiments have 699 positives;
    expression-based experiments have 698.** Fix the label set per experiment
    before scoring, and state which one each result used.

- **Open: the standardized expression matrix has not been delivered.**
  - Where: `kirc_expression_standardized.tsv.gz` (19,297 x 605) exists in the
    collaborator's run but was not among the shared files.
  - Cause: the shared folder targets graph construction, so only the gene-node
    files were handed over.
  - Resolved: open. It can be reconstructed from the raw matrix plus
    `kirc_gene_mapping_all.tsv`, except for the vial choice on 4 patients that
    carry both `01A` and `01B` - that selection lives in the unshared
    `kirc_sample_metadata.tsv`. Request the file rather than guessing.

- **Open: 22,150 KIRC features carry an Entrez ID with no Hetionet Gene node.**
  - Cause: expected. Hetionet v1.0 models 20,945 genes; KIRC covers 60,660
    GENCODE features including pseudogenes, lncRNAs and other non-coding
    biotypes Hetionet never included.
  - Resolved: recorded, no action.

- **Open: Drive holds only part of the project until the first Colab run.**
  - Resolved: open - run `notebooks/OPEN_IN_COLAB.ipynb` once in Colab.

- **Open (2026-10-05): the handoff gene file is the older 19,416-gene version.**
  - Where: `team_handoff_kirc_subgraph/kirc_hetionet_gene_nodes.tsv` vs the 19,425-gene
    file in `data/external/collab/`.
  - Cause: the handoff predates the NCBI verification step. It is an exact subset of the
    final file; the 9 genes it lacks are the NCBI-linked ones (8 merged-Entrez, 1 Ensembl
    cross-reference). Its README therefore quotes 19,416 / 19,289 where the collaborator's
    summary says 19,425 / 19,297.
  - Impact: negligible. Running the handoff script on both files gives identical
    pathway, BP and seed counts; the final file adds 1 gene and 1 `DuG` edge.
  - Resolved: open - pick one as the paper's reference and say so. The final file matches
    the collaborator's own summary.

- **Open (2026-10-05): 486 of the 698 seed genes come only from `DuG`/`DdG`.**
  - Where: kidney-cancer (`DOID:263`) seed set; see `results/seed_expansion/`.
  - Cause: Hetionet derives `DuG`/`DdG` from STARGEO differential-expression
    meta-analyses, i.e. from disease-vs-normal expression, not from curated
    disease-gene association (`DaG`). Seeds built on them are partly an expression result
    already.
  - Evidence it matters: enriched pathways from the `DaG`-only seed (212 genes) and the
    `DuG`/`DdG`-only seed (486 genes) overlap at Jaccard 0.015 (BP 0.01) - they point at
    nearly unrelated biology - and `DaG` alone reproduces 73% of the pathway and 72% of
    the BP result from the full seed.
  - Resolved: open - needs a decision on which edge types define the seed.

- **Note (2026-10-05): all 17 known ccRCC genes checked are already seed genes.**
  - Where: `results/seed_expansion/known_gene_check.tsv`.
  - Consequence: a known-driver list cannot serve as held-out validation against this
    seed; and `BAP1` has 0 pathways, `PBRM1`/`SETD2`/`KDM5C` 2 each, so the chromatin-
    remodelling drivers are essentially unreachable through the pathway route (BP: 17-92).

- **Open (2026-10-05): the handoff zips and large TSVs are not in Drive.**
  - Where: `BML/KIRC_Hetionet_Project/2026-10-05_handoff_and_seed_analysis/`.
  - Cause: the Drive connector accepts file content only as inline text, so the two ~5 MB
    zips and the TSVs above ~100 KB cannot be pushed from this session.
  - Resolved: open - drag the two zips into that folder by hand. `MANIFEST.md` lists the
    sha256 of every file so the upload can be checked.

- **Note (2026-10-06): FDR counts in the funnel tables now use the stricter denominator.**
  - Where: `results/seed_expansion/funnel_*.tsv`, `seed_definition_comparison.tsv`,
    `seed_definition_overlap.tsv` (`scripts/seed_expansion.py`, function `enriched`).
  - Cause: Benjamini-Hochberg had been applied only over contexts that share at least one gene
    with the seed (1,230 pathways, 6,200 BP). Every context is a test (1,822 / 11,381), so the
    p-values are now padded with 1 for the contexts that have no seed gene.
  - Effect: pathways 217 -> 182 (candidates 4,219 -> 4,123); BP 861 -> 653 (11,025 -> 10,847).
    Seed-definition table: all 182 / 653, DaG-only 236 / 877, DuG/DdG-only 33 / 78 enriched
    pathways / BP. The ranking order is identical under both conventions; only the cut-off moves.
  - Resolved: yes - results regenerated; Notion page and Drive `decisions.md`, funnel and
    seed-definition files carry the new figures. The looser counts (217 / 861) are still reported
    next to the stricter ones by `pathway_size_correction.py`.

### 2026-09-21

- [x] Gene ID standardization implemented (Ensembl -> HGNC -> Entrez -> `Gene::Entrez`)
- [x] Hetionet Gene mapping validation implemented and run
- [x] Pathway / Biological Process configuration (`GRAPH_CONTEXT` + `CONTEXT_CONFIG`)
- [x] Common pipeline implementation (`get_context_nodes` / `get_context_edges` / `run_context_experiment`)
- [x] Pathway experiment executed
- [x] Biological Process experiment executed
- [x] Automatic context comparison written to `results/comparison/context_comparison.tsv`

### 2026-09-22

- [x] Gene ID standardization implemented (Ensembl -> HGNC -> Entrez -> `Gene::Entrez`)
- [x] Hetionet Gene mapping validation implemented and run
- [x] Pathway / Biological Process configuration (`GRAPH_CONTEXT` + `CONTEXT_CONFIG`)
- [x] Common pipeline implementation (`get_context_nodes` / `get_context_edges` / `run_context_experiment`)
- [x] Pathway experiment executed
- [x] Biological Process experiment executed
- [x] Automatic context comparison written to `results/comparison/context_comparison.tsv`

### 2026-10-05

- [x] Filed the handoff and seed-analysis results in Drive `BML/KIRC_Hetionet_Project/2026-10-05_handoff_and_seed_analysis` (small files only; see its MANIFEST.md)
- [x] Wrote the experiment-decision summary to Notion and to `analysis/decisions.md` in that folder

### 2026-10-06

- [x] Size-corrected pathway ranking (`scripts/pathway_size_correction.py` -> `results/seed_expansion/pathway_enrichment.tsv`): fold enrichment, hypergeometric p, BH q, plus an annotation-matched null

## 7. Results
Measured on 2026-09-22 from an actual pipeline run.

| Metric | Pathway | Biological Process |
|---|---:|---:|
| Context nodes | 1,822 | 11,381 |
| Gene-context edges | 84,372 | 559,504 |
| KIRC mapped genes | 8,947 | 14,742 |
| Subgraph nodes | 19,456 | 29,631 |
| Subgraph edges | 524,711 | 999,523 |

Result files:

```text
results/
├── pathway/
│   ├── context_nodes.tsv
│   ├── gene_context_edges.tsv
│   └── subgraph_nodes.tsv
├── biological_process/
│   ├── context_nodes.tsv
│   ├── gene_context_edges.tsv
│   └── subgraph_nodes.tsv
└── comparison/
    ├── context_comparison.tsv
    └── context_comparison_detail.tsv
```

Notes:

- `Gene-context edges` counts every edge of the metaedge in Hetionet; `Subgraph edges` counts only those whose gene is KIRC-mapped.
- `Subgraph nodes` / `Subgraph edges` are the NetworkX subgraph after isolate removal: Gene + Disease + the context's nodes, joined by `GiG`, `Gr>G`, `DaG`, `DuG`, `DdG` and the context metaedge. Gene nodes are restricted to the KIRC-mapped `Gene::Entrez` set.

## 8. Problems / Issues
**Standing issues** (re-checked every run; dated entries below are per-run findings.)

- **RESOLVED (2026-09-22): the gene mapping is now the collaborator's standardization.**
  - Where: `data/external/collab/`, loaded by `gene_mapping.load_collaborator_mapping()`.
  - Cause: this pipeline had been deriving its own Ensembl -> Entrez mapping by
    joining HGNC's `ensembl_gene_id` column. That join is wrong for TCGA data:
    TCGA is quantified against a frozen GENCODE v36 build while HGNC tracks the
    current Ensembl release, so every gene whose Ensembl ID was reassigned in
    between was dropped without warning - 33 genes, **SOD2 among them**.
  - Resolved: yes. The collaborator's files are authoritative and read-only.
    Their chain goes Ensembl -> GENCODE v36 `hgnc_id` -> current HGNC -> Entrez,
    which is immune to Ensembl ID drift. All 22 counts in their summary were
    re-verified against the delivered files before adoption.

- **RESOLVED (2026-09-22): `_PAR_Y` identifiers were silently dropped.**
  - Where: `gene_mapping.strip_ensembl_version()`.
  - Cause: the version-stripping regex anchored on `\.\d+$`, but GENCODE writes
    pseudoautosomal duplicates as `ENSG00000002586.20_PAR_Y` - the version is not
    at the end of the string. All 44 such ids survived unstripped and then failed
    every join, landing in `no_hgnc` with no error.
  - Resolved: yes, the suffix is stripped first. The bug only ever affected the
    fallback path, which is no longer used, but it would have recurred.

- **RESOLVED (2026-09-21): the KIRC expression matrix is downloaded directly.**
  - Cause: an earlier run could not reach TCGA/GDC.
  - Resolved: the original notebook's `fetch_xena()` endpoint pair was adopted;
    `gdc.xenahubs.net` is still unreachable here but the S3 fallback serves the
    same files.

- **RESOLVED (2026-09-21): the baseline `Download_and_subgraph.ipynb` is in hand.**
  - Resolved: the original was supplied and the notebook was rebuilt from it;
    14 of its cells are reused verbatim.

- **RESOLVED (2026-09-22): the survival endpoint's encoding is now established.**
  - Where: `TCGA-KIRC.survival.tsv.gz`, loaded by `kirc.load_kirc_survival()`.
  - Cause: the collaborator's own field guide flagged both `OS` and `OS.time` as
    unverified - "Xena's usual convention is 1=death, 0=censored, but confirm this
    build's ETL encoding" and "no unit stated in the metadata". Experiment 4
    (survival validation) inverts entirely if the flag is backwards.
  - Resolved: yes, cross-checked against the clinical table rather than assumed.
    `OS=1` is `vital_status == Dead` for 336/336 rows and `OS=0` is `Alive` for
    608/608 - an exact split. `OS.time` is in days: it equals `days_to_death` for
    all 336 deaths and `days_to_last_follow_up` for all 608 censored rows. The
    conventional reading holds. Recorded in the loader's docstring.

- **Note: `duplicate_target` rows keep their `hetionet_gene_id`.**
  - Where: `kirc_gene_mapping_all.tsv`.
  - Cause: by design - the 26 dropped duplicate features stay in the table for
    traceability and retain the node id they resolved to. Filtering the table on
    `hetionet_gene_id != ""` therefore yields 19,451 rows for 19,425 nodes.
  - Resolved: guarded. `gene_mapping.representative_rows()` selects on
    `mapping_status` instead and asserts uniqueness. This pipeline reads the
    delivered node file, which is already de-duplicated, so it was never exposed.

- **Open: `GSTT1` (`Gene::2952`) has no KIRC expression row.**
  - Where: the 699 genes linked to `Disease::DOID:263` (kidney cancer).
  - Cause: not a defect. GSTT1 is annotated only on a GRCh38 alternate locus, so
    it is absent from the primary-assembly GENCODE quantification TCGA uses. It
    is also a well-known copy-number-variable gene.
  - Resolved: open by design. **Graph-only experiments have 699 positives;
    expression-based experiments have 698.** Fix the label set per experiment
    before scoring, and state which one each result used.

- **Open: the standardized expression matrix has not been delivered.**
  - Where: `kirc_expression_standardized.tsv.gz` (19,297 x 605) exists in the
    collaborator's run but was not among the shared files.
  - Cause: the shared folder targets graph construction, so only the gene-node
    files were handed over.
  - Resolved: open. It can be reconstructed from the raw matrix plus
    `kirc_gene_mapping_all.tsv`, except for the vial choice on 4 patients that
    carry both `01A` and `01B` - that selection lives in the unshared
    `kirc_sample_metadata.tsv`. Request the file rather than guessing.

- **Open: 22,150 KIRC features carry an Entrez ID with no Hetionet Gene node.**
  - Cause: expected. Hetionet v1.0 models 20,945 genes; KIRC covers 60,660
    GENCODE features including pseudogenes, lncRNAs and other non-coding
    biotypes Hetionet never included.
  - Resolved: recorded, no action.

- **Open: Drive holds only part of the project until the first Colab run.**
  - Resolved: open - run `notebooks/OPEN_IN_COLAB.ipynb` once in Colab.

- **Open (2026-10-05): the handoff gene file is the older 19,416-gene version.**
  - Where: `team_handoff_kirc_subgraph/kirc_hetionet_gene_nodes.tsv` vs the 19,425-gene
    file in `data/external/collab/`.
  - Cause: the handoff predates the NCBI verification step. It is an exact subset of the
    final file; the 9 genes it lacks are the NCBI-linked ones (8 merged-Entrez, 1 Ensembl
    cross-reference). Its README therefore quotes 19,416 / 19,289 where the collaborator's
    summary says 19,425 / 19,297.
  - Impact: negligible. Running the handoff script on both files gives identical
    pathway, BP and seed counts; the final file adds 1 gene and 1 `DuG` edge.
  - Resolved: open - pick one as the paper's reference and say so. The final file matches
    the collaborator's own summary.

- **Open (2026-10-05): 486 of the 698 seed genes come only from `DuG`/`DdG`.**
  - Where: kidney-cancer (`DOID:263`) seed set; see `results/seed_expansion/`.
  - Cause: Hetionet derives `DuG`/`DdG` from STARGEO differential-expression
    meta-analyses, i.e. from disease-vs-normal expression, not from curated
    disease-gene association (`DaG`). Seeds built on them are partly an expression result
    already.
  - Evidence it matters: enriched pathways from the `DaG`-only seed (212 genes) and the
    `DuG`/`DdG`-only seed (486 genes) overlap at Jaccard 0.015 (BP 0.01) - they point at
    nearly unrelated biology - and `DaG` alone reproduces 73% of the pathway and 72% of
    the BP result from the full seed.
  - Resolved: open - needs a decision on which edge types define the seed.

- **Note (2026-10-05): all 17 known ccRCC genes checked are already seed genes.**
  - Where: `results/seed_expansion/known_gene_check.tsv`.
  - Consequence: a known-driver list cannot serve as held-out validation against this
    seed; and `BAP1` has 0 pathways, `PBRM1`/`SETD2`/`KDM5C` 2 each, so the chromatin-
    remodelling drivers are essentially unreachable through the pathway route (BP: 17-92).

- **Open (2026-10-05): the handoff zips and large TSVs are not in Drive.**
  - Where: `BML/KIRC_Hetionet_Project/2026-10-05_handoff_and_seed_analysis/`.
  - Cause: the Drive connector accepts file content only as inline text, so the two ~5 MB
    zips and the TSVs above ~100 KB cannot be pushed from this session.
  - Resolved: open - drag the two zips into that folder by hand. `MANIFEST.md` lists the
    sha256 of every file so the upload can be checked.

- **Note (2026-10-06): the FDR counts in the funnel tables use a less conservative denominator.**
  - Where: `results/seed_expansion/funnel_*.tsv` and `seed_definition_comparison.tsv`
    (`scripts/seed_expansion.py`, function `enriched`), and the figures quoted from them
    (217 pathways / 4,219 candidates; 861 BP / 11,025 candidates).
  - Cause: Benjamini-Hochberg was applied only over contexts that share at least one gene with
    the seed (1,230 pathways, 6,200 BP). Counting every context as a test (1,822 / 11,381) is
    stricter. Both are used in practice; the first matches common ORA tools.
  - Effect: pathways 217 -> 182 (candidates 4,219 -> 4,123); BP 861 -> 653 (11,025 -> 10,847).
    The ranking order is identical under both; only the cut-off moves.
  - Resolved: open - pick one convention for the paper. `pathway_size_correction.py` uses the
    stricter one and reports both counts.

### 2026-09-21

- **3 Ensembl IDs appear on more than one row (one Ensembl -> several Entrez).**
  - Where: gene_mapping.validate_gene_mapping()
  - Cause: Ensembl/Entrez identifier systems are not 1:1
  - Resolved: recorded, not auto-corrected
- **22,110 Entrez IDs have no Hetionet Gene node.**
  - Where: gene_mapping.validate_gene_mapping()
  - Cause: Ensembl/Entrez identifier systems are not 1:1
  - Resolved: recorded, not auto-corrected

### 2026-09-22

- none recorded in this run

## 9. Decisions
- Gene graph join key = `Gene::Entrez` (`hetionet_gene_id`); gene symbols are annotation only.
- Pathway edge = `GpPW`; Biological Process edge = `GpBP`. Nothing else is read for a context.
- Pathway and Biological Process share one pipeline parameterised by `GRAPH_CONTEXT`, never two code paths.
- Results are written per context under `results/<context>/`, so the two runs cannot overwrite each other.
- The gene mapping is the collaborator's standardization run, used as-is. This pipeline does not re-derive it and never writes to those files.
- Ensembl -> Entrez is bridged by GENCODE v36 `hgnc_id`, not by Ensembl ID and not by gene symbol.
- `graph` (19,425) and `expression` (19,297) are separate gene universes; every result states which one it used.
- Hetionet edges are pulled from the git-lfs media endpoint and checksum-verified; the plain raw URL serves a 133-byte LFS pointer.

## 10. Next Steps
- [ ] Candidate selection
- [ ] Variance calculation
- [ ] Pathway-based scoring
- [ ] Biological Process-based scoring
- [ ] Candidate ranking
- [ ] Machine learning baseline
- [ ] GNN model
- [ ] Biomarker selection / disease-based scoring

## 11. Change Log
### 2026-09-21

- Added `config.py` with `GRAPH_CONTEXT` / `CONTEXT_CONFIG`
- Added Gene ID standardization (Ensembl -> HGNC -> Entrez -> Hetionet)
- Added gene mapping validation and issue reporting
- Added the shared Gene -> Context pipeline and dual-context driver
- Added per-context result output and automatic comparison
- Added README progress tracking and Google Drive backup

### 2026-09-22

- Documented the collaborator/repo boundary and the handoff file set in section 3
- Verified shared sources by SHA-256 against the collaborator's summary: Hetionet nodes, Hetionet edges and the KIRC matrix are byte-identical
- Aligned TSV reading on `keep_default_na=False`, matching the collaborator's scripts

### 2026-10-05

- Added `notebooks/run_handoff_in_colab.ipynb`: runs the teammate's `apply_kirc_to_hetionet_subgraph.py` from the uploaded zip in Colab, checks the result against the handoff README's expected values and packs `output.zip`

### 2026-10-06

- Added `scripts/pathway_size_correction.py` and `results/seed_expansion/pathway_enrichment.tsv` (1,822 pathways) and `size_correction_summary.tsv`
