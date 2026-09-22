# KIRC-Hetionet Project

## Current Status
Stage: Gene ID + Graph Context Pipeline

Status: In Progress

Last updated: 2026-09-22

Completed:
- Gene ID standardization (collaborator's run adopted as authoritative)
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
| Dataset | Source | Local path |
|---|---|---|
| Hetionet v1.0 nodes (47,031) | `raw.githubusercontent.com/hetio/hetionet` | `data/raw/hetionet-v1.0-nodes.tsv` |
| Hetionet v1.0 edges (2,250,197) | git-lfs **media** endpoint, sha256-verified | `data/raw/hetionet-v1.0-edges.sif.gz` |
| HGNC complete set (approved) | `storage.googleapis.com/public-download-files/hgnc` | `data/raw/hgnc_complete_set.txt` |
| TCGA-KIRC expression, `log2(TPM+1)` | UCSC Xena GDC hub (S3 fallback) | `data/raw/TCGA-KIRC.star_tpm.tsv.gz` |
| TCGA-KIRC clinical | UCSC Xena GDC hub (S3 fallback) | `data/raw/TCGA-KIRC.clinical.tsv.gz` |
| Prior-run gene coverage (reference) | shared collaborator artifact | `data/external/hetionet_genes_absent_from_kirc.tsv` |

Hetionet node kinds (measured): Gene 20,945 · Biological Process 11,381 ·
Side Effect 5,734 · Molecular Function 2,884 · Pathway 1,822 · Compound 1,552 ·
Cellular Component 1,391 · Symptom 438 · Anatomy 402 · Pharmacologic Class 345 ·
Disease 137.

TCGA-KIRC matrix (measured): 60,660 genes x 610 samples — 537 Primary Tumor,
72 Solid Tissue Normal, 1 Additional New Primary. Row index is a
version-suffixed Ensembl ID (`ENSG00000141510.17`).

Download notes:

- The Hetionet edge file is git-lfs. `raw.githubusercontent.com` returns a
  133-byte LFS pointer, so the **media** endpoint is used and the payload is
  checksum-verified.
- Xena's official host (`gdc.xenahubs.net`) is unreachable on some networks;
  `fetch_xena()` falls back to `gdc-hub.s3.us-east-1.amazonaws.com`, which
  serves the same files. Both endpoints are tried in order.

`data/raw/` is git-ignored and re-downloaded on demand.
`results/*/gene_context_edges.tsv` and `results/*/subgraph_edges.tsv` are
git-ignored for size and regenerated by the pipeline.

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

- Replaced the in-house gene mapping with the collaborator's standardization (`data/external/collab/`), loaded read-only and verified on load
- Fixed `_PAR_Y` version stripping (44 ids were silently unmapped)
- Split the gene set into `graph` (19,425) and `expression` (19,297) universes, selectable via `--gene-universe`
- Recomputed both context experiments on the corrected gene set
