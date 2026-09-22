# KIRC 유전자 ID 표준화 요약 (2026-09-15)

## 1. 사용한 파일과 SHA-256

| 파일 | SHA-256 |
|---|---|
| TCGA-KIRC.star_tpm.tsv (압축 해제본) | `00ca38c7a8e969905f545ce16d93d9452652813d368fed911ad8f3aa8d597464` |
| TCGA-KIRC.sample_labels.tsv | `1e3106ee01fc87fd560e9656e7d2b64fc5a3cbad68fe43eb3873d4cc7b8fb29a` |
| TCGA-KIRC.survival.tsv | `ffb1326672599c1a9c563834d266f0c7b28a7fa7c0e2af9f5710a721800b3767` |
| gencode.v36.complete.gtf.gz (GENCODE v36, 2020-09-30) | `81359865bab048361cbeaf91c9441aac4c700b63978e6691d579d78bf2322d30` |
| hgnc_complete_set.txt (2026-09-14 취득) | `9e07bb49393c12aca8bc931a4349010404c24994a8a58745c14676b8c2b74492` |
| hetionet-v1.0-nodes.tsv | `3ad9eb86c031df8bc7d5efefa6edc93562182eb9a6a71e5129d3f18b9d42c50e` |
| hetionet-v1.0-edges.sif.gz | `611f0411ac666be4e0270e54bbca67407e14720bfc3099311f0f5eec25c5a947` |

원본 발현행렬 압축본(TCGA-KIRC.star_tpm.tsv.gz)의 SHA-256은 audit.json 기준 `70ab9bb0faa27cc76003626546c2f87658606a240781f591552d88176de3c3b9`.

## 2. ID 변환 경로와 규칙

```text
KIRC Ensembl ID (GENCODE v36, 버전 포함)
  → 1순위 GENCODE gene 줄의 hgnc_id / 2순위 HGNC ensembl_gene_id == 버전 제거 Ensembl ID
  → HGNC 현재 symbol, entrez_id
  → Hetionet Gene::entrez (nodes.tsv에 존재할 때만)
```

- symbol 변경 분류 판정 순서(코드 고정): previous_symbol → alias_symbol → clone_name_replaced → other_symbol_change
- 클론형 이름 정의: GenBank accession형 `^[A-Z]{1,2}\d{5,6}\.\d+$`, BAC 클론명형 `^[A-Za-z]{2,6}\d*-\d+[A-Z]\d+\.\d+$`, WGS accession형 `^[A-Z]{4}\d{8,}\.\d+$`
- HGNC 미연결 분류 판정 순서: clone_name_unresolved → ambiguous_multicopy(같은 GENCODE 이름이 2행 이상) → no_hgnc_id
- 중복(같은 Hetionet 노드로 가는 feature) 규칙 순서: ① _PAR_Y 행이 있으면 X쪽 행 대표 ② 한 행만 빼고 전 표본 0이면 발현되는 행 대표 ③ 현재 HGNC ensembl_gene_id와 일치하는 locus 대표 ④ 그래도 없으면 최대 발현값. 탈락 행은 mapping_status=duplicate_target으로 보존
- 이름만 일치하는 보조 후보는 NCBI Gene 번호 이력(병합 기록)과 Ensembl 교차참조로 같은 유전자가 확정된 경우에만 연결하고, 이름 자체는 근거로 쓰지 않음 (scripts/verify_auxiliary_candidates_ncbi.py)

## 3. 단계별 개수

| 항목 | 수 |
|---|---:|
| 원본 feature | 60,660 |
| GENCODE gene 줄에 hgnc_id 있음 | 39,014 (그중 현재 HGNC에 없는 번호 123) |
| HGNC 연결 | 42,068 (GENCODE 번호 38,891 + Ensembl 역조회 3,177) |
| Entrez ID 확보 | 41,600 |
| symbol_status unchanged / changed / no_hgnc | 38,172 / 3,896 / 18,592 |
| changed 세부: clone_name_replaced / previous_symbol / alias_symbol / other | 2,963 / 748 / 174 / 11 |
| Ensembl ID 변경(HGNC 현재 ID와 다름) | 51 |
| _PAR_Y 행 / 전 표본 0 행 | 44 / 2,660 |
| Hetionet 연결 행 → 고유 노드 | 19,451 → 19,425 (HGNC→Entrez 직접 19,416 + NCBI 검증 추가 9) |
| 중복 탈락 행 (PAR_Y / 전 표본 0 locus / 현재 HGNC locus 규칙) | 26 (17 / 3 / 6) |
| entrez_not_in_hetionet / hgnc_mapped_no_entrez | 22,150 / 468 |
| HGNC 미연결: clone_name_unresolved / ambiguous_multicopy / no_hgnc_id | 17,234 / 1,194 / 163 |
| Hetionet 노드 중 KIRC feature 없음 | 1,520 (현재 HGNC에 없는 Entrez 1,428 + 이미 연결된 노드로 병합된 유령 노드 8 + KIRC에 없는 유전자 84) |
| 이름 일치 보조 후보 NCBI 검증 | 22개 중 연결 9 (병합 기록 8 + Ensembl 교차참조 1), 이미 연결된 노드로 병합된 유령 노드 8, 보류 3, 다른 유전자 2 |
| kidney cancer(DOID:263) 직접 연결 유전자 | 699개 중 698개 연결, GSTT1(Entrez 2952)은 GRCh38 alternate locus에만 주석되어 KIRC 행렬에 없음 |

## 4. 표준화 발현행렬 (5-2)

| 항목 | 값 |
|---|---|
| 유전자(행) | 19,297 = 연결 대표 행 19,425 − 전 표본 0 128 |
| 표본(열) | 605 = selected_one_vial_per_patient_class == True (tumor 533, normal 72) |
| 행 ID | `Gene::Entrez` (symbol은 kirc_hetionet_gene_nodes.tsv 참조) |
| 값 | 원본 log2(TPM+1) 문자열 그대로. 재정규화·재변환 없음 |
| 605개 표본 안에서 전 표본 0인 유전자 | 2 (DEFB109C, DEFA1; 610 기준으로는 발현 있음. 학습 fold 내 상수 필터에 맡김) |
| 결측·무한대·음수 셀 | 0 |
| 종양·정상 쌍이 있는 환자 | 72명 → 분할 시 patient_id를 group으로 사용 |

## 5. 산출물

| 파일 | 내용 |
|---|---|
| kirc_gene_mapping_all.tsv | 60,660행 전체 매핑표 (36열) |
| kirc_expression_standardized.tsv.gz | 5-2 발현행렬 |
| kirc_sample_metadata.tsv | 5-3 표본 메타데이터 610행 |
| kirc_mapping_exceptions.tsv | 5-4 예외 사유가 있는 행 41,391개 |
| kirc_duplicate_gene_review.tsv | 중복 26그룹 52행 처리 결과 |
| kirc_auxiliary_name_candidates.tsv | 이름 일치 후보 22개 (연결 적용 전 상태) |
| kirc_auxiliary_name_candidates_verified.tsv | 후보 22개의 NCBI 검증 결과와 조치 |
| kirc_hetionet_gene_nodes.tsv | 2번 담당자 전달: KIRC 연결 노드 19,425 |
| hetionet_genes_absent_from_kirc.tsv | 2번 담당자 전달: 발현값 없는 노드 1,520 (병합된 유령 노드 8 포함) |

## 6. 결정 사항과 미결 사항

결정 사항

- 최종 발현행렬은 선택 표본 605개(종양 533, 정상 72)를 사용한다. 같은 환자의 중복 `01B` vial 4개와 `05A` 표본 1개는 제외한다. 유전자 풀은 그래프 연결 대표 행에서 전 표본 0을 제외한 것(A안)이며, baseline 비교도 같은 풀에서 뽑는 것을 전제로 한다.
- Hetionet 중복 26그룹은 `_PAR_Y` 17그룹, 전 표본 0 locus 3그룹, 현재 HGNC locus 6그룹 규칙으로 처리 완료했다.
- `LRTOMT`, `ARMCX5-GPRASP2` 2그룹은 규칙상 현재 HGNC locus를 대표로 선택했다. 탈락 행은 `mapping_status = duplicate_target`으로 보존하고, 두 행 모두 notes에 `higher_expressed_locus_dropped`를 기록했다.
- 전 표본 0 여부는 원본 코호트 610개 표본을 기준으로 판정한다.
- 이름만 일치하는 후보 22개는 NCBI 번호 이력으로 검증했다. 옛 Entrez가 KIRC 행의 현재 Entrez로 병합된 8개와 Ensembl 교차참조가 일치한 RSC1A1은 연결했고, 이미 연결된 노드로 병합된 유령 노드 8개는 absent 목록에 병합 대상을 기록했으며, HGC6.3·C20orf197은 다른 위치로 판정, TIAF1(MYO18A로 병합)·FAM197Y1·MAGEA10-MAGEA5는 번호 근거가 없어 보류했다.
- `DEFB109C`, `DEFA1`은 610개 중 각각 1개 표본에서 발현이 있어 유지했다. 선택된 605개 표본에서는 모두 0이므로 학습 fold에서도 상수 feature로 제거되어 모델 결과에 영향을 주지 않으며, 행렬 단계에서 별도로 제거하지 않는다. 상수 feature 필터는 각 training fold 안에서 적용한다(모델링 단계 규칙).

미결 사항

- (1) 전체 발현 유전자에서 뽑는 baseline(B안) 추가 여부 (2) BRCA·LUAD 확장 시 all_samples_zero와 중복 규칙 ②만 코호트별 재계산

## 7. 재현

```bash
.venv/bin/python scripts/build_kirc_gene_mapping.py          # 2단계
.venv/bin/python scripts/build_kirc_hetionet_mapping.py      # 3~4단계 (후보 산출)
.venv/bin/python scripts/verify_auxiliary_candidates_ncbi.py # 4-2 NCBI 검증 (인터넷 필요; 결과 파일이 있으면 생략 가능)
.venv/bin/python scripts/build_kirc_hetionet_mapping.py      # 3~4단계 재실행 (검증 결과 적용)
.venv/bin/python scripts/build_kirc_standardized_outputs.py  # 5단계
```
