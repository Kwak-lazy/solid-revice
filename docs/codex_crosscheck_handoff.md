# KIRC–Hetionet 지금까지의 작업 정리 (codex 교차검증용)

작성일 2026-11-09. 저장소 `Kwak-lazy/solid-revice`, 브랜치 `claude/determined-lovelace-mj90py`.
이 문서의 수치는 모두 세션에서 직접 계산한 값이고, 출처 파일과 재현 방법을 같이 적었습니다.
**재현하지 못하는 것과 확인하지 못한 것은 8장에 모았습니다. 먼저 읽어 주세요.**

## 1. 목적과 범위

- 목표: 신장암(KIRC, `Disease::DOID:263`)에서 지식 그래프(Hetionet)로 후보 gene을 고르고, 발현 variance와 섞어 평가하는 실험의 준비 단계.
- 이 단계에서 하지 않는 것: variance 계산, gene 순위 모델, ML, GNN. (발현 행렬 `19,297 × 605`와 `kirc_sample_metadata.tsv`가 아직 안 옴.)
- 협업자와의 경계: gene ID 표준화(협업자 산출물)는 그대로 쓰고, 이 저장소는 그래프 구성과 seed·pathway 분석을 맡는다. `data/external/collab/`는 읽기 전용.

## 2. 데이터 출처

| 대상 | 내용 | 검증 |
|---|---|---|
| Hetionet v1.0 | 47,031 nodes / 2,250,197 edges. edges 파일은 Git LFS라서 `media.githubusercontent.com`에서 받음 | sha256 `611f0411…` (README 3장) |
| TCGA-KIRC (UCSC Xena) | 60,660 features × 610 samples, `log2(TPM+1)`, GENCODE v36. 생존: OS 1=사망 0=중도절단, OS.time 일 | 임상표와 대조 |
| gene ID 매핑 | 협업자 표준화를 그대로 사용: Ensembl(GENCODE v36) → `hgnc_id` → 현재 HGNC → Entrez → `Gene::Entrez` | 협업자 파일 개수 검증 |
| gene universe | 그래프 19,425개, 발현 19,297개(서로 다른 집합). 핸드오프 파일은 구버전 19,416/19,289 | 구버전은 부분집합, pathway·seed 수치는 같음 |

- 팀원 핸드오프 스크립트 `apply_kirc_to_hetionet_subgraph.py`를 직접 돌린 결과가 받은 `output.zip` 5개 파일과 **바이트 동일**함. 그래프: 31,417 nodes / 1,082,475 edges (`GpPW`와 `GpBP` 모두 포함).
- 원본 노트북의 base subgraph(31,611 nodes / 1,084,689 edges)와는 엣지 집합까지 일치.
- 이 문서의 분석은 모두 위 **핸드오프 그래프** 기준. 그래프 gene 18,080개, Pathway 주석이 있는 gene 8,928개, Pathway 1,822개, BP 11,381개.

## 3. 저장소 구성

- 파이프라인: `config.py`, `src/kirc_hetionet/*`, `scripts/run_pipeline.py` (`GRAPH_CONTEXT`로 Pathway/BP 전환, 결과는 `results/<context>/`).
- 분석 스크립트(핸드오프 `output/`를 `--graph-dir`로 받음):
  - `scripts/seed_expansion.py`: hop 확장, 필터 funnel, seed 정의 비교, 알려진 gene 점검
  - `scripts/pathway_size_correction.py [--seed DaG|all]`: 크기 보정 순위
  - `scripts/candidate_filter.py`: D5 필터(DaG seed, Pathway 전용)
- 결과: `results/seed_expansion/*.tsv` (파일 목록은 아래 7장).
- README는 진행 로그. **2026-11-09에 README 중복 오류를 고쳤음(8장 참고).**
- 외부 문서: Notion 결정 페이지 `3f07eb0b-6885-81d6-8a9b-d659fd2f235f`(1~8장), Drive `BML/KIRC_Hetionet_Project/2026-10-05_handoff_and_seed_analysis`(`MANIFEST.md`에 파일별 sha256).

## 4. 지금까지 정해진 것

사용자가 결정한 것:
1. **seed = `DaG` 연결 gene만** (212개, 그중 Pathway 주석이 있는 gene 165개). 이유: `DuG`/`DdG`는 STARGEO 발현 메타분석에서 온 엣지라서 발현 기반 실험과 순환이 생김. 전체 698개(`DaG`+`DuG`+`DdG`)는 민감도 분석으로만 사용.
2. **경로는 Pathway만 사용. Biological Process는 제외.**
3. 필터 기준(임시 확정): BH q<0.05(전체 1,822개를 검정으로 센 값), fold≥2, 크기 10~500, seed 겹침≥3. gene 단계는 통과 pathway 2개 이상이 지지하는 gene만.

아직 정해지지 않은 것: 식 버전(원저 DWPC냐 단순화 식이냐), w, β, 평가 지표(AUROC냐 C-index냐), 실험 구성, 모델. (6장)

## 5. 핵심 수치 (출처 파일은 7장)

**seed와 도달 범위 (전체 seed 698 기준, `hop_expansion.tsv`)**
- seed: `DaG` 212 · `DuG` 250 · `DdG` 250, 합집합 698. 그중 Pathway 주석 460, BP 주석 626.
- seed가 닿는 Pathway 1,230(전체 1,822), BP 6,200(전체 11,381). 이웃 gene(`GiG`/`Gr>G`) 8,123.
- Pathway 경유 후보 gene 8,349개(주석 gene 8,928개의 94%). 필터 없이는 후보가 사실상 전부라 필터링이 필수.
- 알려진 ccRCC gene 17개(`known_gene_check.tsv`)는 모두 이미 seed. 정답 목록으로 못 씀.

**FDR 수치 (전체 검정 수를 센 엄격한 기준)**
- 전체 seed: 유의 Pathway 182개(후보 4,123), BP 653개(후보 10,847). `DaG`만: 236 / 877. `DuG`/`DdG`만: 33 / 78.
- 이전에 쓰던 217/4,219, 861/11,025는 seed와 겹치는 context만 검정으로 센 느슨한 기준이었고, 지금은 쓰지 않음. 두 기준 모두 `size_correction_summary*.tsv`에 같이 있음.

**크기 보정 (DaG seed, `size_correction_summary_DaG.tsv`)**
- Spearman(크기, raw seed 개수) 0.629. raw 상위 20개 크기 중앙값 275, 보정 후 70.
- q<0.05 236개 → D5 필터 후 **221개**. 후보 gene 3,066개, 통과 pathway 2개 이상 지지 1,651개 (`d5_filter_summary.tsv`).
- 크기 구간별 통과 수: 10~50 107, 51~100 68, 101~200 42, 201~500 4. q 통과인데 걸린 것: 크기 10 미만 12개, 500 초과 3개, 겹침 3 미만 10개(일부 중복).
- annotation 구간(9개)을 맞춘 null(2,000회, seed 20261006)의 기대값이 균등 가정의 중앙값 1.87배(상위 20개).
- 알려진 pathway의 순위(보정/raw): VEGF 14/19, mTOR 24/60, HIF-1-α TF network 71/81, HIF-2-α TF network 13/69.

**통계 방법 일치 확인**
- 협업자 설명 이미지의 수식(`p_j = Σ_{i=k_j}^{min(K,n_j)} C(K,i)C(N−K,n_j−i)/C(N,n_j)`, N=8,928, K=seed 수, n_j=pathway 크기)을 `math.comb`로 구현해 1,822개 pathway 전부에서 `scipy.stats.hypergeom.sf`로 구한 p와 비교: 최대 상대차 4e-6(저장 시 6자리 반올림), p<0.05 개수 동일(전체 seed 379, DaG 363).
- 협업자 쪽 seed 수 460은 **전체 seed(698)** 에서 나오는 값. DaG만이면 165. 협업자의 0930 랩미팅 페이지는 seed를 "DaG"로 정의하고 "약 711개"라 적었는데 이 수는 재현되지 않음(전체 Hetionet에서 DaG 213, 합집합 699, 핸드오프 그래프에서 212/698).
- 협업자 `kirc_seed_pathway_count_frequency.tsv`: 전체 seed(698)에서 232칸 모두 일치(연결 Pathway 0개 238명, 최대 231, 연결 총합 6,244). DaG만이면 39칸이 다르고, 0개가 47명(22%), 연결 총합 3,399.

## 6. 4번 실험(KG + variance) 정리

**설계(0930 랩미팅 페이지와 협업자 노트)**: 모든 gene g에 `s_final = β·P(s_KG) + (1−β)·P(Var)`(P는 백분위 순위)를 주고 상위 k개를 평가. β=0이면 variance만, β=1이면 KG만. `s_KG`는 직접 연결(GaD) 항과 경로 항(GpPWpGaD)을 α로 섞음. α, β는 nested CV로 선택, 평가는 5-fold × 5회 반복.

**DWPC 식 — 근거가 있는 것**
- `hetio/hetnetpy`의 `pathtools.py`(24~45행)를 원문으로 확인: 경로의 **엣지마다** 양 끝 노드의 해당 엣지 종류 차수를 하나씩 모아 모두 `-w` 제곱해 곱하고, 경로별 역수를 합산.
- `D → s → p → g` 경로에 적용: 엣지 D–s(DaG)는 `deg_DaG(D)·deg_DaG(s)`, s–p(GpPW)는 `deg(s)·deg(p)`, p–g(GpPW)는 `deg(p)·deg(g)`. D의 차수는 상수라 제외. pathway 차수가 **두 번** 들어감.
- 0930 Notion 페이지의 식 `[deg(g)·deg(p)]^-w · [deg(p)·deg(s)]^-w · [deg_DaG(s)]^-w`가 이 규칙과 같음. 실제 경로 5개에서 코드 규칙과 소수점 6자리까지 일치함을 확인.
- 협업자 노트의 식 `(deg(s)·deg(P)·deg(g))^-w`(노드마다 한 번씩)는 **원저 규칙을 단순화한 다른 점수**. pathway 차수가 한 번, deg_DaG(s)는 없음. 크기 20일 때 경로 하나의 가중치는 한 번 곱하면 0.302, 두 번이면 0.091(w=0.4).
- 차수의 뜻: 그 노드에 붙은 해당 엣지 종류의 엣지 수. deg(s)·deg(g)=속한 pathway 수, deg(p)=pathway 크기(핸드오프 그래프 기준, 예: `MAPK1` 231, `Signaling by VEGF` 115), deg_DaG(s)=연관된 질병 수(중앙값 3, 최대 46).

**예비 측정(저장소에 스크립트 없음, 8장 참고)** — DaG seed, 주석 있는 gene 8,928개 기준
- 경로 점수가 0이 아닌 gene 8,561개(그래프 gene의 47%), 그중 seed 165, 비seed 8,396.
- 점수와 gene의 pathway 수 상관 / 상위 200개 pathway 수 중앙값: 원저 규칙 w=0 0.83/68, w=0.4 0.61/54.5, w=0.5 0.48/13.5, w=0.7 0.27/4, w=1 0.14/4. 협업자 식 w=0.4 0.70/56, w=0.5 0.60/55, w=0.7 0.30/5.
- 숨긴 seed 복원 AUROC(주석 있는 DaG seed 165개를 5-fold로 20%씩 숨기고, 보이는 seed로만 점수 계산, 숨긴 seed vs 비seed 주석 gene 비교, 보이는 seed는 평가에서 제외, 3~5회 반복, rng 20261006): pathway 개수만 0.649, 원저 규칙 w=0.4/0.5/0.7이 0.715/0.717/0.705, 협업자 식 0.689/0.706/0.715.
- 해석: AUROC는 식·w 사이 차이가 작고, 허브 편향은 크게 다름. 따라서 판단 기준은 AUROC가 아니라 허브 편향과 기준선(0.649) 대비 증가분.

**KG를 variance와 섞는 이유(가설)**: variance는 이 환자 데이터의 실제 변동을 보지만 질병 관련성은 모르고, KG는 문헌 지식을 주지만 이 데이터를 보지 않고 허브·주석 편향이 있으며 pathway 없는 gene은 구분하지 못함. 섞으면 서로의 약점을 보완한다는 것이 가설이고, β sweep이 이를 직접 시험함. variance 쪽 성능은 발현 행렬이 없어 측정하지 못함.

**결정이 필요한 것(추천 포함)**
1. 식 버전: 원저 DWPC(deg(p) 두 번, `deg_DaG(s)` 포함)를 추천. 협업자 식을 쓰면 다른 이름을 붙이고 w를 0.7 근처까지 올려 허브 편향 확인.
2. 차수의 정의(엣지 종류별 여부), 3. w(원저 규칙 0.5, 협업자 식 0.7, 0.4는 민감도), 4. 직접 연결(GaD) 항의 처리(숨긴 seed 평가에서는 정보가 없어 α를 튜닝할 수 없음 → 포함 규칙으로 분리 추천),
5. β sweep {0, 0.25, 0.5, 0.75, 1}과 선택 기준(생존 C-index), 6. 점수 0인 gene 처리(하위 동률, 주석 없는 gene 선택 비율 보고),
7. 4번이 3번의 enrichment와 같은 정보를 이중으로 쓰지 않게 분리, 8. 검증은 기준선 대비 증가분으로, seed를 숨길 때 그 gene의 DaG 엣지와 "경로의 s" 역할을 같이 제거,
9. 협업자 문서의 미정의 항목(3-A, 3-B, "seed 교란 두 방향", "동일 모델")과 seed 수(460 vs 165) 확인.

## 7. 파일 → 주장 대응 (교차검증 체크리스트)

| 확인할 주장 | 값 | 파일 / 재현 |
|---|---|---|
| seed 698 = DaG 212 ∪ DuG 250 ∪ DdG 250 | 698 | `hop_expansion.tsv`, `seed_definition_comparison.tsv` |
| 주석 있는 seed | 전체 460, DaG 165 | `size_correction_summary.tsv`, `size_correction_summary_DaG.tsv` (`seeds_in_universe_n`) |
| background N | 8,928 | 위 두 파일 `universe_genes_N` |
| 유의 Pathway 수(q<0.05, 전체 1,822 검정) | 전체 182, DaG 236 | `pathway_enrichment*.tsv`의 `q` 열 |
| D5 필터 통과 | DaG 221, 전체 163 | `d5_filter_summary.tsv` |
| D5 후보 gene | DaG 3,066 → 지지≥2 1,651 | `d5_filter_summary.tsv`, `d5_candidates_pathway.tsv`(1,651행) |
| 상위 지지 gene | MAP2K1 106, PIK3R1 101, HRAS 97, AKT1 96, GRB2 91 | `d5_candidates_pathway.tsv` |
| 전체/DaG q 통과 집합의 Jaccard | 0.467(q<0.05), 0.488(D5 통과) | `pathway_enrichment.tsv`와 `_DaG.tsv` 비교 |
| 무작위 seed에서 q<0.05 | 평균 0.01~0.03개 | 스크립트 없음(8장) |

재현 명령(핸드오프 `output/`가 있을 때):
```
python scripts/seed_expansion.py --graph-dir <핸드오프>/output
python scripts/pathway_size_correction.py --graph-dir <핸드오프>/output --seed DaG
python scripts/candidate_filter.py --graph-dir <핸드오프>/output
```
전체 seed 결과(`pathway_enrichment.tsv`)는 `--seed all`로 재생성하면 바이트 동일.

## 8. 재현하지 못하는 것과 확인하지 못한 것

**먼저 알아야 할 것**
- **핸드오프 `output/`(nodes.tsv, edges.tsv)가 저장소에 없음.** 위 재현 명령은 이 파일이 있어야 돌아감. 파일은 Drive에도 없고(Drive 도구가 텍스트 inline만 받음), 팀원이 준 zip 2개의 `MANIFEST.md`에 sha256만 있음. 스크립트를 돌리려면 zip을 받아 `apply_kirc_to_hetionet_subgraph.py`를 실행해 `output/`를 만들어야 함.
- **4번 실험의 예비 측정(DWPC 허브 편향, 숨긴 seed AUROC, D5 임계값 sweep, 무작위 seed 검증, DuG/DdG 라벨 복원)은 저장소에 스크립트가 없음.** 세션 임시 폴더에서만 돌렸고 컨테이너가 바뀌어 없어짐. 6장의 알고리즘 설명만으로 독립 재구현해 비교하는 것이 가장 의미 있는 교차검증이 됨.
- 이전 세션에서 만든 standalone Colab 노트북(`KIRC_Hetionet_standalone.ipynb`)은 사용자에게 전달만 했고 저장소에 없음.

**확인하지 못한 것**
- Hetionet 논문 본문(eLife·PLoS)은 접근이 막혀 읽지 못함. DWPC 정의는 `hetnetpy`·`hetmatpy` 코드로 확인했고, `w=0.4`를 논문이 썼다는 협업자 설명은 확인하지 못함(`hetmatpy` 기본값은 0.5).
- 협업자 이미지의 "차수"가 노드 전체 엣지 수인지 엣지 종류별인지는 확정할 수 없음.
- 협업자 쪽 p 값과 상위 pathway 목록은 받은 적이 없어 pathway 단위 대조는 못 함. 받은 것은 수식 이미지와 칸별 빈도 파일뿐.
- 생존 C-index는 발현 행렬이 없어 측정하지 못함.
- Notion 0919·0923 페이지는 연결에서 열리지 않았음.

**세션 중 바로잡은 것(같은 실수를 피하려고)**
- FDR 기준: 처음엔 seed와 겹치는 context만 검정으로 세어 217/861이 나왔고, 전체 context 기준(182/653)으로 바로잡아 Notion·Drive·저장소를 모두 갱신.
- BP "후보 1,801개"는 q<0.01 기준 값이었고 q<0.05 기준은 3,687개. (BP는 이후 범위에서 제외.)
- "무관해도 91개가 통과" → 실측 평균 39개(대부분의 pathway는 겹침이 없어 p=1).
- 협업자와 내 결과가 "같다"는 것은 p-value 계산과 순위이지 seed 정의는 아님(460 vs 165).
- "원저 DWPC 식은 제가 정리한 것"이라고 한 표현을 정정: 0930 Notion 페이지의 식이며, 원저 코드와 일치함을 이번에 확인.
- **README 중복 오류(내 실수)**: 2026-10-06 수정(`f66f666`)에서 삽입 위치를 잘못 잡아 197줄이 복제되었고 이후 3개 커밋이 그 위에 쌓였음. 2026-11-09에 깨끗한 판본(`fccb7ef`)에서 의도한 변경만 다시 적용해 바로잡음(552줄). 중복으로 인해 잘못된 위치에 들어갔던 2026-10-06 노트들도 section 8로 옮김.
