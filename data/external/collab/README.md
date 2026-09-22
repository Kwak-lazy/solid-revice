# KIRC–Hetionet 그래프 구축용 공유 파일

이 폴더에는 TCGA-KIRC 유전자 발현 feature를 Hetionet Gene 노드(`Gene::Entrez ID`)와 연결한 결과가 들어 있습니다.

## 파일

- `kirc_hetionet_gene_nodes.tsv`: KIRC와 연결된 Hetionet Gene 노드 19,425개. 그래프 구축의 핵심 입력이며 `hetionet_gene_id` 열로 Hetionet 노드와 조인합니다.
- `hetionet_genes_absent_from_kirc.tsv`: Hetionet에는 있지만 KIRC 발현값이 없는 Gene 노드 1,520개와 그 사유입니다. 그래프에는 남길 수 있지만 발현 feature로는 사용할 수 없습니다.
- `kirc_gene_mapping_all.tsv`: KIRC 전체 60,660개 feature의 Ensembl–HGNC–Entrez–Hetionet 매핑표입니다. ID와 예외를 역추적할 때 사용합니다.
- `kirc_gene_standardization_summary.md`: 사용한 참조자료, 매핑·중복 처리 규칙, 단계별 개수와 검산 결과입니다.

## 사용 시 주의

- 이름(symbol)이 아니라 `Gene::<Entrez ID>` 형식의 `hetionet_gene_id`로 조인합니다.
- `kirc_hetionet_gene_nodes.tsv`의 19,425개 중 `all_samples_zero == True`인 128개는 최종 발현행렬에 포함되지 않습니다.
- 그래프와 실제 발현 feature를 함께 사용할 때는 `all_samples_zero == False`인 19,297개를 사용합니다.

