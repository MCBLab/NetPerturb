# Relatório: avaliação do scMalignantFinder para AML

**Data:** outubro de 2026
**Ferramenta avaliada:** [scMalignantFinder](https://github.com/Jonyyqn/scMalignantFinder) (Yu, Li & Chen, *Communications Biology* 2025)
**Dataset de validação:** vanGalen et al. 2019 (`benchmarks/vangalen.seurat.rds`), 44.823 células, 27.899 genes
**Critério de aceite definido pelo usuário:** ≥90% de concordância com a classificação malignant/normal original do vanGalen (`PredictionRefined`), excluindo as 1.133 células "unclear"

## 1. Objetivo

O vanGalen et al. usou, na própria metodologia do artigo, um método de clusterização para separar células malignas de normais em AML (coluna `PredictionRefined`: malignant/normal/unclear). A pergunta era se o **scMalignantFinder**, uma ferramenta de classificação malignant/normal por assinatura de expressão treinada em dados pan-câncer, reproduz essa separação com concordância suficiente (≥90%) para ser adotada em outros datasets do projeto (Naldini, Lasry, BeneytoCalabuig).

## 2. Resultado resumido

| # | Configuração | Agreement | AUC-ROC | Sensibilidade (malignant) | Especificidade (normal) |
|---|---|---|---|---|---|
| 1 | Modelo pré-treinado (pan-câncer, pronto do Zenodo) | 57,81% | 0,592 | 23,28% | 88,85% |
| 2 | Retreino: Naldini+Beneyto, LogisticRegression | 70,26% (74,0% c/ limiar ótimo) | 0,817 | 92,03% | 50,68% |
| 3 | Retreino: Naldini+Beneyto, LogisticRegression, `class_weight="balanced"` | 71,34% | 0,820 | 91,18% | 53,50% |
| 4 | **Retreino: Naldini+Beneyto+Lasry, LogisticRegression (melhor resultado)** | **80,92%** | **0,878** | **84,95%** | **77,30%** |
| 5 | Retreino: 3 cohorts, RandomForest | 75,05% | 0,832 | 86,82% | 64,44% |
| 6 | Retreino: 3 cohorts, XGBoost | 62,76% | 0,697 | 77,26% | 49,74% |
| 7 | Retreino: 3 cohorts, LogisticRegression, lista de genes específica de AML (em vez da lista pan-câncer) | 79,51% | 0,877 | 79,87% | 79,18% |

**Nenhuma configuração atingiu os 90% definidos como critério.** O melhor resultado foi a configuração #4: 80,92% de concordância, AUC 0,878.

## 3. O que foi testado, em ordem

### 3.1 Linha de base: modelo pré-treinado

Primeiro teste, usando o modelo e a lista de 2.707 genes publicados pelos autores (`model.joblib` + `ordered_feature.tsv`, [Zenodo 17888140](https://zenodo.org/records/17888140)), sem nenhum retreino. `norm_type=True` (normalização por contagem total, sem log — é assim que o próprio pacote espera o dado, confirmado lendo o código-fonte).

**Resultado: 57,81% de concordância, AUC 0,592** — praticamente aleatório. O modelo erra sistematicamente para o lado de "normal": 76,7% das células malignas do vanGalen foram chamadas de normais.

Antes de desistir, três hipóteses de ajuste rápido foram eliminadas:
- **Normalização**: o código do pacote confirma que treino e predição usam só `normalize_total`, sem log — já era o que estava sendo feito.
- **Cobertura de genes**: 2.697 dos 2.707 genes do modelo existem no vanGalen (99,6%) — não é o problema.
- **Calibração do limiar de corte**: a AUC de 0,592 já descarta essa hipótese — nenhum limiar levaria a uma acurácia muito melhor que o sorteio, porque a própria probabilidade não carrega sinal discriminante suficiente.

**Conclusão desta etapa:** o modelo pan-câncer pré-treinado não generaliza para blastos de AML. A hipótese mais provável é que ele foi treinado majoritariamente em tumores sólidos (assinaturas de EMT, proliferação epitelial etc.), e blastos de AML se parecem transcricionalmente com progenitores/HSC normais — não têm esse tipo de assinatura.

### 3.2 Retreino com dados reais de AML

O pacote permite treinar o classificador do zero (`train_h5ad_path` + rótulos `Raw_annotation` ∈ {Normal, Malignant, Tumor}). Em vez de usar o modelo pan-câncer, treinamos com os próprios dados de AML do projeto (Zeng et al. 2025 / BoneMarrowMap), mantendo o vanGalen de fora do treino (held-out).

**Fonte dos dados de treino e mapeamento de rótulos.** A coluna `NormalvsMalignant` do `scAML_21cohorts_Annotations_byCell.csv` **não usa o mesmo vocabulário em todos os datasets** — cada estudo tem sua própria convenção, provavelmente porque cada grupo definiu "malignidade" com um método diferente (genotipagem, clusterização, etc.). Foi preciso inspecionar cada um (cruzando com a coluna `Diagnosis`) antes de montar os rótulos de treino:

| Dataset | Vocabulário original | Mapeamento usado | Malignant | Normal | Total |
|---|---|---|---|---|---|
| Naldini_2023 | `MUT`/`WT`/`ND`/`NoCall`/NA — **são chamadas de genotipagem por célula**, não rótulos de malignidade | `MUT`→Malignant, `WT`→Normal; `ND`/`NoCall`/NA descartados (indeterminados) | 17.583 | 2.860 | 20.443 de 107.243 |
| BeneytoCalabuig_2023 | `leukemic`/`healthy`/`Normal`/`unsure`/NA | `leukemic`→Malignant; `healthy` (célula não-maligna dentro de amostra AML) e `Normal` (amostra de doador sadio) → Normal; `unsure`/NA descartados | 37.268 | 27.345 | 64.613 de 88.602 |
| Lasry_2022 | `malignant`/`Control`/`microenvironment` | `malignant`→Malignant; `Control` (amostra saudável) e `microenvironment` (célula não-maligna dentro de amostra AML/AMKL) → Normal | 123.147 | 86.910 | 210.057 de 210.057 (100%, nenhuma célula sem rótulo) |

As contagens brutas de Naldini e BeneytoCalabuig vieram dos objetos completos já baixados do Zenodo. As do Lasry vieram da reconstrução a partir do GEO (GSE185381) feita anteriormente neste projeto — por sorte, os 47 prefixos de amostra já extraídos cobriam **100%** das células do Lasry, não só o subconjunto Adulto+AML usado nos benchmarks do nf_scrank, então não foi necessário baixar nada novo do GEO.

#### 3.2.1 Naldini + Beneyto (85.056 células de treino)

Primeiro retreino, com os dois datasets menores. **Resultado: 70,26% de concordância, AUC 0,817** — salto grande em relação ao modelo pré-treinado (AUC de 0,592 para 0,817), confirmando que o problema era mesmo o domínio de treino, não normalização ou features.

O viés mudou de direção: agora o modelo super-chama "malignant" (sensibilidade 92%, especificidade só 51%). Duas tentativas de corrigir isso:
- **Ajuste de limiar de corte** (em vez do padrão 0,5): usando o limiar ótimo por Youden's J (0,882), a acurácia sobe para 73,98% (sensibilidade 77,4%, especificidade 70,9%).
- **`class_weight="balanced"`** no treino: patch aplicado em `init_core_model()` (não exposto como parâmetro pela API do pacote). Resultado: 71,34%, AUC 0,820 — melhora marginal, dentro do ruído. A proporção de treino (64% malignant / 36% normal) não era desbalanceada o suficiente para esse ajuste fazer diferença grande.

#### 3.2.2 Naldini + Beneyto + Lasry (295.113 células de treino) — melhor resultado

Adicionar o Lasry como terceiro dataset de treino (sem precisar baixar nada novo, como descrito acima) deu o maior salto do projeto: **80,92% de concordância, AUC 0,878**. Sensibilidade 84,95%, especificidade 77,30% — o mais equilibrado de todos os testes. Nesta configuração, o limiar padrão (0,5) já é essencialmente o ótimo (Youden's J devolveu 0,498), ou seja, não há mais ganho fácil ajustando o corte.

### 3.3 Comparação de algoritmos (3 cohorts de treino)

Com o mesmo conjunto de treino (295.113 células, lista de genes pan-câncer original), testamos os três algoritmos que o pacote suporta nativamente:

| Algoritmo | Agreement | AUC-ROC |
|---|---|---|
| **LogisticRegression** | **80,92%** | **0,878** |
| RandomForest (parâmetros padrão) | 75,05% | 0,832 |
| XGBoost (parâmetros padrão) | 62,76% | 0,697 |

A regressão logística regularizada generalizou melhor entre datasets com efeitos de lote/plataforma diferentes do que os métodos baseados em árvore sem ajuste de hiperparâmetros — padrão esperado quando se mistura datasets com forte efeito de lote e não se faz tuning.

### 3.4 Lista de genes específica de AML

Última tentativa: em vez de usar a lista de 2.707 genes pan-câncer dos autores, calculamos DEGs (`scanpy.tl.rank_genes_groups`, Wilcoxon, malignant vs. normal) **separadamente em cada um dos 3 cohorts de treino**, e mantivemos só os genes consistentes (mesma direção de regulação) em pelo menos 2 dos 3 — para evitar que a lista capturasse só efeito de lote de um único dataset. Resultado: 1.415 genes (1.378 efetivamente usados, após interseção com o painel de genes do vanGalen).

**Resultado: 79,51% de concordância, AUC 0,877** — essencialmente igual à lista pan-câncer genérica (80,92%/0,878), sem ganho real. Isso indica que o salto de desempenho já tinha vindo todo do retreino com população de AML real; trocar a lista de genes (mantendo o mesmo treino) não agrega mais sinal.

## 4. Interpretação e conclusão

1. O modelo **pré-treinado** do scMalignantFinder não é adequado para AML — foi treinado predominantemente em biologia de tumor sólido e não reconhece a assinatura de blasto leucêmico.
2. **Retreinar com dados reais de AML resolve a maior parte do problema** (AUC de 0,59 para 0,88), mas não o suficiente para atingir 90% de concordância com o `PredictionRefined` do vanGalen.
3. Nem balanceamento de classes, nem troca de algoritmo (RandomForest/XGBoost), nem uma lista de genes específica de AML melhoraram o melhor resultado (LogisticRegression, 3 cohorts, lista pan-câncer: 80,92%).
4. O teto de ~80-81% provavelmente reflete limitações estruturais que não são fáceis de contornar com os dados que já temos:
   - **Efeito de lote entre datasets**: Naldini, Beneyto, Lasry e vanGalen vêm de protocolos/tecnologias diferentes; um classificador por expressão absorve parte desse ruído como sinal.
   - **Natureza do `PredictionRefined`**: não sabemos ao certo se a classificação original do vanGalen veio de genotipagem direta de mutações (como no Naldini) ou de um método transcricional — se for genotipagem, é um padrão-ouro de natureza diferente de qualquer classificador por expressão, e pode nunca bater 90% mesmo com um modelo "perfeito".

## 5. Recomendação

Com o esforço de ajuste "barato" esgotado, **a ferramenta não atingiu o critério de 90% definido para ser estendida aos outros datasets do projeto** (Naldini, Lasry, BeneytoCalabuig). As alternativas que ficaram de fora deste teste, por serem de escopo maior:
- Correção explícita de efeito de lote entre os datasets de treino antes de treinar (ex.: Harmony/Scanorama na etapa de normalização).
- Incluir mais cohorts de treino do Zeng et al. 2025 (ficamos em 3 dos 21 disponíveis, limitados pelos datasets que já tínhamos com contagens brutas baixadas).
- Investigar a origem exata do `PredictionRefined` do vanGalen, para confirmar se é um alvo atingível por classificação de expressão.

Os datasets já usados no projeto (Naldini, Lasry, BeneytoCalabuig) têm seus próprios rótulos malignant/normal curados (`NormalvsMalignant`), então essa separação já está disponível sem depender do scMalignantFinder.

## 6. Problemas técnicos resolvidos (apêndice)

Para referência futura, caso o ambiente Python precise ser reconstruído ou outro pacote pesado precise ser instalado neste cluster:

1. **`miniconda3` em `/home/lgdqamorim/scratch/` estava corrompido** (símlinks de `python`/`pip` apontando para `python3.13`, que não existe ali). Solução: usar o módulo do cluster `softwares/python/3.10.5-gnu8` + `python3 -m venv` para criar um ambiente limpo em `/home/lgdqamorim/scratch/scmalignantfinder/venv`.
2. **O módulo do cluster define `PYTHONPATH`** apontando para o `site-packages` global dele, e essa entrada fica **antes** da do venv no `sys.path`. Qualquer pacote que já exista ali (pandas antigo, numpy antigo, etc.) sombreia silenciosamente a versão instalada no venv. Esse mesmo venv também não tem `pip` próprio — só funciona via esse `PYTHONPATH`. Solução: manter `PYTHONPATH` setado para todo comando `pip`, e fazer `unset PYTHONPATH` só na hora de **rodar** os scripts Python.
3. **O resolvedor novo do pip** (`pip`’s dependency resolver) entrou em backtracking catastrófico (mais de uma hora sem terminar) tentando resolver as dependências pesadas do scMalignantFinder (pyscenic, squidpy, dask, todas com versão exata fixada). Solução: `pip install --use-deprecated=legacy-resolver`.
4. Passar pacotes extras sem versão fixada (`anndata scipy pandas`) junto com `scMalignantFinder` no mesmo comando fez o resolvedor antigo dar prioridade a eles sobre os pins exatos do próprio pacote, instalando `pandas 1.4.3` em vez do `2.1.4` necessário. Solução: instalar **só** `scMalignantFinder`.
5. Pacotes pequenos já visíveis via o `PYTHONPATH` global (`h5py`, `requests`, `jinja2`, `six`, `certifi` etc.) foram marcados por pip como "already satisfied" e nunca copiados para dentro do venv — e desapareceram quando `PYTHONPATH` foi removido para rodar os scripts. Solução: `pip install --ignore-installed`.
6. `setuptools>=81` removeu o módulo `pkg_resources`, que pacotes mais antigos da árvore (`pyscenic`/`ctxcore`/`loompy`) ainda importam diretamente. Solução: fixar `setuptools<81`.
7. O `__init__.py` do próprio `scMalignantFinder` importa incondicionalmente seu submódulo `spatial` (usado só para identificação de região maligna em dados espaciais — feature que nunca usamos), que arrasta `squidpy → spatialdata → ome-zarr`, cuja cadeia de dependências exige `zarr>=3` enquanto outra parte da árvore trava em `zarr<3` — um conflito real do lado deles, numa parte do pacote que não usamos. Solução: patch local (só no nosso venv) tornando esse import opcional (`try/except ImportError`), já que `classifier.py` (a única parte usada) não depende de `spatial` nem de squidpy.

## 7. Arquivos gerados

Todos em `benchmarks/scmalignantfinder/`:

- `pretrained_model/` — modelo e lista de genes originais, baixados do Zenodo.
- `vangalen_counts.mtx` / `_genes.tsv` / `_barcodes.tsv` — contagens brutas do vanGalen exportadas do Seurat.
- `vangalen_prediction_refined.csv` — rótulo original do vanGalen (`PredictionRefined`), usado como padrão de comparação em todos os testes.
- `vangalen_scmalignantfinder_predictions.csv` — predição do modelo pré-treinado (teste #1).
- `train_{Naldini_2023,BeneytoCalabuig_2023,Lasry_2022}_{counts.mtx,genes.tsv,barcodes.tsv,labels.csv}` — contagens e rótulos de treino extraídos por dataset.
- `combine_training_aml.h5ad` — treino combinado (3 cohorts, 295.113 células).
- `combined_tumor_up_down_degs.txt` — lista de genes pan-câncer original (2.707 genes), usada nos testes #2–#6.
- `aml_specific_degs.txt` — lista de genes calculada nos nossos dados (1.415 genes), usada no teste #7.
- `vangalen_scmalignantfinder_predictions_retrained_*.csv` — predições de cada configuração de retreino testada.
- `vangalen_comparison_core.csv` — tabela unida (predição pré-treinada × rótulo original) usada na primeira comparação.
