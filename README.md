# HotZones: Variance-Guided Framework for Individual-Specific Epigenetic Networks

> **Paper**: *Variance-Guided HotZones Sharpen Patient Geometry in Individual-Specific Epigenetic Networks*
> **Authors**: Gaia Righetti, Kavya Singh, Giulia Yousefi, Federico Melograna, Fabio Stella, Kristel Van Steen
> **Contact**: kavya.singh@uliege.be

---

## Overview

HotZones is a computational framework for identifying compact, high-variance gene modules from Individual-Specific Networks (ISNs) built from epigenomic data. Rather than analysing entire networks edge-by-edge, HotZones focuses on specific network regions where inter-individual variation is highest and most coherent-the "hot zones" of biological heterogeneity.

Applied to DNA methylation profiles of 2,547 TCGA patients (with HotZone detection performed on a representative 50-sample subset), the framework identifies a 135-gene HotZone signature that yields:

- **2.3× improvement** in PCA variance capture (40.3% vs. 18.8%)
- **Higher-quality clustering** (silhouette 0.64 vs. 0.54) with more distinct patient subgroups (11 vs. 3)
- **1.3–1.4× superior pathway enrichment density** compared to random gene-set controls

---

## Repository Structure

```
HotZones/
├── README.md
├── 1_preprocessing/
│   ├── IMPUTATION.ipynb                # Missing value imputation (KNN)
│   └── COFOUNDER_CC.ipynb              # Confounding correction via polynomial regression
├── 2_mapping/
│   ├── mapping.ipynb                   # CpG probe-to-gene mapping
│   └── gene2DNAm_map.txt               # Reference annotation mapping file
├── 3_pca_summarization/
│   ├── PCA.ipynb                       # Per-gene PC1 summarisation across CpGs
│   └── DNAm_PCA_analysis_log.txt       # Run log: 447 genes, 2,546 samples
├── 4_isn_construction/
│   └── lioness.py                      # PANDA + LIONESS individual-specific network inference

── PipelineHot.ipynb               # HotZone detection + all validation analyses
│                                       # (PCA, t-SNE, clustering, heatmaps, random controls)
cohort_analysis/
└── Data-analysis.ipynb             # Subsample representativeness: PCA overlay and tumour stage comparison
```

---

## Pipeline Overview

The full pipeline runs in six sequential steps. Each step's output feeds directly into the next.

```
Raw TCGA DNAm data (450K array)
        │
        ▼
[Step 1a] KNN Imputation  ─────────────────────  IMPUTATION.ipynb
        │  → DNAm_processed_imputed.csv          (2,547 samples × 2,000 CpGs)
        ▼
[Step 1b] Confounding Correction  ─────────────  COFOUNDER_CC.ipynb
        │  → DNA_corrected.csv
        ▼
[Step 2] CpG-to-Gene Mapping  ─────────────────  mapping.ipynb
        │  → methylation_with_probeids.csv
        ▼
[Step 3] Per-Gene PCA Summarisation  ──────────  PCA.ipynb
        │  → DNAm_imputed_corr_1PC_filt_pcgene_python.csv  (447 genes × 2,546 individuals)
        ▼
[Step 4] ISN Construction: PANDA + LIONESS  ────  lioness.py
        │  → labeled_lioness_data.csv  (one network per patient)
        ▼
[Step 5] HotZone Detection + Full Validation  ─  PipelineHot.ipynb
        │  → filtered_df.csv (135-gene HotZone interactions)
        │  → PCA, t-SNE, clustering, heatmap, silhouette plots
        ▼
[Step 6] Cohort-Level Analysis  ───────────────  Data-analysis.ipynb
           → PCA overlay (full cohort vs. subsample, Supplementary Figure 36)
           → Tumour stage distribution comparison (Supplementary Figure 37)
```

---

## Step-by-Step Instructions

### Prerequisites

All notebooks were developed and tested in **Google Colab**. They can be run directly in Colab (recommended) or adapted to a local Jupyter environment.

**Dependencies — install all before running:**

```bash
# Core scientific stack
pip install pandas numpy scipy scikit-learn matplotlib seaborn

# Network inference (Step 4)
pip install netZooPy

# Graph algorithms and community detection (Step 5)
pip install leidenalg igraph networkx

# Interactive visualisation (Step 5)
pip install plotly
```

> **Note on Step 5 (`PipelineHot.ipynb`)**: The notebook's opening cells also run `pip install hdbscan umap-learn GraphRicciCurvature graphtools`. These packages are installed as setup but are **not imported or used** in the notebook's actual analysis code. You do not need them to reproduce the results, but they will be installed if you run the setup cells as written.

---

### Step 1a — Missing Value Imputation

**Notebook**: `1_preprocessing/IMPUTATION.ipynb`

**Purpose**: Handles missing CpG beta values in the raw TCGA methylation matrix using K-Nearest Neighbours imputation.

**Input**: `TCGA_DNAm_processed.csv`
- Rows: 2,547 samples (patients), no header
- Columns: 2,000 pre-selected CpG sites represented as methylation beta values

**What it does**:
- Loads the methylation matrix (confirmed shape from run log: 2,547 rows × 2,000 columns)
- Applies `KNNImputer` with k=15 neighbours from scikit-learn
- Reports missing value counts before and after imputation
- Saves output without index or header (`index=False, header=False`)

**Output**: `DNAm_processed_imputed.csv`
- Shape: 2,547 samples × 2,000 CpG sites
- Confirmed from run: 0 missing values before and after (data was already complete in this run)

**Key parameter**:
| Parameter | Value | Description |
|-----------|-------|-------------|
| `k_neighbors` | 15 | Number of neighbours for KNN imputation |

---

### Step 1b — Confounding Correction

**Notebook**: `1_preprocessing/COFOUNDER_CC.ipynb`

**Purpose**: Removes sample-order confounding from the imputed methylation matrix using polynomial regression residualisation based on sample index.

**Input**: `DNAm_processed_imputed.csv`

**What it does**:
- Creates polynomial covariates (degree 3) derived from sample index (`np.arange(n_samples)`)
- For each of the 2,000 columns, fits a `LinearRegression` model using the covariates
- Replaces observed values with residuals (observed − predicted), removing covariate-driven signal
- Prints progress every 1,000 columns

**Output**: `DNA_corrected.csv`
- Same shape as input (2,547 × 2,000); values are residuals centred near zero
- Saved with `index_label='TargetID'`

**Key parameter**:
| Parameter | Value | Description |
|-----------|-------|-------------|
| `degree` | 3 | Polynomial degree for covariate generation |

---

### Step 2 — CpG-to-Gene Mapping

**Notebook**: `2_mapping/mapping.ipynb`
**Output mapping file**: `2_mapping/gene2DNAm.map`

**Purpose**: Assigns authentic Illumina probe IDs to the imputed methylation matrix columns, then maps CpG probes to their target genes using RefGene annotations.

**Inputs**:
- `DNAm_processed_imputed.csv` — imputed methylation matrix (output of Step 1a)
- `humanmethylation450_15017482_v1-2.csv` — Illumina 450K array manifest (download separately; see Data Availability)

**What it does**:
- Reads the manifest and extracts valid probe IDs from the `IlmnID` column (skipping first 7 header rows)
- Assigns probe IDs to the columns of the methylation matrix
- Maps each probe to its gene symbol; handles multi-gene probes by splitting on `;`

**Outputs**:
- `methylation_with_probeids.csv` — methylation matrix with probe IDs as column names
- `gene2DNAm.map` — tab-separated file mapping each gene to its associated probe IDs

> **Note**: The Illumina manifest file is not included due to its size. Download it from (https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html) and place it in the `2_mapping/` folder before running.

---

### Step 3 — Per-Gene PCA Summarisation

**Notebook**: `3_pca_summarization/PCA.ipynb`
**Log**: `3_pca_summarization/DNAm_PCA_analysis_log.txt`

**Purpose**: Reduces multiple CpG measurements per gene to a single representative value per patient using the first principal component (PC1).

**Input**: `methylation_with_probeids.csv`

**What it does**:
- Groups CpG probes by gene symbol
- For each gene, runs PCA across all patients on its CpG sub-matrix
- Retains PC1 score as that gene's methylation value for each patient

**Output**: `DNAm_imputed_corr_1PC_filt_pcgene_python.csv`
- Shape: 2546 individuals × 447 genes
- Each entry = PC1 score of that gene in that individual

**Confirmed run statistics** (from `DNAm_PCA_analysis_log.txt`):
```
Total genes: 447
CpGs per gene — Min: 2, Max: 29, Median: 3, Mean: 3.41
Processing time: ~5 seconds
```
Note: The log reports matrix dimensions as [447, 2546] (genes × individuals), 
while the saved output is transposed to individuals × genes.

---

### Step 4 — Individual-Specific Network Construction

**Script**: `4_isn_construction/lioness.py`

**Purpose**: Constructs one co-methylation network per patient using PANDA followed by LIONESS. This is the core ISN inference step.

**Input**: `DNAm_imputed_corr_1PC_filt_pcgene_python.csv`

**What it does**:
1. Reads the patient × gene methylation matrix (index_col=0)
2. Computes a gene–gene Pearson correlation matrix directly from the patient × gene data.
3. Constructs motif data from this correlation matrix (positive correlation → edge weight 1)
4. Transposes the methylation matrix to gene × patient format for PANDA input
5. Runs **PANDA** (`netZooPy`) with these parameters:
   - `computing='cpu'`
   - `remove_missing=False`
   - `keep_expression_matrix=True`
   - `save_memory=False`
   - `modeProcess='legacy'`
6. Runs **LIONESS** on the PANDA result to derive per-patient edge weights via leave-one-out interpolation:

   `e_ij(s) = N · (e_ij(all) − e_ij(all−s)) + e_ij(all−s)`

7. Renames sample columns as `Sample_1`, `Sample_2`, ..., `Sample_N`
8. Exports the labeled result

**Output**: `labeled_lioness_data.csv`
- Columns: `tf`, `gene`, `Sample_1`, ..., `Sample_N`
- Each row = one gene–gene edge; each sample column = that edge's weight for that patient

**Dependencies**:
```bash
pip install netZooPy matplotlib pandas numpy
```

> **Note**: LIONESS is computationally intensive. At least 16 GB RAM is recommended. The script runs in CPU mode (`computing='cpu'`).

---

### Step 5 — HotZone Detection and Validation

**Notebook**: `5_hotzone_detection_and_validation/PipelineHot.ipynb`

**Purpose**: The central analysis notebook. Performs HotZone detection and all downstream validation: PCA, t-SNE, K-means clustering, heatmaps, silhouette analysis, and random gene-set controls.

**Input**: `labeled_lioness_data.csv`

---

#### Part A — HotZone Detection

**1. Edge variance and pre-filtering**
- Computes row-wise variance across all sample columns: `df.iloc[:,2:].var(axis=1)`
- Sorts edges by variance descending; retains the top 10%: `iloc[:int(len(df)*0.1)]`

**2. Progressive sparsification**
Quantile thresholds applied to the pre-filtered top-10% edges:

```python
quantile_levels = [0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.97, 0.98]
```

At each threshold, edges with variance below that quantile are removed.

**3. Leiden community detection**
At each threshold, a weighted undirected NetworkX graph is built and converted to igraph, then Leiden is applied:
```python
leidenalg.find_partition(G_ig, leidenalg.RBConfigurationVertexPartition, weights='weight')
random_seed = 42   # set via random.seed(42)
```
Resolution γ=1.0 (default for `RBConfigurationVertexPartition`). Gene-to-cluster assignments are stored at every threshold.

**4. Sankey diagram**
Tracks gene flow between communities across consecutive thresholds to visualise stability and fragmentation.

**5. Quality metrics per community**
- **Mean edge variance**: average `Row_Variance` of edges connected to genes in the community
- **Fiedler value (λ₂)**: computed via `scipy.linalg.eigvalsh` on the normalised Laplacian (`scipy.sparse.csgraph.laplacian` with `normed=True`); second smallest eigenvalue

**6. HotZone extraction**
At threshold Q98.0, genes from surviving communities are extracted based on top-variance edges (top ~2%). However, the saved filtered_df.csv contains all original edges whose TF and gene both belong to the HotZone gene set, not only the Q98-filtered edges.

```python
threshold_to_extract = 'Q98.0'   # communities 0, 1, 2
```
Results saved as `filtered_df.csv`. The three communities contain 57, 51, and 27 genes respectively, merged into the 135-gene HotZone set.

---

#### Part B — Validation Analyses

All analyses run on five gene configurations: All genes, HotZone genes (135), Community 0 (57), Community 1 (51), Community 2 (27).

**PCA**
- Data is used directly without explicit standardisation before PCA (no StandardScaler applied in the PCA steps)
- Projected to 2 principal components; variance explained reported per configuration

**t-SNE** (`sklearn.manifold.TSNE`)
- `n_components=2`, `perplexity=10`, `random_state=42`
- Applied to: All genes, HotZone genes, Community 0, Community 1, Community 2

**K-means clustering on PCA embeddings** (`sklearn.cluster.KMeans`, `random_state=42`)

| Gene set | n_clusters |
|----------|------------|
| HotZone genes (PCA) | 11 |
| All genes (PCA) | 3 |
| Community 0 (PCA) | 2 |
| Community 1 (PCA) | 6 |
| Community 2 (PCA) | 2 |

Metrics computed per configuration: silhouette score, Calinski-Harabasz score, Davies-Bouldin score.

**t-SNE + K-means** (`perplexity=10`, `random_state=42`)

| Gene set | n_clusters |
|----------|------------|
| HotZone genes (t-SNE) | 8 |
| All genes (t-SNE) | 2 |
| Communities 0, 1, 2 (t-SNE) | 2 each |

**Silhouette curves are computed with varying ranges depending on the analysis:**
- Some PCA-based evaluations use K_range = range(2, 20) (k = 2–19)
- Some comparison plots use max_k = 10 (k = 2–10)
- The final random-control comparison uses max_k = 20 (k = 2–20)

**Heatmaps**
- `StandardScaler` applied per matrix
- Z-score heatmaps plotted with `seaborn.heatmap`, `cmap='vlag'`, axis tick labels removed

**Sankey comparison plot**
Shows how the 50 individuals distribute across cluster assignments in: All genes → Hot genes → Community 0 → Community 1 → Community 2

**Random gene-set controls**
- `n_iterations = 10` independent random subsets, each of `n_genes = 135` genes
- Sampled from all gene names in the full dataset
- Silhouette scores computed up to `max_k=20` for each and overlaid on the comparison plot

---

**Outputs**:
- `filtered_df.csv`- interaction matrix containing all original gene–gene (TF–target) edges where both genes belong to the 135-gene HotZone set, across all samples
- Sankey plot: community evolution across variance thresholds
- PCA scatter plots for all five gene configurations
- t-SNE scatter plots (plain and with K-means overlay) for all configurations
- Silhouette score plots: individual curves (k=2–19) and combined comparison with random controls (k=2–20)
- Z-score heatmaps for All genes and HotZone genes; community-specific matrices are also generated for downstream analyses
- Sankey comparison plot of K-means cluster assignments across gene configurations

**Dependencies**:
```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn
pip install leidenalg igraph networkx plotly
```

---
### Step 6 — Cohort-Level Analysis

**Notebook**: `6_cohort_analysis/Data-analysis.ipynb`

**Purpose**: Verifies that the 50-individual analysis subsample used in Step 5 adequately represents the full TCGA cohort (n = 2,546) in terms of molecular variance structure and clinical metadata composition, as reported in Supplementary Figures 36 and 37 of the paper.

**Inputs**:
- `methylation_full.csv` — full cohort methylation matrix; must contain column `submitter_id`
- `methylation_sub.csv` — 50-sample subsample methylation matrix; must contain column `submitter_id`
- `metadata_full.csv` — patient clinical metadata; must contain `submitter_id` and a grouping column


**What it does**:
- Loads and aligns all three files by `submitter_id`
- Standardises methylation data using StandardScaler: fitted on the full cohort and applied to both full cohort and subsample to ensure comparability
- Computes total molecular variance and PC1 variance explained in both the full cohort and subsample to quantify variance preservation (Supplementary Table 1)
- **PCA overlay** (`n_components=2`): fits PCA on the full cohort; projects the subsample onto the same PCA space; confirms that the subsample spans the same molecular axes of variance as the full cohort (Supplementary Figure 36)
- **Tumour stage distribution comparison**: compares the proportion of tumour stages in the subsample against the full cohort to confirm the subsample does not introduce sampling bias with respect to clinical metadata (Supplementary Figure 37)

**Outputs**:
- PCA scatter plot: full cohort and subsample projected onto the same axes (Supplementary Figure 36)
- Bar plot: tumour stage proportions in full cohort vs. subsample (Supplementary Figure 37)

**Dependencies**:
```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn
```

---

## Data Availability

| Data | Source | Access |
|------|--------|--------|
| TCGA DNA methylation | The Cancer Genome Atlas | [GDC Portal](https://portal.gdc.cancer.gov/) |
| Illumina 450K array manifest (`humanmethylation450_15017482_v1-2.csv`, 192 MB) | Illumina Support | [Download CSV](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html) |
| HotZone gene list (135 genes) | This study | See Supplementary Table 4 |

---

## Citation

If you use this code or framework in your work, please cite:

> Righetti G., Singh K., Yousefi G., Melograna F., Stella F., Van Steen K. (2025). *Variance-Guided HotZones Sharpen Patient Geometry in Individual-Specific Epigenetic Networks.*

---

## License

This code is provided for academic and research use. Please contact the corresponding author for licensing terms.ks
