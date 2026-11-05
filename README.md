# HotZones: Variance-Guided Framework for Individual-Specific Epigenetic Networks

> **Paper**: *Variance-Guided HotZones Sharpen Patient Geometry in Individual-Specific Epigenetic Networks*
> **Authors**: Kavya Singh*, Gaia Righetti*, Giulia Tremolada, Fabio Stella, Kristel Van Steen
> **Contact**: kavya.singh@uliege.be

---

## Overview

HotZones is a computational framework for identifying compact, high-variance gene modules from Individual-Specific Networks (ISNs) built from epigenomic data. Rather than analysing a whole network edge by edge, HotZones focuses on the regions where inter-individual variation is highest and most coherent: the "hot zones" of biological heterogeneity.

The code in this repository applies the framework to DNA methylation profiles of 2,546 TCGA patients:

1. CpG-level methylation is summarised to one value per gene and patient (PC1 across the gene's CpGs).
2. One ISN per patient is inferred with PANDA + LIONESS on the full cohort.
3. HotZone discovery is run on a 50-patient subsample: edges are ranked by their variance across patients, progressively sparsified, and clustered with Leiden community detection. Communities are described by mean edge variance and Fiedler value, and the Q = 0.98 level gives the HotZone gene set.
4. The HotZone set is validated against three core criteria: geometry-based relevance (PCA, K-means with silhouette, t-SNE, heatmaps), biological coherence (pathway enrichment against the complementary genes and a size-matched random control), and replication-based robustness (a stability analysis over 20 independent 50-patient subsamples). Two further analyses are exploratory and are reported in the Supplementary Information: differential wiring between the two largest HotZone-restricted patient clusters (Section S5.1) and clinical-endpoint associations (Section S5.2). A cohort-representativeness check of the 50-patient subsample is also included. For this proof-of-principle analysis, validity rests primarily on the three core criteria, not on mechanistic or clinical confirmation.

The repository follows the order of the analyses, so each step can be traced from raw data to its outputs. 

---

## Repository Structure

```
HotZones/
├── README.md
├── 1_preprocessing/
│   ├── imputation.ipynb                       # KNN imputation of CpG beta values (k = 15)
│   └── co_founder_correction.ipynb              # Diagnostic only: variance explained by a degree-3 polynomial of sample order
├── 2_mapping/
│   └── mapping.ipynb                          # Assigns Illumina 450K probe IDs to CpG columns; builds the gene-to-CpG map
├── 3_pca_summarization/
│   └── PCA.ipynb                              # Per-gene PC1 of the CpGs; restriction to HGNC protein-coding genes
├── 4_isn_construction/
│   └── lioness.py                             # Full-cohort PANDA + LIONESS (all 2,546 patients)
├── 5_hotzone_detection_and_geometry_validation/
│   └── FINALPipelineHot.ipynb                 # Variance graph, sparsification, Leiden, Sankey, mean variance, Fiedler;
│                                               # PCA / K-means / silhouette / t-SNE / heatmaps on the HotZone set
├── 6_pathway_enrichment/
│   └── FinalChapter.ipynb                     # Hypergeometric pathway enrichment, cumulative resampling curves, Fisher variant
├── 7_stability_and_differential_wiring/
│   ├── subsample.py                           # Replicate step 1: draw a seeded 50-patient subsample
│   ├── select_lioness_columns.py              # Replicate step 2: pick that replicate's 50 columns from the full LIONESS output
│   ├── hotzone_run.py                         # Replicate step 3: variance graph + Leiden + cluster statistics at Q = 0.98
│   ├── consistency_analysis.py                # Replicate step 4: cross-replicate Jaccard similarity and core-gene extraction
│   ├── HotZone_Stability_Analysis_KVS.ipynb   # R notebook: stability analysis (core criterion) AND exploratory differential wiring
│   └── cluster_job_submission_notes.docx      # HPC/cluster job submission order for the 20-replicate array job
└── 8_clinical_and_cohort_validation/
    └── Data_analysis-NPJ.ipynb                # Cohort representativeness, sample-ID mapping, exploratory clinical-endpoint tests
```

---

## Pipeline Overview

```
Raw TCGA DNAm data (450K array; 2,000 pre-selected CpG columns)
        │
        ▼
[Step 1a] KNN imputation (k = 15),             2,546 individuals  ──  1_preprocessing/imputation.ipynb
        │  → DNAm_processed_imputed.csv        (2,546 x 2,000)
        │
        ├─► [Step 1b] Polynomial sample-order diagnostic (not used downstream)  ──  1_preprocessing/co_founder_correction.ipynb
        │       → DNA_corrected.csv
        ▼
[Step 2] Probe-ID assignment + gene-to-CpG mapping  ──  2_mapping/mapping.ipynb
        │  → methylation_with_probeids.csv, gene2DNAm.map
        ▼
[Step 3] Per-gene PC1 + HGNC protein-coding filter  ──  3_pca_summarization/PCA.ipynb
        │  → DNAm_PCA_results.csv (all genes with >= 2 CpGs)
        │  → DNAm_imputed_corr_1PC_filt_pcgene_python.csv  (2,546 patients x 359 genes)
        ▼
[Step 4] Full-cohort ISNs: PANDA + LIONESS  ──  4_isn_construction/lioness.py
        │  → lioness_labeled_full.csv   (tf, gene, one column per patient, real TCGA IDs)
        │
        ├────────────────────────────────────────────┐
        ▼                                            ▼
[Step 5] Main discovery subsample (50 patients)     [Step 7] 20 replicate subsamples
   labeled_lioness_data.csv (Sample_1..Sample_50)      subsample.py → select_lioness_columns.py →
        │                                              hotzone_run.py (per run_N)
        ▼                                              → consistency_analysis.py
   FINALPipelineHot.ipynb                              → HotZone_Stability_Analysis_KVS.ipynb (R)
   → filtered_df.csv + figures
        │
        ▼
[Step 6] FinalChapter.ipynb
   → enrichment-curve figures
        │
        ▼
[Step 8] Data_analysis-NPJ.ipynb: cohort representativeness,
         HotZone-cluster vs clinical-endpoint tests
```

---

## Step-by-Step Instructions

### Prerequisites

Notebooks were developed and run in **Google Colab**; the scripts in `7_stability_and_differential_wiring/` were run on an HPC/Slurm cluster (see `cluster_job_submission_notes.docx`). All can be adapted to a local Jupyter/R environment. Colab-specific calls (`google.colab.files`, `drive.mount`) must be replaced by local file paths when running locally.

```bash
# Core scientific stack
pip install pandas numpy scipy scikit-learn matplotlib seaborn statsmodels scikit-bio

# Network inference (Step 4)
pip install netZooPy

# Graph algorithms and community detection (Steps 5 and 7)
pip install leidenalg igraph networkx

# Interactive visualisation (Step 5)
pip install plotly

# Imported by the pathway notebook (Step 6); the enrichment tests themselves use scipy
pip install gseapy
```

The R notebook in Step 7 installs its own packages on first run: `ggplot2`, `vegan`, `circlize`, `ComplexHeatmap`, `biomaRt`, `UpSetR`, `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA` (CRAN/Bioconductor, via `BiocManager`). In Colab, set *Runtime > Change runtime type > R* before running it.

---

### Step 1a: Missing-value imputation

**Notebook**: `1_preprocessing/imputation.ipynb`

**Input**: `TCGA_DNAm_processed.csv`: raw samples × 2,000 pre-selected CpG columns (methylation beta values).

**What it does**:
1. Uses the 2,546 individuals analysed in all downstream analyses
2. Asserts the resulting cohort size.
3. Applies scikit-learn `KNNImputer` (k = 15) and reports the numbers of missing values before and after. Observed values are left unchanged; if the input has no missing values, the notebook says so and imputation changes nothing.

**Output**: `DNAm_processed_imputed.csv` (2,546 × 2,000). 

---

### Step 1b: Sample-order diagnostic (not used downstream)

**Notebook**: `1_preprocessing/co_founder_correction.ipynb`

**Input**: `DNAm_processed_imputed.csv`

**What it does**: For each CpG column, fits a degree-3 polynomial of the sample order (row position) by linear regression, subtracts the fit, and reports how much variance the polynomial explains (mean, median and maximum R² across CpGs).

**Output**: `DNA_corrected.csv` (the residuals). 

---

### Step 2: Probe-ID assignment and gene-to-CpG mapping

**Notebook**: `2_mapping/mapping.ipynb` (two cells, run in order)

**Inputs**:
- `DNAm_processed_imputed.csv` (output of Step 1a; shape is asserted to be 2,546 × 2,000)
- `humanmethylation450_15017482_v1-2.csv`: Illumina 450K manifest (download separately; see Data Availability; not included because of its large size)

**Cell 1, probe-ID assignment**: reads the manifest (`skiprows=7`) and takes the unique `IlmnID` values in manifest order. Output: `methylation_with_probeids.csv` (index kept).

**Cell 2, gene mapping**: maps every probe to its gene symbols through the manifest's `UCSC_RefGene_Name` field; probes with several genes (separated by `;`) are added to each gene (duplicates removed), and probes without a gene are skipped. Output: `gene2DNAm.map` (tab-separated: gene, then its probe IDs).

---

### Step 3: Per-gene PCA summarisation

**Notebook**: `3_pca_summarization/PCA.ipynb`

**Inputs**: `methylation_with_probeids.csv` (read with `index_col=0`; 2,546 samples asserted) and `gene2DNAm.map`.

**What it does**:
1. Keeps genes with **at least 2 CpGs** present in the matrix.
2. For each gene (parallelised with `multiprocessing.Pool(10)`): standardises the gene's CpGs (`StandardScaler`), runs PCA (`svd_solver='full'`, `random_state=1`) across all patients, and keeps **PC1** as the gene's representative methylation value per patient. The variance explained by PC1 and the number of PCs needed for 90% of the variance are also computed per gene.
3. Restricts the genes to **HGNC approved protein-coding genes** (`hgnc_complete_set.txt`, downloaded automatically from the Google Storage or EBI mirror; if both fail, the notebook asks for a manual upload). The final matrix is asserted to be 2,546 × 359.

**Outputs**: `DNAm_PCA_results.csv` (PC1 of all genes with ≥ 2 CpGs), `DNAm_imputed_corr_1PC_filt_pcgene_python.csv` (2,546 patients × 359 genes; used by all later steps) and `DNAm_PCA_analysis.log`.

---

### Step 4: Full-cohort individual-specific network construction

**Script**: `4_isn_construction/lioness.py`

Runs **once**, on the complete 2,546-patient cohort.

**Input**: `DNAm_imputed_corr_1PC_filt_pcgene_python.csv`

**What it does**:

1. Computes the gene–gene Pearson correlation matrix from the patient × gene matrix.
2. Builds a binary prior scaffold: every gene pair with positive correlation gets weight 1, otherwise 0.
3. Runs **PANDA** (`netZooPy`; `computing="cpu"`, `remove_missing=False`, `keep_expression_matrix=True`, `save_memory=False`, `modeProcess="legacy"`).
4. Runs **LIONESS** on the PANDA object to obtain per-patient edge weights.
5. Renames the columns to `tf`, `gene` followed by the real TCGA submitter IDs.

**Output**: `lioness_labeled_full.csv`: one row per gene–gene edge (all ordered pairs of the 359 genes), columns `tf`, `gene`, plus one column per patient.

**This single output feeds two branches**:
- the **main discovery analysis** (Step 5) through `labeled_lioness_data.csv`, which contains the 50 columns of `lioness_labeled_full.csv` belonging to the discovery subsample (IDs in `methylation_sub.csv`).
- the **20-replicate stability analysis** (Step 7) through `select_lioness_columns.py`.


**Dependencies**: `netZooPy`, `pandas`, `numpy`. LIONESS is computationally intensive (all 359² gene pairs × 2,546 patients); ≥ 16 GB RAM is recommended; it runs in CPU mode.

---

### Step 5: HotZone detection and geometry-based validation

**Notebook**: `5_hotzone_detection_and_geometry_validation/FINALPipelineHot.ipynb`

**Input**: `labeled_lioness_data.csv`: columns `tf`, `gene`, `Sample_1 … Sample_50` (the main 50-patient discovery subsample; see the note in Step 4).

#### Part A: HotZone detection

1. **Edge variance and pre-filtering**: row-wise variance across the sample columns (`Row_Variance`); edges sorted in decreasing order and the **top 10%** kept.
2. **Progressive sparsification**: quantile thresholds `[0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.97, 0.98]`, computed **within the pre-filtered top 10%** of edges (labelled `Q0.0 … Q98.0`).
3. **Leiden community detection at each threshold**: a weighted undirected graph is built from the edges at or above the threshold (edge weight = `Row_Variance`) and partitioned with `leidenalg.find_partition(..., RBConfigurationVertexPartition, weights='weight')` (default resolution parameter 1.0). `random.seed(42)` is set at the start; no explicit `seed` argument is passed to `find_partition` (see Reproducibility notes). Leiden is therefore not reproducible run to run, and community sizes can differ slightly between runs.
4. **Sankey plot** (Plotly) tracking how genes move between clusters across thresholds.
5. **Quality metrics per community and threshold**:
   - *mean edge variance*: mean `Row_Variance` of the top-10% edges that have at least one endpoint in the community;
   - *Fiedler value*: second-smallest eigenvalue of the normalised Laplacian (`scipy`) of the weighted subgraph formed by the community's own edges at that threshold (set to 0 if the community has fewer than 2 edges).
6. **HotZone extraction at Q = 0.98**: all genes assigned to a community at `Q98.0` form the HotZone gene set; `filtered_df` keeps every edge of the full table whose two endpoints are both in that set.

**Output written by the notebook**: `filtered_df.csv`.

#### Part B: Geometry-based validation

PCA is run directly on the edge values (`PCA(n_components=2)`, no standardisation). The same analysis is repeated for five feature sets:

| Feature set | Edges used |
|---|---|
| All genes | every edge of the full table |
| HotZone set | `filtered_df` (edges inside the pooled HotZone gene set) |
| Community 0, 1, 2 | edges inside each Q = 0.98 community separately |

For each feature set the notebook produces: the 2-D PCA scatter; K-means (`random_state=42`) with the number of clusters fixed in the cell after a **silhouette sweep over k = 2 … 19**; and the silhouette, Calinski–Harabasz and Davies–Bouldin scores of the chosen solution. Further analyses:

- Sankey plot linking the K-means cluster labels of the five feature sets;
- Z-score heatmaps (`StandardScaler`) of the individual × edge matrices;
- t-SNE (`perplexity=10`, `random_state=42`) followed by K-means.

> The chosen K-means *k* for each feature set is typed into the notebook after inspecting the silhouette sweep. The notebook also compares the HotZone set with 10 random sets of 135 genes drawn from all 359 genes (silhouette sweep, k = 2 to 20, raw edge space). The size-matched pathway comparison is in Step 6.

**Dependencies**: `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `leidenalg`, `igraph`, `networkx`, `plotly`.

---

### Step 6: Pathway enrichment (biological validation)

**Notebook**: `6_pathway_enrichment/FinalChapter.ipynb`

**Inputs**:
- `labeled_lioness_data.csv`: the union of the `tf` and `gene` columns defines the **background gene set**.
- The **discovery (HotZone) gene list**, which is typed directly into the notebook (cell "discovery_genes"). The **complementary set** is the background minus the discovery genes. 

| Database | File |
|---|---|
| Reactome (MSigDB C2:CP) | `c2.cp.reactome.v2025.1.Hs.symbols.gmt.txt` |
| WikiPathways curated (MSigDB C2:CP) | `c2.cp.wikipathways.v2025.1.Hs.symbols.gmt.txt` |
| MSigDB canonical pathways (C2:CP) | `c2.cp.v2025.1.Hs.symbols.gmt.txt` |

**What it does**:
1. **Hypergeometric test** (`scipy.stats.hypergeom`) of the overlap between a gene subset and each pathway, with the background as population. A pathway counts as significant at **p < 0.05 (unadjusted)**.
2. **Metric**: the number of genes of the tested subset that belong to at least one significant pathway.
3. **Cumulative resampling curves**: random subsets of 10%, 20%, … 100% of the gene list (`steps=10`; 50 random draws per fraction), compared for *discovery vs complementary* and *discovery vs a size-matched random control* drawn from the complementary genes. Curves show the mean ± 1 SD, and ratio curves (discovery/complementary) are also computed. Python's `random` module is seeded with 42 inside the curve functions and with 10, 12 or 20 before the size-matched draws, depending on the cell.
4. **Fisher variant**: for the MSigDB canonical collection (C2:CP), the cumulative curves are repeated with a one-sided Fisher exact test (`alternative='greater'`) instead of the hypergeometric test.

**Outputs**: enrichment-curve figures for each collection.

**Dependencies**: `pandas`, `numpy`, `scipy`, `matplotlib` (`gseapy` is imported by the notebook but not used for the tests).

---

### Step 7: Stability analysis and differential wiring

**Folder**: `7_stability_and_differential_wiring/`

#### 7a. Python replicate pipeline (20 independent 50-patient subsamples)

Run from a directory that contains `DNAm_imputed_corr_1PC_filt_pcgene_python.csv` (Step 3) and `lioness_labeled_full.csv` (Step 4):

```bash
for i in $(seq 1 20); do
  python subsample.py            --run_id $i
  python select_lioness_columns.py --run_id $i
  python hotzone_run.py          --run_id $i
done
python consistency_analysis.py --n_runs 20
```

(On the HPC cluster each `run_id` was submitted as an array job; see `cluster_job_submission_notes.docx`.)

| Script | What it does | Writes |
|---|---|---|
| `subsample.py --run_id N` | Draws 50 patients without replacement from the gene-level PC1 matrix with `np.random.seed(run_id)` | `run_N/methylation_sub.csv`, `run_N/selected_ids.csv` (column `submitter_id`) |
| `select_lioness_columns.py --run_id N` | Keeps `tf`, `gene` and the 50 selected patient columns of `lioness_labeled_full.csv`| `run_N/lioness_labeled.csv` |
| `hotzone_run.py --run_id N` | Same variance-graph + Leiden procedure as Step 5 (top 10% edges, the same eight quantile thresholds, mean variance and Fiedler value per community); extracts the `Q98.0` HotZone genes. `random.seed(42)` is set but does not fix Leiden | `run_N/hotzone_genes.csv`, `run_N/hotzone_gene_clusters.csv`, `run_N/hotzone_cluster_stats.csv` (cluster, n_genes, mean_edge_variance, fiedler_value) |
| `consistency_analysis.py --n_runs 20` | Reads all `run_N/hotzone_genes.csv`; pairwise Jaccard similarity; per-gene presence frequency; core genes present in **≥ 80%** of the runs | `jaccard_similarity_matrix.csv`, `jaccard_heatmap.png`, `gene_frequency_across_runs.csv`, `gene_consistency_barplot.png`, `core_stable_genes.csv` |

#### 7b. R notebook: `HotZone_Stability_Analysis_KVS.ipynb`

The notebook implements the stability analysis (replication-based robustness, a core validation criterion), followed by the differential-wiring analysis. It runs in R (Colab: *Runtime > R*).

**Stability analysis**:
- bar chart of cluster sizes per replicate;
- gene × run binary presence matrix and consistency score; clustered heatmap with a "runs present" annotation; UpSetR overlap plot;
- histogram of gene consistency with the 80% threshold (`ceiling(0.8 × n_runs)`);
- Jaccard dissimilarity MDS and mean/SD of pairwise Jaccard similarity;
- extraction of genes present in ≥ 80% of runs (`highly_consistent_genes_80pct.csv`) and their annotation via `biomaRt` (Ensembl; synonym fallback; needs internet) (`highly_consistent_fallback_annotated_genes.csv`);
- Simple Matching Coefficient (SMC) and Matthews Correlation Coefficient (MCC) matrices over the 358 genes selected in at least one run, their MDS plots, and the mean/SD of pairwise SMC and MCC;
- Ward hierarchical clustering on 1 − SMC with an outlier cut-off (mean + 2.5 SD of the dendrogram heights);
- variance–cohesion landscape: mean edge variance versus Fiedler value of every cluster of every run, with 2-D density contours and cluster size as point size.

**Differential-wiring analysis (exploratory; Supplementary Section S5.1)** (same notebook, section "Differential analysis"). Inputs, expected in the working directory:

| File | Content | Produced by |
|---|---|---|
| `edge_matrix.csv` | edges × patients matrix; row names are `gene1__gene2` | **not generated by a notebook in this repository** (see note below) |
| `edges_df.csv` | edge list with columns `gene1`, `gene2` | same |
| `hotzone_genes.csv` | HotZone gene list of the main analysis, column `gene` | same |
| `cluster_labels_hot_11.csv` | patient-to-cluster labels (first column = patient, column `cluster`) of the HotZone K-means solution | same |
| `master_pathway_significant_genes.csv` | pathway-validated genes, column `gene` | 

The analysis (1) restricts the edges to those with both genes in the HotZone set; (2) averages edge weights per patient cluster; (3) flags clusters with fewer than 2 patients as singletons (compared descriptively with all other patients) and selects the two largest remaining clusters; (4) tests each edge between these two clusters with a Wilcoxon rank-sum test (mean difference), an F-test of variances (var.test), and a permutation test of the log2 variance ratio (B = 1000; no seed is set in the notebook), all with Benjamini–Hochberg FDR; (5) runs GO (biological process) and Reactome enrichment (clusterProfiler, ReactomePA) on the genes of the 100 edges with the largest absolute log2 variance ratio (a ranking, not an FDR-filtered set), using no p-value cut-off; the same is done separately for the top 100 edges of the permutation test, but only the variance-based enrichment tables are saved to disk; and (6) repeats the three edge tests restricted to edges with both genes in master_pathway_significant_genes.csv (within the HotZone-restricted edges, same two clusters).

**Outputs**: `res_mean.csv`, `res_var.csv`, `res_perm.csv`, `singleton_<cluster>_vs_rest.csv`, `go_enrichment_results.csv`, `reactome_enrichment_results.csv`, `ego_object.rds`, `react_object.rds`, `tier3_res_mean.csv`, `tier3_res_var.csv`, `tier3_res_perm.csv`, `tier3_summary.csv`, `tier3_significant_mean.csv`, `tier3_significant_var.csv`, `tier3_significant_perm.csv`.

> **Note on inputs**: `FINALPipelineHot.ipynb` writes only `filtered_df.csv`. The four files `edge_matrix.csv`, `edges_df.csv`, `hotzone_genes.csv` and `cluster_labels_hot_11.csv` listed above must be exported from the objects of Step 5 before this part is run; the export code is not part of the deposited notebook. Note that `hotzone_genes.csv` here is the main-analysis list, not the per-replicate `run_N/hotzone_genes.csv`.
>
> **Note on cell order**: the Jaccard "manuscript text" cell uses `n_genes`, which is defined later in the SMC/MCC cell, and the first `enrichGO` cell uses `top_genes`, which is defined in the following cell. Run the SMC/MCC cell before the Jaccard text cell, and the `top_genes` definition before the first enrichment call.

**Dependencies (Python)**: `pandas`, `numpy`, `networkx`, `igraph`, `leidenalg`, `matplotlib`, `seaborn`.
**Dependencies (R)**: `ggplot2`, `vegan`, `circlize`, `ComplexHeatmap`, `biomaRt`, `UpSetR`, `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`.

---

### Step 8: Cohort representativeness and exploratory clinical-endpoint analysis

**Notebook**: `8_clinical_and_cohort_validation/Data_analysis-NPJ.ipynb` (Colab; mounts Google Drive in the first cell)

The notebook contains two blocks of analyses that share its set-up cells. A helper cell converts an uploaded `.xlsx` file to `.csv`.

**Block 1: is the 50-patient subsample representative of the full cohort?**
- **Inputs**: `methylation_full.csv` (2,546 patients × 359 genes), `methylation_sub.csv` (the 50 discovery patients × 359 genes) and `metadata_full.csv` (tumour stage, age, race, gender, cancer type).
- **Variance and PCA**: both matrices are standardised with a scaler fitted on the full cohort. The notebook reports how much of the total variance the subsample keeps (102.1%), the PC1 variance share in each set (45.2% full cohort vs 46.9% subsample), and a PCA plot with the subsample projected onto the full-cohort axes.
- **Tumour stage composition**: stage counts and percentages in the full cohort and the subsample are tabulated and plotted (`cancer_type_distribution_full_vs_sub.csv`; despite the file name, this table contains tumour stage).

**Block 2: HotZone-cluster versus clinical-endpoint tests (exploratory; Supplementary Section S5.2)** 
- **Input**: `hotzone_clusters_with_metadata.csv`: HotZone-restricted cluster assignments joined to clinical metadata (columns `cluster`, `tumor_stage` coded 1–4, `test_cancers` coded 1–6, `gender` coded 0 = female / 1 = male). This file is assembled outside the notebook.
- **Tests**: Kruskal–Wallis for tumor stage across clusters (median stage per cluster saved); chi-square test of cancer type versus cluster, with a check of the expected-count assumption; global chi-square plus pairwise Fisher exact tests with Benjamini–Hochberg FDR for gender versus cluster. Each test is shown as a stacked bar plot and a count heatmap.
- **Outputs**: `median_stage_per_cluster.csv`, `tumor_stage_per_cluster.png`, `cancer_type_per_cluster.png`, `gender_per_cluster.png`, `fisher_gender_pairwise_fdr.csv`, `statistical_summary_final.csv` (the gender row of this summary reports the global chi-square statistic and p-value).

**Dependencies**: `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `statsmodels`, `scikit-bio`.

---

## Reproducibility notes

- **Leiden is not seeded.** Leiden clustering (Steps 5 and 7) is not seeded; `random.seed(42)` does not control `leidenalg`, so community sizes can differ slightly between runs.
- **Seeds.** K-means and t-SNE use `random_state=42`; the per-gene PCA uses `random_state=1`; replicate subsampling uses `np.random.seed(run_id)`. The permutation test in the R notebook is not seeded.
- **Probe-ID assignment** in Step 2 is positional (first 2,000 unique manifest IDs).
- **Pathway significance** in Step 6 uses unadjusted p < 0.05.
- **Package versions** are not pinned; record them (`pip freeze`, `sessionInfo()`) when archiving a release.

---

## Data Availability

## Data Availability

TCGA DNA methylation data were obtained via the NCI Genomic Data Commons portal (https://portal.gdc.cancer.gov/).

| Data | Source | Access |
|------|--------|--------|
| Illumina 450K array manifest (`humanmethylation450_15017482_v1-2.csv`, ~192 MB) | Illumina Support | [Download](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html) |
| HGNC gene table (`hgnc_complete_set.txt`) | HUGO Gene Nomenclature Committee | Downloaded automatically by `PCA.ipynb` |
| Pathway gene sets (Reactome, WikiPathways, MSigDB C2:CP v2025.1) | MSigDB | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) |
| HotZone gene list | This study | Supplementary Table; also typed into `FinalChapter.ipynb` |

---

## Citation

If you use this code or framework in your work, please cite:

> > **Authors**: Kavya Singh*, Gaia Righetti*, Giulia Tremolada, Fabio Stella, Kristel Van Steen. *Variance-Guided HotZones Sharpen Patient Geometry in Individual-Specific Epigenetic Networks.* (submitted).

---

## License

This code is provided for academic and research use. Please contact the corresponding author for licensing terms.
