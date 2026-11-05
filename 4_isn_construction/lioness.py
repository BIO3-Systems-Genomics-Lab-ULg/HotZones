import pandas as pd
import numpy as np
from netZooPy.lioness import Lioness
from netZooPy.lioness.analyze_lioness import AnalyzeLioness
from netZooPy.panda import Panda

# -------------------------
# This script runs ONCE, on the full 2,546-person cohort.
# It does NOT take a --run_id argument, since there is only
# one full-cohort LIONESS build (not one per replicate).
# -------------------------

IN_FILE  = "DNAm_imputed_corr_1PC_filt_pcgene_python.csv"
OUT_FILE = "lioness_labeled_full.csv"

print("Starting full-cohort LIONESS pipeline...")

exp_data = pd.read_csv(IN_FILE, index_col=0)
print(f"Full cohort shape: {exp_data.shape}")

# Sanity checks
assert exp_data.shape[0] == 2546, f"Expected 2546 samples, got {exp_data.shape[0]}"
assert exp_data.isnull().sum().sum() == 0, "Missing values detected"
print(f"Genes: {exp_data.shape[1]}, Samples: {exp_data.shape[0]}")
print(f"First 3 sample IDs: {exp_data.index[:3].tolist()}")

selected_ids = exp_data.index.tolist()

prior_network = exp_data.T
corr_matrix = exp_data.corr()
gene_names = corr_matrix.columns
pairs = pd.MultiIndex.from_product([gene_names, gene_names], names=["source", "target"])
motif_data = pd.DataFrame(index=pairs)

# original binary prior — kept as is
motif_data["weight"] = (corr_matrix.values.flatten() > 0).astype(int)
motif_data.reset_index(inplace=True)

print("Running PANDA...")
panda_obj = Panda(
    prior_network,
    motif_data,
    None,
    computing="cpu",
    remove_missing=False,
    keep_expression_matrix=True,
    save_memory=False,
    modeProcess="legacy"
)

print("Running LIONESS...")
lioness_obj = Lioness(panda_obj)

analyze_obj = AnalyzeLioness(lioness_obj)
labeled_data = analyze_obj.lioness_results.copy()
print(f"LIONESS output shape: {labeled_data.shape}")

# Preserve real TCGA sample IDs
labeled_data.columns = ["tf", "gene"] + selected_ids
print(f"First 5 columns: {labeled_data.columns[:5].tolist()}")

labeled_data.to_csv(OUT_FILE, index=False)
print(f"Full-cohort LIONESS output saved: {OUT_FILE}")
print("Done.")
