import pandas as pd
import numpy as np
import argparse
import os
from netZooPy.lioness import Lioness
from netZooPy.lioness.analyze_lioness import AnalyzeLioness
from netZooPy.panda import Panda

parser = argparse.ArgumentParser()
parser.add_argument("--run_id", type=int, required=True)
args = parser.parse_args()

RUN_ID  = args.run_id
RUN_DIR = f"run_{RUN_ID}"

# Change 1: ensure directory exists
os.makedirs(RUN_DIR, exist_ok=True)

print(f"[Run {RUN_ID}] Starting LIONESS pipeline...")

sub_path = os.path.join(RUN_DIR, "methylation_sub.csv")
exp_data = pd.read_csv(sub_path, index_col=0)
print(f"[Run {RUN_ID}] Subsample shape: {exp_data.shape}")

# Change 2: sanity checks
assert exp_data.shape[0] == 50, f"Expected 50 samples, got {exp_data.shape[0]}"
assert exp_data.isnull().sum().sum() == 0, "Missing values detected"
print(f"[Run {RUN_ID}] Genes: {exp_data.shape[1]}, Samples: {exp_data.shape[0]}")
print(f"[Run {RUN_ID}] First 3 sample IDs: {exp_data.index[:3].tolist()}")

selected_ids = exp_data.index.tolist()
ids_path = os.path.join(RUN_DIR, "selected_ids.csv")
pd.Series(selected_ids, name="submitter_id").to_csv(ids_path, index=False)
print(f"[Run {RUN_ID}] Saved {len(selected_ids)} sample IDs")

prior_network = exp_data.T
corr_matrix = exp_data.corr()
gene_names = corr_matrix.columns
pairs = pd.MultiIndex.from_product([gene_names, gene_names], names=["source", "target"])
motif_data = pd.DataFrame(index=pairs)

# your original binary prior — kept as you want
motif_data["weight"] = (corr_matrix.values.flatten() > 0).astype(int)
motif_data.reset_index(inplace=True)

print(f"[Run {RUN_ID}] Running PANDA...")
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

print(f"[Run {RUN_ID}] Running LIONESS...")
lioness_obj = Lioness(panda_obj)

analyze_obj = AnalyzeLioness(lioness_obj)
labeled_data = analyze_obj.lioness_results.copy()
print(f"[Run {RUN_ID}] LIONESS output shape: {labeled_data.shape}")

# Change 3: preserve real TCGA sample IDs
labeled_data.columns = ["tf", "gene"] + selected_ids
print(f"[Run {RUN_ID}] First 5 columns: {labeled_data.columns[:5].tolist()}")

out_path = "/massstorage/URT/GEN/BIO3/Kavya/run_1/lioness_labeled.csv"
labeled_data.to_csv(out_path, index=False)
print(f"[Run {RUN_ID}] LIONESS output saved: {out_path}")
print(f"[Run {RUN_ID}] Done.")
