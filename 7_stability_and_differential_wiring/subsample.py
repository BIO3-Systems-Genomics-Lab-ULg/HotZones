import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--run_id", type=int, required=True)
args = parser.parse_args()

RUN_ID  = args.run_id
SEED    = RUN_ID
N       = 50
IN_FILE = "DNAm_imputed_corr_1PC_filt_pcgene_python.csv"
OUT_DIR = f"run_{RUN_ID}"

os.makedirs(OUT_DIR, exist_ok=True)

print(f"[Run {RUN_ID}] Loading full methylation matrix...")
methyl = pd.read_csv(IN_FILE, index_col=0)
print(f"[Run {RUN_ID}] Full matrix shape: {methyl.shape}")

np.random.seed(SEED)
all_samples = methyl.index.tolist()
selected_samples = np.random.choice(all_samples, size=N, replace=False).tolist()

print(f"[Run {RUN_ID}] Selected {N} samples (seed={SEED})")
print(f"[Run {RUN_ID}] First 5 IDs: {selected_samples[:5]}")

methyl_sub = methyl.loc[selected_samples]
sub_path = os.path.join(OUT_DIR, "methylation_sub.csv")
methyl_sub.to_csv(sub_path)
print(f"[Run {RUN_ID}] Saved subsample: {sub_path} shape: {methyl_sub.shape}")

ids_path = os.path.join(OUT_DIR, "selected_ids.csv")
pd.Series(selected_samples, name="submitter_id").to_csv(ids_path, index=False)
print(f"[Run {RUN_ID}] Saved selected IDs: {ids_path}")
