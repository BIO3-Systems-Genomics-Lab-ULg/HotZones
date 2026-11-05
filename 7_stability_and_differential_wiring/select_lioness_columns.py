import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--run_id", type=int, required=True)
args = parser.parse_args()

RUN_ID  = args.run_id
RUN_DIR = f"run_{RUN_ID}"

print(f"[Run {RUN_ID}] Selecting 50 columns from full-cohort LIONESS output...")

# Load the ONE precomputed full-cohort LIONESS output (2,546 sample columns)
full_lioness = pd.read_csv("lioness_labeled_full.csv")

# Load this replicate's selected 50 sample IDs (produced by subsample.py)
ids_path = os.path.join(RUN_DIR, "selected_ids.csv")
selected_ids = pd.read_csv(ids_path)["submitter_id"].tolist()

# Check all selected IDs actually exist as columns in the full LIONESS output
missing_ids = [sid for sid in selected_ids if sid not in full_lioness.columns]
if missing_ids:
    raise ValueError(
        f"[Run {RUN_ID}] {len(missing_ids)} selected sample IDs not found in "
        f"lioness_labeled_full.csv: {missing_ids[:5]}..."
    )

# Keep tf, gene, plus only the 50 selected sample columns
columns_to_keep = ["tf", "gene"] + selected_ids
labeled_data = full_lioness[columns_to_keep]

out_path = os.path.join(RUN_DIR, "lioness_labeled.csv")
labeled_data.to_csv(out_path, index=False)
print(f"[Run {RUN_ID}] Saved 50-column subset: {out_path}, shape: {labeled_data.shape}")
