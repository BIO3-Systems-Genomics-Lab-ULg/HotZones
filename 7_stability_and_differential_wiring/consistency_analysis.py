import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--n_runs", type=int, default=20)
args = parser.parse_args()

N_RUNS = args.n_runs
print(f"Consistency analysis across {N_RUNS} runs...")

gene_sets = {}
for run_id in range(1, N_RUNS + 1):
    genes_path = f"run_{run_id}/hotzone_genes.csv"
    if os.path.exists(genes_path):
        genes = pd.read_csv(genes_path)["gene"].tolist()
        gene_sets[run_id] = set(genes)
        print(f"Run {run_id}: {len(genes)} HotZone genes")
    else:
        print(f"Run {run_id}: FILE NOT FOUND - skipping")

available_runs = list(gene_sets.keys())
print(f"Available runs: {len(available_runs)}/{N_RUNS}")

n = len(available_runs)
jaccard_matrix = np.zeros((n, n))
for i, run_i in enumerate(available_runs):
    for j, run_j in enumerate(available_runs):
        set_i = gene_sets[run_i]
        set_j = gene_sets[run_j]
        intersection = len(set_i & set_j)
        union = len(set_i | set_j)
        jaccard_matrix[i, j] = intersection / union if union > 0 else 0.0

jaccard_df = pd.DataFrame(
    jaccard_matrix,
    index=[f"Run{r}" for r in available_runs],
    columns=[f"Run{r}" for r in available_runs]
)
print("Jaccard Similarity Matrix:")
print(jaccard_df.round(3).to_string())
jaccard_df.to_csv("jaccard_similarity_matrix.csv")
print("Saved: jaccard_similarity_matrix.csv")

plt.figure(figsize=(12, 10))
sns.heatmap(jaccard_df, annot=True, fmt=".2f", cmap="YlOrRd",
            vmin=0, vmax=1, linewidths=0.5, annot_kws={"size": 8})
plt.title("Jaccard Similarity of HotZone Genes Across 20 Runs", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("jaccard_heatmap.png", dpi=300, bbox_inches="tight")
print("Saved: jaccard_heatmap.png")

all_genes = set()
for genes in gene_sets.values():
    all_genes.update(genes)

gene_freq = {}
for gene in all_genes:
    count = sum(1 for genes in gene_sets.values() if gene in genes)
    gene_freq[gene] = count

freq_df = pd.DataFrame([
    {"gene": gene, "n_runs": count, "pct_runs": round(100 * count / len(available_runs), 1)}
    for gene, count in gene_freq.items()
]).sort_values("n_runs", ascending=False)

print(f"Total unique genes across all runs: {len(freq_df)}")
print("Top 20 most consistent genes:")
print(freq_df.head(20).to_string(index=False))
freq_df.to_csv("gene_frequency_across_runs.csv", index=False)
print("Saved: gene_frequency_across_runs.csv")

threshold_pct = 80
core_genes = freq_df[freq_df["pct_runs"] >= threshold_pct]
print(f"Core stable genes (>={threshold_pct}% of runs): {len(core_genes)}")
print(core_genes.to_string(index=False))
core_genes.to_csv("core_stable_genes.csv", index=False)
print("Saved: core_stable_genes.csv")

plot_df = freq_df.head(30)
colors = ["#2ecc71" if p >= threshold_pct else "#3498db" for p in plot_df["pct_runs"]]
plt.figure(figsize=(14, 6))
plt.bar(plot_df["gene"], plot_df["pct_runs"], color=colors, edgecolor="white")
plt.axhline(y=threshold_pct, color="red", linestyle="--", linewidth=1.5,
            label=f"{threshold_pct}% threshold (core stable)")
plt.xlabel("Gene", fontsize=12)
plt.ylabel("% Runs Gene Appears In", fontsize=12)
plt.title("HotZone Gene Consistency Across 20 Subsampling Runs", fontsize=14, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.legend()
plt.tight_layout()
plt.savefig("gene_consistency_barplot.png", dpi=300, bbox_inches="tight")
print("Saved: gene_consistency_barplot.png")

mean_jaccard = jaccard_df.values[np.triu_indices(n, k=1)].mean()
std_jaccard  = jaccard_df.values[np.triu_indices(n, k=1)].std()
print("===== SUMMARY =====")
print(f"Total runs completed:           {len(available_runs)}/{N_RUNS}")
print(f"Mean Jaccard similarity:        {mean_jaccard:.3f}")
print(f"Std Jaccard similarity:         {std_jaccard:.3f}")
print(f"Total unique genes detected:    {len(freq_df)}")
print(f"Core stable genes (>={threshold_pct}%):      {len(core_genes)}")

