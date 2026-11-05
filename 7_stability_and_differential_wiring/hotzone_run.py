# =============================================================
# hotzone_run.py
# Runs HotZone detection on LIONESS output
# Extracts Q98 HotZone genes and saves them
# Usage: python hotzone_run.py --run_id 1
# =============================================================

import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg
import argparse
import random
import os

# -------------------------
# Parse run ID
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--run_id", type=int, required=True)
args = parser.parse_args()

RUN_ID     = args.run_id
RUN_DIR    = f"run_{RUN_ID}"
RANDOM_SEED = 42  # fixed for HotZone clustering reproducibility

random.seed(RANDOM_SEED)

print(f"[Run {RUN_ID}] Starting HotZone detection...")

# -------------------------
# Load LIONESS output
# -------------------------
lioness_path = os.path.join(RUN_DIR, "lioness_labeled.csv")
df = pd.read_csv(lioness_path)
print(f"[Run {RUN_ID}] LIONESS data shape: {df.shape}")

# -------------------------
# STEP 1: Compute row variance and keep top 10%
# -------------------------
df['Row_Variance'] = df.iloc[:, 2:].var(axis=1)
df_sorted = df.sort_values(by="Row_Variance", ascending=False)
top_10_percent = df_sorted.iloc[:int(len(df_sorted) * 0.1)]

print(f"[Run {RUN_ID}] Top 10% edges: {len(top_10_percent)}")

# -------------------------
# STEP 2: Define quantile thresholds
# -------------------------
quantile_levels = [0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.97, 0.98]
quantiles = [top_10_percent['Row_Variance'].quantile(q) for q in quantile_levels]
threshold_labels = [f"Q{int(q*1000)/10}" for q in quantile_levels]

# -------------------------
# STEP 3: Leiden clustering at each threshold
# -------------------------
gene_cluster_map = {}

for i, threshold in enumerate(quantiles):
    label = threshold_labels[i]
    filtered = top_10_percent[top_10_percent['Row_Variance'] >= threshold]

    # Build graph
    G_nx = nx.Graph()
    for _, row in filtered.iterrows():
        G_nx.add_edge(row['tf'], row['gene'], weight=row['Row_Variance'])

    edges = [(u, v, float(d['weight'])) for u, v, d in G_nx.edges(data=True)]
    G_ig = ig.Graph.TupleList(edges, edge_attrs=['weight'])

    # Leiden clustering
    partition = leidenalg.find_partition(
        G_ig,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight'
    )

    cluster_map = {}
    for cluster_id, cluster in enumerate(partition):
        for idx in cluster:
            gene = G_ig.vs[idx]['name']
            cluster_map[gene] = cluster_id

    gene_cluster_map[label] = cluster_map
    print(f"[Run {RUN_ID}] Threshold {label}: {len(set(cluster_map.values()))} clusters, {len(cluster_map)} genes")

# -------------------------
# STEP 3b: Compute mean edge-variance and Fiedler value per cluster
# (replicates Gaia's original Cells 31 and 33 exactly)
# -------------------------
from collections import defaultdict
from scipy.sparse.csgraph import laplacian
from scipy.linalg import eigvalsh

# --- Mean variance per cluster (Cell 31 logic) ---
cluster_variances = defaultdict(dict)

for threshold in threshold_labels:
    gene_cluster_map_for_threshold = gene_cluster_map[threshold]

    for cluster_id in set(gene_cluster_map_for_threshold.values()):
        genes_in_cluster = [gene for gene, assigned_cluster in gene_cluster_map_for_threshold.items() if assigned_cluster == cluster_id]

        mask = top_10_percent['tf'].isin(genes_in_cluster) | top_10_percent['gene'].isin(genes_in_cluster)
        cluster_variance_values = top_10_percent[mask]['Row_Variance']

        mean_variance = cluster_variance_values.mean() if not cluster_variance_values.empty else 0
        cluster_variances[threshold][cluster_id] = mean_variance

# --- Fiedler value per cluster (Cell 33 logic) ---
fiedler_values = defaultdict(dict)

for threshold in threshold_labels:
    gene_cluster_map_for_threshold = gene_cluster_map[threshold]
    filtered_df = top_10_percent[
        top_10_percent['Row_Variance'] >= quantiles[threshold_labels.index(threshold)]
    ]

    for cluster_id in set(gene_cluster_map_for_threshold.values()):
        cluster_genes = [gene for gene, cid in gene_cluster_map_for_threshold.items() if cid == cluster_id]

        cluster_edges = filtered_df[
            filtered_df['tf'].isin(cluster_genes) & filtered_df['gene'].isin(cluster_genes)
        ]

        if len(cluster_edges) < 2:
            fiedler_values[threshold][cluster_id] = 0
            continue

        G_sub = nx.Graph()
        for _, row in cluster_edges.iterrows():
            G_sub.add_edge(row['tf'], row['gene'], weight=row['Row_Variance'])

        A = nx.to_numpy_array(G_sub, weight='weight')
        L = laplacian(A, normed=True)
        eigenvalues = eigvalsh(L)
        fiedler = eigenvalues[1] if len(eigenvalues) > 1 else 0
        fiedler_values[threshold][cluster_id] = fiedler    

# -------------------------
# STEP 4: Extract Q98 HotZone genes
# -------------------------
threshold_to_extract = 'Q98.0'
hotzone_genes = list(gene_cluster_map[threshold_to_extract].keys())

print(f"[Run {RUN_ID}] HotZone genes at {threshold_to_extract}: {len(hotzone_genes)}")

# -------------------------
# STEP 5: Save HotZone genes
# -------------------------
genes_path = os.path.join(RUN_DIR, "hotzone_genes.csv")
pd.Series(hotzone_genes, name="gene").to_csv(genes_path, index=False)
print(f"[Run {RUN_ID}] HotZone genes saved to: {genes_path}")

# Also save full gene-cluster mapping at Q98
cluster_df = pd.DataFrame([
    {"gene": gene, "cluster": cluster_id}
    for gene, cluster_id in gene_cluster_map[threshold_to_extract].items()
])
cluster_path = os.path.join(RUN_DIR, "hotzone_gene_clusters.csv")
cluster_df.to_csv(cluster_path, index=False)
print(f"[Run {RUN_ID}] Gene-cluster mapping saved to: {cluster_path}")

# Save cluster-level variance and Fiedler stats at Q98 (for the 20-replicate stability analysis)
q98_label = 'Q98.0'
stats_rows = []
for cluster_id in set(gene_cluster_map[q98_label].values()):
    n_genes = list(gene_cluster_map[q98_label].values()).count(cluster_id)
    stats_rows.append({
        "cluster": cluster_id,
        "n_genes": n_genes,
        "mean_edge_variance": cluster_variances[q98_label][cluster_id],
        "fiedler_value": fiedler_values[q98_label][cluster_id]
    })

stats_df = pd.DataFrame(stats_rows)
stats_path = os.path.join(RUN_DIR, "hotzone_cluster_stats.csv")
stats_df.to_csv(stats_path, index=False)
print(f"[Run {RUN_ID}] Cluster variance/Fiedler stats saved to: {stats_path}")
