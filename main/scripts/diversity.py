import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from qiime2 import Artifact
from skbio.stats.ordination import OrdinationResults


# ---- Load alpha diversity artifacts ----
shannon = qiime2.Artifact.load("core-metrics-results/shannon_vector.qza").view(pd.Series)
evenness = qiime2.Artifact.load("core-metrics-results/evenness_vector.qza").view(pd.Series)
chao1 = qiime2.Artifact.load("core-metrics-results/chao1_vector.qza").view(pd.Series)
simpson = qiime2.Artifact.load("core-metrics-results/simpson_vector.qza").view(pd.Series)

# ---- Load metadata ----
metadata = pd.read_csv("metadata.tsv", sep="\t", index_col=0)

# ---- Combine alpha diversity into one DataFrame ----
alpha_df = pd.concat([
    shannon.rename("Shannon"),
    evenness.rename("Evenness"),
    chao1.rename("Chao1"),
    simpson.rename("Simpson")
], axis=1)

# ---- Merge alpha diversity with metadata ----
merged = alpha_df.join(metadata)
merged['Group'] = merged['Group'].str.strip()

# ---- Output folder ----
os.makedirs("alpha_diversity_plots", exist_ok=True)
sns.set(style="whitegrid")

# ---- Define comparisons in flexible style ----
comparisons = [
    ("Group", None),                # all groups
    ("Butyrate", None),  # only C subgroups
    ("Disease", None)      # only E subgroups
]

# ---- Function to plot alpha diversity ----
def plot_alpha_diversity(df, x_col, levels=None, title=None, outfile=None):
    # Subset if levels provided
    if levels is not None:
        df = df[df[x_col].isin(levels)]
    
    # Melt for seaborn
    melted = df.melt(
        id_vars=[x_col],
        value_vars=["Shannon", "Evenness", "Chao1", "Simpson"],
        var_name="Metric",
        value_name="Diversity"
    )
    
    # Plot
    g = sns.catplot(
        data=melted,
        x=x_col, y="Diversity",
        col="Metric",
        kind="box",
        col_wrap=2,
        sharey=False,
        height=4, aspect=1.2
    )
    g.map_dataframe(sns.stripplot, x=x_col, y="Diversity", color="black", alpha=0.5)
    plt.subplots_adjust(top=0.85)
    
    if title:
        g.figure.suptitle(title)
    if outfile:
        g.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

# ---- Loop over comparisons ----
for x_col, levels in comparisons:
    name = f"{x_col}_{'_'.join(levels) if levels else 'all'}"
    title = f"Alpha Diversity - {name.replace('_', ' ')}"
    outfile = f"alpha_diversity_plots/alpha_diversity_{name}.png"
    
    plot_alpha_diversity(
        df=merged,
        x_col=x_col,
        levels=levels,
        title=title,
        outfile=outfile
    )


# Load PCoA results
pcoa_results = {
    "Bray-Curtis": Artifact.load("core-metrics-results/bray_curtis_pcoa_results.qza").view(OrdinationResults),
    "Jaccard": Artifact.load("core-metrics-results/jaccard_pcoa_results.qza").view(OrdinationResults),
    "Weighted Unifrac": Artifact.load("core-metrics-results/weighted_unifrac_pcoa_results.qza").view(OrdinationResults),
    "Unweighted Unifrac": Artifact.load("core-metrics-results/unweighted_unifrac_pcoa_results.qza").view(OrdinationResults),
}

metadata = pd.read_csv("metadata.tsv", sep="\t", index_col=0)

# Strip whitespace from all string/object columns
for col_name in metadata.select_dtypes(include=['object']).columns:
    metadata[col_name] = metadata[col_name].str.strip()
n_pcs = 19

comparisons = [
    ("Group", None),
    ("Butyrate", None),
    ("Disease", None)
]

output_dir = "beta_diversity_plots"
os.makedirs(output_dir, exist_ok=True)

for col, filter_values in comparisons:
    if col not in metadata.columns:
        print(f"Warning: Column '{col}' not found in metadata. Skipping.")
        continue

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for ax, (distance_metric, pcoa_res) in zip(axes, pcoa_results.items()):
        coords = pcoa_res.samples
        df = coords.merge(metadata, left_index=True, right_index=True)
        df.rename(columns={i: f'PC{i+1}' for i in range(n_pcs)}, inplace=True)

        # Filter if needed
        df_subset = df.copy()
        if filter_values is not None:
            df_subset = df_subset[df_subset[col].isin(filter_values)]

        if df_subset.empty:
            print(f"Warning: No data for {col} in {distance_metric}. Skipping plot.")
            ax.set_title(f"{distance_metric} (no data)")
            ax.axis('off')
            continue

        # Ensure PC1 and PC2 exist
        if "PC1" not in df_subset.columns or "PC2" not in df_subset.columns:
            print(f"Warning: PC1 or PC2 missing for {distance_metric}. Skipping.")
            ax.set_title(f"{distance_metric} (no PC1/PC2)")
            ax.axis('off')
            continue

        unique_groups = df_subset[col].dropna().unique()
        palette = sns.color_palette(n_colors=len(unique_groups))
        group_to_color = dict(zip(unique_groups, palette))

        # Scatter plot
        sns.scatterplot(
            x="PC1",
            y="PC2",
            hue=col if not df_subset[col].isnull().all() else None,
            data=df_subset,
            s=100,
            alpha=0.8,
            palette=group_to_color if not df_subset[col].isnull().all() else None,
            ax=ax
        )

        # Draw ellipses
        for group, data_subset in df_subset.groupby(col):
            if len(data_subset) < 2:
                continue
            centroid_x = data_subset["PC1"].mean()
            centroid_y = data_subset["PC2"].mean()
            width = data_subset["PC1"].std() * 2
            height = data_subset["PC2"].std() * 2
            ellipse = Ellipse(
                (centroid_x, centroid_y),
                width=width, height=height,
                edgecolor=group_to_color[group],
                facecolor='none', lw=2, alpha=0.7
            )
            ax.add_patch(ellipse)
            ax.scatter(centroid_x, centroid_y, marker='x', color='black', s=120, zorder=10)

        ax.set_title(distance_metric)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    # Combine legend outside grid
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title=col, bbox_to_anchor=(1.05, 0.5), loc='center left')

    subset_label = "all" if filter_values is None else "_".join(filter_values)
    fig.suptitle(f"PCoA comparison ({col} = {subset_label})", fontsize=16)

    plt.tight_layout(rect=[0, 0, 0.85, 0.95])

    filename = f"PCoA_comparison_{col}_{subset_label}.png"
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=300)
    plt.close()
