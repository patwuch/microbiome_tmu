import os
import pandas as pd
import biom
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

output_dir = "taxa_barplots"
os.makedirs(output_dir, exist_ok=True)

# ===========================
# USER SETTINGS
# ===========================
feature_table_biom = "exported/exported-feature-table/feature-table.biom"
taxonomy_tsv = "exported/exported-taxonomy/taxonomy.tsv"
top_n = 20            # Number of taxa to show in legend

# ===========================
# LOAD DATA
# ===========================
table = biom.load_table(feature_table_biom)
df = pd.DataFrame(table.matrix_data.toarray().T, 
                  index=table.ids(axis='sample'), 
                  columns=table.ids(axis='observation'))

taxonomy = pd.read_csv(taxonomy_tsv, sep='\t', index_col=0)

# 👇 FIX FOR AttributeError: 'float' object has no attribute 'lower'
# Ensure the Taxon column is treated as a string and missing values are empty strings.
taxonomy['Taxon'] = taxonomy['Taxon'].astype(str).fillna('').str.strip()

level_dict = {
    "Kingdom": 0,
    "Phylum": 1,
    "Class": 2,
    "Order": 3,
    "Family": 4,
    "Genus": 5,
}

# Define the order of levels for easy lookup
LEVEL_NAMES = list(level_dict.keys())
LEVEL_INDICES = list(level_dict.values())

# Helper function to find the highest-level assignment for a feature
def relabel_unassigned_taxa(row, current_level_index):
    """
    Looks up the next available, higher taxonomic level for unassigned features,
    safely handling taxonomy strings that are shorter than 6 levels.
    """
    
    # Safely split the taxonomy string and pad with 'Unassigned' if needed.
    taxa_list = [t.strip() for t in row['Taxon'].split(';')]
    # Ensure the list is padded up to the maximum number of levels (6)
    while len(taxa_list) < len(LEVEL_NAMES): 
        taxa_list.append("Unassigned")
        
    # 1. Get the current level assignment
    current_label = taxa_list[current_level_index]
    
    # 2. Check if the current label is unassigned (or similar)
    if 'unassigned' in current_label.lower() or current_label.strip() == '':
        # 3. Iterate through higher levels (from current_level_index - 1 down to 0)
        for i in range(current_level_index - 1, -1, -1):
            higher_level_label = taxa_list[i]
            
            # 4. If a higher level is assigned, use it as the prefix
            if 'unassigned' not in higher_level_label.lower() and higher_level_label.strip() != '':
                higher_level_name = LEVEL_NAMES[i]
                # Format: "Unassigned ({HigherLevel} {HigherLabel})"
                return f"Unassigned ({higher_level_name} {higher_level_label})"
        
        # 5. If no higher level is assigned (all are unassigned), return a generic label
        return "Unassigned (All)"
    
    # 6. If the current label is not unassigned, return it as is
    return current_label


# ===========================
# LOOP THROUGH TAXONOMIC LEVELS
# ===========================
for tax_level, level_index in level_dict.items():
    print(f"Processing {tax_level}...")

    # --- LOGIC FOR RELABELING ---
    if tax_level != "Kingdom": # Kingdom level is the highest, no higher level to check
        # Apply the relabeling function to create the new, more descriptive tax column
        taxonomy[tax_level] = taxonomy.apply(
            relabel_unassigned_taxa, 
            axis=1, 
            current_level_index=level_index
        )
    else:
        # Kingdom level is simpler: just get the label and handle Unassigned
        # We rely on the padding in the function below to ensure safe access,
        # but for Kingdom we can still use the simpler logic for readability.
        taxonomy[tax_level] = taxonomy['Taxon'].str.split(';').str[level_index].str.strip().fillna("Unassigned")
        taxonomy[tax_level] = taxonomy[tax_level].apply(
            lambda x: "Unassigned (All)" if 'unassigned' in x.lower() or x.strip() == '' else x
        )
    
    # --- Special case for Genus (Family|Genus) still needed ---
    if tax_level == "Genus":
        # Extract the Family name. We use the safe list access within a lambda.
        taxonomy["Family_prefix"] = taxonomy['Taxon'].apply(
            lambda x: ([t.strip() for t in x.split(';')] + ['Unassigned'] * 6)[level_dict["Family"]]
        )
        
        # We need the *relabelled* genus for the second part of the string
        # Combine Family prefix with the (potentially relabelled) Genus name
        taxonomy[tax_level] = taxonomy.apply(
            lambda row: f"{row['Family_prefix']}|{row[tax_level]}" 
                        if 'unassigned' not in row['Family_prefix'].lower() and row['Family_prefix'].strip() != ''
                        else row[tax_level], # If Family is unassigned, just use the Genus label
            axis=1
        )
    # -----------------------------------------------------------

    # Aggregate counts by taxonomic level
    df_tax = df.groupby(taxonomy[tax_level], axis=1).sum()
    df_tax_norm = df_tax.div(df_tax.sum(axis=1), axis=0)

    # Identify top N taxa by mean relative abundance
    mean_abundance = df_tax_norm.mean(axis=0)
    top_taxa = mean_abundance.sort_values(ascending=False).head(top_n).index

    # Prepare dataframe for plotting
    df_top = df_tax_norm[top_taxa].copy()
    df_top['Other'] = df_tax_norm.drop(columns=top_taxa, errors='ignore').sum(axis=1)

    # Color map
    n_colors = df_top.shape[1]
    cmap = cm.get_cmap('tab20', n_colors)
    colors = [cmap(i) for i in range(n_colors)]

    # Plot manually stacked bars
    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(df_top.shape[0])
    for i, col in enumerate(df_top.columns):
        ax.bar(
            df_top.index,
            df_top[col],
            bottom=bottom,
            color=colors[i],
            label=col,
            width=0.8,
            edgecolor="none",
            linewidth=0,
            antialiased=False
        )
        bottom += df_top[col].values

    # Axis labels and title
    ax.set_ylabel("Relative abundance")
    ax.set_xlabel("Samples")
    plt.xticks(rotation=90)
    plt.title(f"Relative abundance at {tax_level} level")

    # Legend in descending abundance order
    handles, labels = ax.get_legend_handles_labels()
    handles, labels = handles[::-1], labels[::-1]
    ax.legend(
        handles, labels,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize=8,
        title="Taxa",
        title_fontsize=9,
        frameon=False
    )

    plt.tight_layout()
    output_png = f"{output_dir}/taxa_barplot_{tax_level}.png"
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()