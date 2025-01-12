import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os

# TCGA settings
tcga_cancer = "Pancreatic Ductal Adenocarcinoma"
tcga_clinical = os.path.join("data_clinical", "nationwidechildrens.org_clinical_patient_paad.txt")
tcga_snv_dir = 'data_snv_PAAD'
tcga_include = {
    'histologic_diagnosis': ["Pancreas-Adenocarcinoma Ductal Type"],
}

# Plot settings remain the same as before
plt.style.use('ggplot')
plot_format = 'png'
plot_size = (7.0, 3.5)
plot_dpi = 600
plot_fontsize = 6
gene_fontsize = 6

# Mutation type mapping remains the same
mutation_categories = {
    'Missense_Mutation': ('Missense', '#336699'),
    'In_Frame_Del': ('Inframe', '#009999'),
    'In_Frame_Ins': ('Inframe', '#009999'),
    'Splice_Site': ('Critical Site', '#cc9933'),
    'Translation_Start_Site': ('Critical Site', '#cc9933'),
    'Nonstop_Mutation': ('Critical Site', '#cc9933'),
    'Frame_Shift_Del': ('Frameshift', '#ff6600'),
    'Frame_Shift_Ins': ('Frameshift', '#ff6600'),
    'Nonsense_Mutation': ('Nonsense', '#cc0033'),
    'Silent': ('Synonymous', '#d2dae2'),
}

def read_clinical_data(clinical_file, include_histology_dict):
    """Read clinical data and filter based on multiple histology columns.
    
    Args:
        clinical_file (str): Path to clinical data file
        include_histology_dict (dict): Dictionary mapping histology column names to lists of included values
    
    Returns:
        list: List of case IDs meeting the histology criteria
    """
    # Read the file with the first row as header, skip the 2nd and 3rd rows
    df = pd.read_csv(clinical_file, sep='\t', skiprows=[1,2])
    
    # Create mask for each histology column
    masks = []
    for col, values in include_histology_dict.items():
        masks.append(df[col].isin(values))
    
    # Combine masks with OR operation
    final_mask = masks[0]
    for mask in masks[1:]:
        final_mask = final_mask | mask
    
    # Filter based on combined mask
    filtered_df = df[final_mask]
    
    # Get case IDs (first column)
    return filtered_df.iloc[:, 1].tolist()

# Rest of the code remains unchanged
def read_maf_files(directory, case_ids):
    """Read MAF files for specific cases."""
    all_mutations = []
    
    for case_id in case_ids:
        maf_file = os.path.join(directory, f"{case_id}.maf")
        if os.path.exists(maf_file):
            df = pd.read_csv(maf_file, sep='\t', comment='#', low_memory=False)
            required_cols = ['Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode']
            all_mutations.append(df[required_cols])
            if len(df) == 0:
                print(f"{maf_file} is empty!")
        else:
            print(f"{maf_file} not found!")
    
    if not all_mutations:
        raise ValueError("No matching MAF files found for the filtered cases")
        
    return pd.concat(all_mutations, ignore_index=True)

def create_mutation_matrix(mutations_df, top_n_genes=50):
    """Create a mutation matrix for visualization."""
    # Filter out mutations not in our categories
    mutations_df = mutations_df[mutations_df['Variant_Classification'].isin(mutation_categories.keys())]
    
    # Get total number of cases
    total_cases = len(mutations_df['Tumor_Sample_Barcode'].unique())
    
    # Calculate mutation frequency per gene
    gene_freq = mutations_df.groupby('Hugo_Symbol')['Tumor_Sample_Barcode'].nunique()
    gene_freq_pct = (gene_freq / total_cases * 100).round(1)
    top_genes = gene_freq_pct.nlargest(top_n_genes).index
    
    # Create gene x sample matrix
    mutation_matrix = pd.crosstab(
        mutations_df['Hugo_Symbol'],
        mutations_df['Tumor_Sample_Barcode']
    ).loc[top_genes]
    
    # Store mutation types
    mutation_types = defaultdict(dict)
    priority_order = ['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 
                     'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation',
                     'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Silent']
    
    for gene in top_genes:
        gene_mutations = mutations_df[mutations_df['Hugo_Symbol'] == gene]
        for sample in gene_mutations['Tumor_Sample_Barcode'].unique():
            sample_mutations = gene_mutations[gene_mutations['Tumor_Sample_Barcode'] == sample]
            for mut_type in priority_order:
                if mut_type in sample_mutations['Variant_Classification'].values:
                    mutation_types[(gene, sample)] = mut_type
                    break
    
    return mutation_matrix, mutation_types, gene_freq_pct

def sort_samples(mutation_matrix):
    """Sort samples based on their mutation patterns."""
    binary_matrix = (mutation_matrix > 0).astype(int)
    weights = np.zeros(binary_matrix.shape[1])
    for i, gene_row in enumerate(binary_matrix.values):
        weights += gene_row * 2**(binary_matrix.shape[0] - i - 1)
    return np.argsort(-weights)

def create_mutation_landscape(mutations_df, output_file, top_n_genes=30):
    """Generate mutation landscape plot."""
    mutation_matrix, mutation_types, gene_freq_pct = create_mutation_matrix(mutations_df, top_n_genes)
    sorted_cols = sort_samples(mutation_matrix)
    ordered_matrix = mutation_matrix.iloc[:, sorted_cols]
    
    plt.figure(figsize=plot_size)
    ax = plt.gca()
    
    # Plot mutations
    for i, gene in enumerate(ordered_matrix.index[::-1]):
        for j, sample in enumerate(ordered_matrix.columns):
            if ordered_matrix.loc[gene, sample] > 0:
                mut_type = mutation_types.get((gene, sample))
                if mut_type in mutation_categories:
                    category, color = mutation_categories[mut_type]
                    ax.add_patch(plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor=color, edgecolor='none'
                    ))
    
    # Add horizontal lines between genes
    for i in range(len(ordered_matrix.index) - 1):
        ax.axhline(y=i+0.5, color='white', linewidth=1, alpha=0.5, zorder=1)
    
    # Customize plot
    ax.set_xlim(-0.5, len(ordered_matrix.columns) - 0.5)
    ax.set_ylim(-0.5, len(ordered_matrix.index) - 0.5)
    
    # Add gene labels with frequency
    yticks_pos = range(len(ordered_matrix.index))
    gene_labels = [f"{gene} ({gene_freq_pct[gene]}%)" for gene in ordered_matrix.index[::-1]]
    ax.set_yticks(yticks_pos)
    ax.set_yticklabels(gene_labels, fontsize=gene_fontsize)
    
    # Format axes
    ax.yaxis.set_ticks_position('none')
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_xlabel(f'Samples (n={len(ordered_matrix.columns)})', fontsize=plot_fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)
    
    # Add legend
    legend_categories = [
        ('Missense', '#336699'),
        ('Inframe', '#009999'),
        ('Critical Site', '#cc9933'),
        ('Frameshift', '#ff6600'),
        ('Nonsense', '#cc0033'),
        ('Synonymous', '#d2dae2'),
    ]
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                      for label, color in legend_categories]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1),
             loc='upper left', fontsize=plot_fontsize)
    
    plt.title(tcga_cancer, fontsize=plot_fontsize+2)
    plt.tight_layout()
    plt.savefig(output_file, format=plot_format, dpi=plot_dpi, bbox_inches='tight')
    plt.close()

def main():
    # Read and filter clinical data
    print("Reading clinical data...")
    case_ids = read_clinical_data(tcga_clinical, tcga_include)
    print(f"Found {len(case_ids)} cases matching histology criteria")
    
    # Read filtered MAF files
    print("Reading MAF files for filtered cases...")
    mutations_df = read_maf_files(tcga_snv_dir, case_ids)
    
    # Create the plot
    print("Generating mutation landscape plot...")
    output_file = tcga_cancer.replace(' ', '_') + '.' + plot_format
    create_mutation_landscape(mutations_df, output_file, top_n_genes=30)
    print(f"Plot saved as: {output_file}")

if __name__ == "__main__":
    main()