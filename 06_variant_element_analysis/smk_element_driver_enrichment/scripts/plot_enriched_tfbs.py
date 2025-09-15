import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse


def plot_enrichment_results(input_path, comparison_group, label_A, label_B, output_path, threshold=0.05, plot_title="TFBS Enrichment Analysis", plot_xlabel="Odds Ratio", plot_ylabel="Transcription Factors", plot_size=(10,8), plot_dpi=300, plot_format="png"):
    """
    Plots TFBS enrichment results for a specified comparison group using explicit parameters.

    Args:
        input_path (str): Path to the enrichment results file (CSV/TSV).
        comparison_group (str): The specific comparison group to plot.
        label_A (str): Label for group A.
        label_B (str): Label for group B.
        threshold (float): FDR threshold for significance.
        output_prefix (str): Prefix for output plot file.
        plot_title (str): Title for the plot.
        plot_xlabel (str): X-axis label.
        plot_ylabel (str): Y-axis label.
        plot_size (tuple): Figure size.
        plot_dpi (int): DPI for plot.
        plot_format (str): Output format (png, pdf, etc.).
    """
    df = pd.read_csv(input_path, sep=None, engine='python')
    plot_df = df[df['comparison_group'] == comparison_group].copy()
    plot_df['is_significant'] = plot_df['adjusted_pvalue'] < threshold
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=plot_size, dpi=plot_dpi)
    ax = sns.scatterplot(
        data=plot_df,
        x='odds_ratio',
        y='TFBS_name',
        size='active_overlaps',
        hue='is_significant',
        style='is_significant',
        palette={True: 'blue', False: 'grey'},
        sizes=(50, 500),
        legend=False
    )
    ax.set_title(plot_title, fontsize=16, pad=20)
    ax.set_xlabel(plot_xlabel, fontsize=12)
    ax.set_ylabel(plot_ylabel, fontsize=12)
    plt.axvline(x=1, color='red', linestyle='--', linewidth=1)
    significant_handle = plt.scatter([], [], color='blue', label=f'Significant (q < {threshold})', s=100)
    non_significant_handle = plt.scatter([], [], color='grey', label='Not Significant', s=100)
    sizes_legend = [50, 200, 500]
    overlap_handles = [plt.scatter([], [], c='black', alpha=0.5, s=s, label=str(s)) for s in sizes_legend]
    ax.legend(handles=[significant_handle, non_significant_handle], loc='lower right', title='Significance', bbox_to_anchor=(1.25, 0.5))
    ax.legend(handles=overlap_handles, loc='upper right', title='Number of Overlaps', bbox_to_anchor=(1.25, 0.75))
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot TFBS enrichment results.")
    parser.add_argument("input_path", type=str, help="Path to the enrichment results file (CSV/TSV).")
    parser.add_argument("comparison_group", type=str, help="The specific comparison group to plot.")
    parser.add_argument("label_A", type=str, help="Label for group A.")
    parser.add_argument("label_B", type=str, help="Label for group B.")
    parser.add_argument("output_path", type=str, help="Path to save the output plot.")
    parser.add_argument("--threshold", type=float, default=0.05, help="FDR threshold for significance.")
    parser.add_argument("--plot_title", type=str, default="TFBS Enrichment Analysis", help="Title for the plot.")
    parser.add_argument("--plot_xlabel", type=str, default="Odds Ratio", help="X-axis label.")
    parser.add_argument("--plot_ylabel", type=str, default="Transcription Factors", help="Y-axis label.")
    parser.add_argument("--plot_size", type=float, nargs=2, default=[10,8], help="Figure size (width height).")
    parser.add_argument("--plot_dpi", type=int, default=300, help="DPI for plot.")
    parser.add_argument("--plot_format", type=str, default="png", help="Output format (png, pdf, etc.).")
    args = parser.parse_args()
    plot_enrichment_results(
        args.input_path,
        args.comparison_group,
        args.label_A,
        args.label_B,
        args.output_path,
        threshold=args.threshold,
        plot_title=args.plot_title,
        plot_xlabel=args.plot_xlabel,
        plot_ylabel=args.plot_ylabel,
        plot_size=tuple(args.plot_size),
        plot_dpi=args.plot_dpi,
        plot_format=args.plot_format
    )

# example:
# | TFBS_name | odds_ratio | adjusted_pvalue | active_overlaps | comparison_group |
# |---|---|---|---|---|
# | TF_A | 5.2 | 0.0001 | 150 | active_vs_inactive |
# | TF_B | 3.1 | 0.002 | 85 | active_vs_inactive |