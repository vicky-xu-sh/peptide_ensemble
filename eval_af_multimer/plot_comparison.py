import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse


# Function to create box plot
def create_box_plot(metric, data, output_path):
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='source', y=metric, data=data)
    plt.title(f'Box Plot of {metric}')
    plt.xlabel('Source')
    plt.ylabel(metric)
    plt.savefig(output_path)
    plt.close()

# Function to create scatter plot
def create_scatter_plot(metric, data, output_path):
    plt.figure(figsize=(8, 6))
    plt.scatter(data[f'{metric}_baseline'], data[f'{metric}_ensemble'], alpha=0.5)
    plt.plot([data[f'{metric}_baseline'].min(), data[f'{metric}_baseline'].max()],
             [data[f'{metric}_baseline'].min(), data[f'{metric}_baseline'].max()],
             'r--', label='y=x')
    plt.title(f'Scatter Plot of {metric}: Baseline vs Ensemble')
    plt.xlabel(f'{metric} (Baseline)')
    plt.ylabel(f'{metric} (Ensemble)')
    plt.legend()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot metrics comparison plots.")
    parser.add_argument("--csv_dir", 
        type=str,
        default="./",
        help="Directory containing the two summary metrics CSV files.",
    )
    parser.add_argument("--baseline_csv_name",
        type=str,
        default="af-multimer_results_summary_baseline.csv",
        help="Name of the baseline CSV file",
    )
    parser.add_argument("--ensemble_csv_name",
        type=str,
        default="af-multimer_results_summary_ensemble.csv",
        help="Name of the ensemble CSV file",
    )
    parser.add_argument("--output_dir", 
        type=str,
        default="../plots",
        help="Directory to save the output plots.",
    )

    args = parser.parse_args()

    # Set up the workspace
    baseline_file = os.path.join(args.csv_dir, args.baseline_csv_name)
    ensemble_file = os.path.join(args.csv_dir, args.ensemble_csv_name)
    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    baseline_df = pd.read_csv(baseline_file)
    ensemble_df = pd.read_csv(ensemble_file)

    # Rename columns for clarity
    baseline_df['source'] = 'baseline'
    ensemble_df['source'] = 'ensemble'

    # Merge on 'id' to get common IDs
    merged_df = pd.merge(baseline_df, ensemble_df, on='id', suffixes=('_baseline', '_ensemble'))

    # Metrics to plot
    metrics = ['plddt', 'ptm', 'iptm', 'peptide_plddt']

    # Prepare data for box plots: melt the merged data
    box_data = []
    for _, row in merged_df.iterrows():
        for source in ['baseline', 'ensemble']:
            for metric in metrics:
                value = row[f'{metric}_{source}']
                if pd.notna(value):
                    box_data.append({'source': source, 'metric': metric, 'value': value})

    box_df = pd.DataFrame(box_data)

    # Generate plots
    for metric in metrics:
        # Box plot
        metric_box_data = box_df[box_df['metric'] == metric]
        box_output = os.path.join(args.output_dir, f'{metric}_box_plot.png')
        create_box_plot('value', metric_box_data, box_output)
        
        # Scatter plot
        scatter_output = os.path.join(args.output_dir, f'{metric}_scatter_plot.png')
        create_scatter_plot(metric, merged_df, scatter_output)

    print("Plots generated and saved in", args.output_dir)