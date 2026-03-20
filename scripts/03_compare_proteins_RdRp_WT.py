#!/usr/bin/env python3
"""
Multi-WT Protein Disorder Comparison and Ranking
Compares disorder predictions across all WT protein variants and ranks them

This script:
1. Aggregates disorder data from all WT proteins
2. Ranks WT proteins by intrinsic disorder level
3. Generates comparison figures and statistics
4. Creates summary ranking table

Usage:
    python3 03_compare_proteins_WT.py
    python3 03_compare_proteins_WT.py --output_dir comparison_results_WT
"""
import argparse
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle


class ProteinDisorderComparison:
    """Compare disorder predictions across multiple WT proteins."""
    
    def __init__(self, output_dir: str = "protein_comparison_WT"):
        """Initialize comparison analysis with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        self.plots_dir = self.output_dir / "plots"
        self.data_dir = self.output_dir / "data"
        self.plots_dir.mkdir(exist_ok=True)
        self.data_dir.mkdir(exist_ok=True)
        
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.size'] = 10
        
        print(f"Output directory: {self.output_dir}")
        print(f"  - Plots: {self.plots_dir}")
        print(f"  - Data: {self.data_dir}")
    
    def find_protein_data(self, disorder_analysis_dir: str = "disorder_analysis_WT") -> dict:
        """
        Find all WT protein analysis results.
        
        Parameters
        ----------
        disorder_analysis_dir : str
            Path to the disorder_analysis_WT directory containing all WT protein results
        
        Returns
        -------
        proteins : Dict[str, Dict]
            Dictionary with WT protein names and their analysis data
        """
        analysis_path = Path(disorder_analysis_dir)
        
        if not analysis_path.exists():
            print(f"✗ Analysis directory not found: {analysis_path}")
        print(f"  Expected: protein_analysis/")
            return {}
        
        proteins = {}
        
        # Find all protein subdirectories
        protein_dirs = sorted([d for d in analysis_path.iterdir() if d.is_dir()])
        
        for protein_dir in protein_dirs:
            # Try both possible filenames
            data_file = protein_dir / 'data' / 'disorder_summary_statistics.csv'
            if not data_file.exists():
                data_file = protein_dir / 'data' / 'disorder_statistics.csv'
            
            if data_file.exists():
                try:
                    df = pd.read_csv(data_file)
                    proteins[protein_dir.name] = {'data': df}
                    print(f"  ✓ Loaded data for {protein_dir.name}")
                except Exception as e:
                    print(f"  ✗ Error loading {protein_dir.name}: {e}")
        
        return proteins
    
    def extract_disorder_metrics(self, proteins: dict) -> pd.DataFrame:
        """
        Extract disorder metrics from all WT proteins (AlphaFold2 focus).
        
        Parameters
        ----------
        proteins : Dict[str, Dict]
            Dictionary of WT protein analysis data
        
        Returns
        -------
        df : pd.DataFrame
            DataFrame with disorder metrics for all WT proteins
        """
        metrics_list = []
        
        for protein_name, data in sorted(proteins.items()):
            df = data['data']
            
            # Extract metrics from the statistics file
            # The CSV has rows for each tool, we want the AlphaFold2 row
            for idx, row in df.iterrows():
                if 'AlphaFold2' in str(row.get('Tool', '')):
                    # Extract values - handle the disorder percentage
                    disorder_str = str(row.get('Disorder Percentage', '0%')).strip().rstrip('%')
                    try:
                        disorder_pct = float(disorder_str)
                    except:
                        disorder_pct = 0
                    
                    # Extract residue counts from the format "744 / 1727"
                    residue_str = str(row.get('Disordered Residues (>0.5)', '0 / 0'))
                    parts = residue_str.split('/')
                    disordered_count = int(parts[0].strip()) if len(parts) > 0 else 0
                    total_residues = int(parts[1].strip()) if len(parts) > 1 else 1
                    
                    mean_disorder = float(row.get('Mean Disorder', 0))
                    std_disorder = float(row.get('Std Dev', 0))
                    
                    metrics_list.append({
                        'Protein': protein_name,
                        'Disorder %': disorder_pct,
                        'Mean Disorder': mean_disorder,
                        'Std Dev': std_disorder,
                        'Disordered Residues': disordered_count,
                        'Total Residues': total_residues
                    })
                    break
        
        return pd.DataFrame(metrics_list)
    
    def extract_multi_tool_metrics(self, proteins: dict) -> pd.DataFrame:
        """
        Extract disorder metrics from all tools for all WT proteins.
        
        Parameters
        ----------
        proteins : Dict[str, Dict]
            Dictionary of WT protein analysis data
        
        Returns
        -------
        df : pd.DataFrame
            DataFrame with multi-tool metrics for all WT proteins
        """
        metrics_list = []
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        for protein_name, data in sorted(proteins.items()):
            df = data['data']
            
            # Extract metrics for each tool
            tool_means = {}
            for idx, row in df.iterrows():
                tool_name = str(row.get('Tool', ''))
                mean_disorder = float(row.get('Mean Disorder', 0))
                
                if 'AlphaFold2' in tool_name:
                    tool_means['AlphaFold2 (pLDDT)'] = mean_disorder
                elif 'IUPred3' in tool_name:
                    tool_means['IUPred3'] = mean_disorder
                elif 'DISOPRED3' in tool_name:
                    tool_means['DISOPRED3'] = mean_disorder
                elif 'SPOT' in tool_name:
                    tool_means['SPOT-Disorder'] = mean_disorder
            
            # Calculate consensus (mean of all tools)
            if tool_means:
                consensus = np.mean(list(tool_means.values()))
                tool_std = np.std(list(tool_means.values()))
                tool_agreement = 1.0 - tool_std  # Higher = better agreement
                
                metrics_list.append({
                    'Protein': protein_name,
                    'AlphaFold2 (pLDDT)': tool_means.get('AlphaFold2 (pLDDT)', 0),
                    'IUPred3': tool_means.get('IUPred3', 0),
                    'DISOPRED3': tool_means.get('DISOPRED3', 0),
                    'SPOT-Disorder': tool_means.get('SPOT-Disorder', 0),
                    'Consensus': consensus,
                    'Tool_Agreement': tool_agreement
                })
        
        return pd.DataFrame(metrics_list)
    
    def rank_proteins(self, metrics_df: pd.DataFrame) -> pd.DataFrame:
        """
        Rank WT proteins by disorder level.
        
        Parameters
        ----------
        metrics_df : pd.DataFrame
            DataFrame with disorder metrics
        
        Returns
        -------
        ranked_df : pd.DataFrame
            Ranked DataFrame with disorder rankings
        """
        # Sort by disorder percentage (descending)
        ranked_df = metrics_df.sort_values('Disorder %', ascending=False).reset_index(drop=True)
        ranked_df['Rank'] = range(1, len(ranked_df) + 1)
        
        # Reorder columns
        ranked_df = ranked_df[['Rank', 'Protein', 'Disorder %', 'Mean Disorder', 
                               'Std Dev', 'Disordered Residues', 'Total Residues']]
        
        return ranked_df
    
    def print_ranking_table(self, ranked_df: pd.DataFrame) -> None:
        """Print ranking table to console."""
        print("\n" + "="*100)
        print("WT PROTEIN DISORDER RANKING")
        print("="*100)
        print("\nRanked from MOST to LEAST Intrinsic Disorder:\n")
        print(ranked_df.to_string(index=False))
        
        # Add summary statistics
        print("\n" + "="*100)
        print("SUMMARY STATISTICS")
        print("="*100)
        if len(ranked_df) > 0:
            print(f"\nHighest disorder:  {ranked_df.iloc[0]['Protein']} ({ranked_df.iloc[0]['Disorder %']:.1f}%)")
            print(f"Lowest disorder:   {ranked_df.iloc[-1]['Protein']} ({ranked_df.iloc[-1]['Disorder %']:.1f}%)")
            print(f"Mean disorder across WT proteins: {ranked_df['Disorder %'].mean():.1f}%")
            if len(ranked_df) > 1:
                print(f"Std dev across WT proteins: {ranked_df['Disorder %'].std():.1f}%")
    
    def plot_disorder_ranking_bars(self, ranked_df: pd.DataFrame) -> None:
        """Create bar plot ranking WT proteins by disorder."""
        if len(ranked_df) == 0:
            print("  ⚠ No data to plot")
            return
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Create smoother color gradient from red (high disorder) to green (low disorder)
        norm = plt.Normalize(vmin=ranked_df['Disorder %'].min(), vmax=ranked_df['Disorder %'].max())
        colors = plt.cm.RdYlGn_r(norm(ranked_df['Disorder %']))
        
        # Create bars
        bars = ax.barh(range(len(ranked_df)), ranked_df['Disorder %'], color=colors, edgecolor='black', linewidth=1.5)
        
        # Create clean labels with rank and protein name
        clean_labels = [f"{int(row['Rank'])}. {row['Protein']}" for _, row in ranked_df.iterrows()]
        
        # Customize axes
        ax.set_yticks(range(len(ranked_df)))
        ax.set_yticklabels(clean_labels, fontsize=11, fontweight='bold')
        ax.set_xlabel('Intrinsic Disorder (%)', fontsize=12, fontweight='bold')
        ax.set_title('WT Protein Fragment Ranking by Intrinsic Disorder\n(AlphaFold2 pLDDT-based prediction)',
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add value labels on bars
        for i, (idx, row) in enumerate(ranked_df.iterrows()):
            ax.text(row['Disorder %'] + 0.5, i, f"{row['Disorder %']:.1f}%", 
                   va='center', fontsize=10, fontweight='bold')
        
        # Add grid
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, ranked_df['Disorder %'].max() * 1.15)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'wt_protein_disorder_ranking.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'wt_protein_disorder_ranking.png'}")
        plt.close()
    
    def plot_disorder_vs_residues(self, ranked_df: pd.DataFrame) -> None:
        """Create scatter plot: disorder % vs number of residues for WT proteins."""
        if len(ranked_df) == 0:
            print("  ⚠ No data to plot")
            return
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot with smooth color gradient based on disorder %
        norm = plt.Normalize(vmin=ranked_df['Disorder %'].min(), vmax=ranked_df['Disorder %'].max())
        scatter = ax.scatter(ranked_df['Total Residues'], ranked_df['Disorder %'],
                           s=500, c=ranked_df['Disorder %'], cmap='RdYlGn_r', norm=norm,
                           edgecolors='black', linewidth=2, alpha=0.8)
        
        # Add protein labels
        for idx, row in ranked_df.iterrows():
            ax.annotate(row['Protein'], 
                       (row['Total Residues'], row['Disorder %']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=9, fontweight='bold')
        
        # Customize
        ax.set_xlabel('Protein Length (Residues)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Intrinsic Disorder (%)', fontsize=12, fontweight='bold')
        ax.set_title('WT Protein Disorder vs Fragment Size',
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add colorbar with proper label
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Disorder (%)', fontsize=11, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'wt_disorder_vs_residues.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'wt_disorder_vs_residues.png'}")
        plt.close()
    
    def plot_multi_tool_comparison(self, multi_tool_df: pd.DataFrame) -> None:
        """Create comparison plot showing all 4 tools for each WT protein."""
        if len(multi_tool_df) == 0:
            print("  ⚠ No data to plot")
            return
        
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        # Create bar plot comparing mean disorder across tools
        fig, ax = plt.subplots(figsize=(12, 6))
        
        x = np.arange(len(multi_tool_df))
        width = 0.2
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        
        for i, tool in enumerate(tools):
            values = multi_tool_df[tool].values if tool in multi_tool_df.columns else [0]*len(multi_tool_df)
            ax.bar(x + i*width, values, width, label=tool, color=colors[i], edgecolor='black', linewidth=1)
        
        ax.set_xlabel('WT Protein Fragment', fontsize=12, fontweight='bold')
        ax.set_ylabel('Mean Disorder Score', fontsize=12, fontweight='bold')
        ax.set_title('Multi-Tool Disorder Prediction Comparison\nAcross All WT Protein Variants',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x + width * 1.5)
        ax.set_xticklabels(multi_tool_df['Protein'], rotation=45, ha='right', fontweight='bold')
        ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'wt_multi_tool_comparison_bars.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'wt_multi_tool_comparison_bars.png'}")
        plt.close()
    
    def save_ranking_table(self, ranked_df: pd.DataFrame) -> None:
        """Save WT ranking table as CSV."""
        csv_path = self.data_dir / 'wt_protein_disorder_ranking.csv'
        ranked_df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")
    
    def save_multi_tool_comparison(self, multi_tool_df: pd.DataFrame) -> None:
        """Save multi-tool comparison data as CSV."""
        csv_path = self.data_dir / 'wt_protein_multi_tool_comparison.csv'
        multi_tool_df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")
    
    def print_multi_tool_summary(self, multi_tool_df: pd.DataFrame) -> None:
        """Print multi-tool analysis summary to console."""
        if len(multi_tool_df) == 0:
            return
        
        print("\n" + "="*100)
        print("MULTI-TOOL WT DISORDER PREDICTION COMPARISON")
        print("="*100)
        
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        print("\n📊 Mean Disorder Scores by Tool (Across All WT Proteins):\n")
        
        tool_stats = {}
        for tool in tools:
            if tool in multi_tool_df.columns:
                values = multi_tool_df[tool].values
                tool_stats[tool] = {
                    'mean': np.mean(values),
                    'std': np.std(values),
                    'min': np.min(values),
                    'max': np.max(values)
                }
        
        for tool, stats in tool_stats.items():
            print(f"{tool}:")
            print(f"  Mean: {stats['mean']:.4f}")
            print(f"  Range: [{stats['min']:.4f}, {stats['max']:.4f}]")
            print(f"  Std Dev: {stats['std']:.4f}\n")
        
        if 'Tool_Agreement' in multi_tool_df.columns:
            print(f"Inter-Tool Agreement (Consensus Quality):")
            print(f"  Mean Agreement Score: {multi_tool_df['Tool_Agreement'].mean():.4f}")
            print(f"  Range: [{multi_tool_df['Tool_Agreement'].min():.4f}, {multi_tool_df['Tool_Agreement'].max():.4f}]")
    
    def run_comparison(self, disorder_analysis_dir: str = "disorder_analysis_WT") -> None:
        """Run full comparison analysis on all WT proteins."""
        print(f"\n{'='*80}")
        print("Locating WT Protein Analysis Results")
        print(f"{'='*80}\n")
        
        proteins = self.find_protein_data(disorder_analysis_dir)
        
        if not proteins:
            print("\n✗ No WT protein analysis results found")
            return
        
        print(f"\n✓ Found {len(proteins)} WT proteins")
        
        # Extract metrics and rank
        print("\n" + "="*80)
        print("Extracting Disorder Metrics")
        print("="*80)
        
        metrics_df = self.extract_disorder_metrics(proteins)
        multi_tool_df = self.extract_multi_tool_metrics(proteins)
        
        if len(metrics_df) == 0:
            print("✗ No metrics extracted")
            return
        
        ranked_df = self.rank_proteins(metrics_df)
        
        # Print rankings
        self.print_ranking_table(ranked_df)
        self.print_multi_tool_summary(multi_tool_df)
        
        # Create plots
        print("\n" + "="*80)
        print("Generating Comparison Figures")
        print("="*80)
        
        self.plot_disorder_ranking_bars(ranked_df)
        self.plot_disorder_vs_residues(ranked_df)
        self.plot_multi_tool_comparison(multi_tool_df)
        
        # Save data
        print("\n" + "="*80)
        print("Saving Data Files")
        print("="*80)
        
        self.save_ranking_table(ranked_df)
        self.save_multi_tool_comparison(multi_tool_df)
        
        print(f"\n✓ WT Protein comparison analysis complete!")
        print(f"  Results saved in: {self.output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Multi-WT Protein Disorder Comparison and Ranking',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: analyze disorder_analysis_WT/ and save to protein_comparison_WT/
  python3 03_compare_proteins_WT.py

  # Custom output directory
  python3 03_compare_proteins_WT.py --output_dir my_wt_comparison

  # Custom analysis directory
  python3 03_compare_proteins_WT.py --analysis_dir /path/to/disorder_analysis_WT
        """
    )
    
    parser.add_argument('--analysis_dir', default='disorder_analysis_WT',
                       help='Path to disorder_analysis_WT directory (default: disorder_analysis_WT)')
    parser.add_argument('--output_dir', default='protein_comparison_WT',
                       help='Output directory for WT comparison results (default: protein_comparison_WT)')
    
    args = parser.parse_args()
    
    # Run comparison
    comparison = ProteinDisorderComparison(output_dir=args.output_dir)
    comparison.run_comparison(disorder_analysis_dir=args.analysis_dir)


if __name__ == '__main__':
    main()
