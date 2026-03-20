#!/usr/bin/env python3
"""
Multi-Protein Disorder Comparison and Ranking
Compares disorder predictions across all 12 protein fragments and ranks them

This script:
1. Aggregates disorder data from all 12 proteins
2. Ranks proteins by intrinsic disorder level
3. Generates comparison figures and statistics
4. Creates summary ranking table

Usage:
    python3 03_compare_proteins.py
    python3 03_compare_proteins.py --output_dir comparison_results
"""

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Tuple
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import seaborn as sns


class ProteinDisorderComparison:
    """Compare disorder predictions across multiple proteins."""
    
    def __init__(self, output_dir: str = "protein_comparison"):
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
    
    def find_protein_data(self, disorder_analysis_dir: str = "protein_analysis") -> Dict[str, Dict]:
        """
        Find all protein analysis results.
        
        Parameters
        ----------
        disorder_analysis_dir : str
            Path to the disorder_analysis directory containing all protein results
        
        Returns
        -------
        proteins : Dict[str, Dict]
            Dictionary with protein names and their analysis data
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
            protein_name = protein_dir.name
            csv_file = protein_dir / "data" / "disorder_summary_statistics.csv"
            
            if csv_file.exists():
                try:
                    df = pd.read_csv(csv_file)
                    proteins[protein_name] = {
                        'dir': protein_dir,
                        'stats': df,
                        'csv_path': csv_file
                    }
                except Exception as e:
                    print(f"⚠ Error reading {protein_name}: {e}")
        
        return proteins
    
    def extract_disorder_metrics(self, proteins: Dict[str, Dict]) -> pd.DataFrame:
        """
        Extract disorder metrics from all proteins and all tools.
        
        Parameters
        ----------
        proteins : Dict[str, Dict]
            Dictionary of protein analysis data
        
        Returns
        -------
        df : pd.DataFrame
            DataFrame with disorder metrics for all proteins (AlphaFold2 focus)
        """
        metrics_list = []
        
        for protein_name, data in sorted(proteins.items()):
            stats_df = data['stats']
            
            # Extract AlphaFold2 (pLDDT) row
            af_row = stats_df[stats_df['Tool'] == 'AlphaFold2 (pLDDT)']
            
            if not af_row.empty:
                af_row = af_row.iloc[0]
                
                # Parse disorder percentage (remove % and convert to float)
                disorder_pct = af_row['Disorder Percentage']
                if isinstance(disorder_pct, str):
                    disorder_pct = float(disorder_pct.rstrip('%'))
                
                # Parse disordered residues count
                disorder_residues = af_row['Disordered Residues (>0.5)']
                if isinstance(disorder_residues, str):
                    total_residues = int(disorder_residues.split('/')[1])
                    disordered_count = int(disorder_residues.split('/')[0])
                else:
                    disordered_count = 0
                    total_residues = 0
                
                # Parse mean and std dev
                mean_disorder = float(af_row['Mean Disorder'])
                std_disorder = float(af_row['Std Dev'])
                
                metrics_list.append({
                    'Protein': protein_name,
                    'Mean Disorder': mean_disorder,
                    'Std Dev': std_disorder,
                    'Disorder %': disorder_pct,
                    'Disordered Residues': disordered_count,
                    'Total Residues': total_residues,
                })
        
        return pd.DataFrame(metrics_list)
    
    def extract_multi_tool_metrics(self, proteins: Dict[str, Dict]) -> pd.DataFrame:
        """
        Extract disorder metrics from all tools for all proteins.
        
        Parameters
        ----------
        proteins : Dict[str, Dict]
            Dictionary of protein analysis data
        
        Returns
        -------
        df : pd.DataFrame
            DataFrame with multi-tool metrics for all proteins
        """
        metrics_list = []
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        for protein_name, data in sorted(proteins.items()):
            stats_df = data['stats']
            
            row_dict = {'Protein': protein_name}
            
            for tool in tools:
                tool_row = stats_df[stats_df['Tool'] == tool]
                
                if not tool_row.empty:
                    tool_row = tool_row.iloc[0]
                    
                    # Parse mean disorder
                    mean_disorder = float(tool_row['Mean Disorder'])
                    
                    # Parse disorder percentage
                    disorder_pct = tool_row['Disorder Percentage']
                    if isinstance(disorder_pct, str):
                        disorder_pct = float(disorder_pct.rstrip('%'))
                    
                    # Parse std dev
                    std_dev = float(tool_row['Std Dev'])
                    
                    row_dict[f'{tool}_Mean'] = mean_disorder
                    row_dict[f'{tool}_Pct'] = disorder_pct
                    row_dict[f'{tool}_StdDev'] = std_dev
            
            # Calculate consensus (mean across tools)
            tool_means = [row_dict.get(f'{tool}_Mean') for tool in tools if f'{tool}_Mean' in row_dict]
            if tool_means:
                row_dict['Consensus_Mean'] = np.mean(tool_means)
                row_dict['Tool_Agreement'] = 1 - np.std(tool_means)  # Higher = more agreement
            
            metrics_list.append(row_dict)
        
        return pd.DataFrame(metrics_list)
    
    def rank_proteins(self, metrics_df: pd.DataFrame) -> pd.DataFrame:
        """
        Rank proteins by disorder level.
        
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
        print("PROTEIN DISORDER RANKING")
        print("="*100)
        print("\nRanked from MOST to LEAST Intrinsic Disorder:\n")
        print(ranked_df.to_string(index=False))
        
        # Add summary statistics
        print("\n" + "="*100)
        print("SUMMARY STATISTICS")
        print("="*100)
        print(f"\nHighest disorder:  {ranked_df.iloc[0]['Protein']} ({ranked_df.iloc[0]['Disorder %']:.1f}%)")
        print(f"Lowest disorder:   {ranked_df.iloc[-1]['Protein']} ({ranked_df.iloc[-1]['Disorder %']:.1f}%)")
        print(f"Mean disorder across all proteins: {ranked_df['Disorder %'].mean():.1f}%")
        print(f"Std dev across proteins: {ranked_df['Disorder %'].std():.1f}%")
    
    def plot_disorder_ranking_bars(self, ranked_df: pd.DataFrame) -> None:
        """Create bar plot ranking proteins by disorder."""
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Create smoother color gradient from red (high disorder) to green (low disorder)
        # Normalize disorder values to 0-1 range for colormap
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
        ax.set_title('Protein Fragment Ranking by Intrinsic Disorder\n(AlphaFold2 pLDDT-based prediction)',
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add value labels on bars
        for i, (idx, row) in enumerate(ranked_df.iterrows()):
            ax.text(row['Disorder %'] + 0.5, i, f"{row['Disorder %']:.1f}%", 
                   va='center', fontweight='bold', fontsize=10)
        
        # Add grid
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(0, ranked_df['Disorder %'].max() * 1.15)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'protein_disorder_ranking.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'protein_disorder_ranking.png'}")
        plt.close()
    
    def plot_disorder_vs_residues(self, ranked_df: pd.DataFrame) -> None:
        """Create scatter plot: disorder % vs number of residues."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
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
        ax.set_title('Protein Disorder vs Fragment Size',
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add colorbar with proper label
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Disorder (%)', fontsize=11, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'disorder_vs_residues.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'disorder_vs_residues.png'}")
        plt.close()
    
    def plot_comparison_heatmap(self, ranked_df: pd.DataFrame) -> None:
        """Create heatmap comparing key metrics."""
        # Prepare data for heatmap
        heatmap_data = ranked_df[['Protein', 'Disorder %', 'Mean Disorder', 'Std Dev']].copy()
        heatmap_data.set_index('Protein', inplace=True)
        
        # Normalize columns for better visualization
        heatmap_normalized = (heatmap_data - heatmap_data.mean()) / heatmap_data.std()
        
        fig, ax = plt.subplots(figsize=(10, 12))
        
        sns.heatmap(heatmap_normalized.T, annot=True, fmt='.2f', cmap='RdYlGn_r',
                   cbar_kws={'label': 'Standardized Score'}, ax=ax,
                   linewidths=1, linecolor='black')
        
        ax.set_xlabel('Protein Fragment', fontsize=12, fontweight='bold')
        ax.set_ylabel('Disorder Metric', fontsize=12, fontweight='bold')
        ax.set_title('Normalized Disorder Metrics Heatmap\n(Standardized Z-scores)',
                    fontsize=14, fontweight='bold', pad=20)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'disorder_metrics_heatmap.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'disorder_metrics_heatmap.png'}")
        plt.close()
    
    def plot_distribution_comparison(self, ranked_df: pd.DataFrame) -> None:
        """Create distribution comparison plot."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Disorder Metrics Distribution Across All Proteins',
                    fontsize=14, fontweight='bold', y=1.00)
        
        # 1. Disorder % histogram
        ax = axes[0, 0]
        ax.hist(ranked_df['Disorder %'], bins=8, edgecolor='black', color='steelblue', alpha=0.7)
        ax.axvline(ranked_df['Disorder %'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {ranked_df['Disorder %'].mean():.1f}%")
        ax.set_xlabel('Disorder %', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Distribution of Disorder Percentage')
        ax.legend()
        ax.grid(alpha=0.3)
        
        # 2. Mean disorder histogram
        ax = axes[0, 1]
        ax.hist(ranked_df['Mean Disorder'], bins=8, edgecolor='black', color='coral', alpha=0.7)
        ax.axvline(ranked_df['Mean Disorder'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {ranked_df['Mean Disorder'].mean():.3f}")
        ax.set_xlabel('Mean Disorder Score', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Distribution of Mean Disorder Scores')
        ax.legend()
        ax.grid(alpha=0.3)
        
        # 3. Std Dev histogram
        ax = axes[1, 0]
        ax.hist(ranked_df['Std Dev'], bins=8, edgecolor='black', color='lightgreen', alpha=0.7)
        ax.axvline(ranked_df['Std Dev'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {ranked_df['Std Dev'].mean():.3f}")
        ax.set_xlabel('Standard Deviation', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Distribution of Disorder Variability (Std Dev)')
        ax.legend()
        ax.grid(alpha=0.3)
        
        # 4. Residue count histogram
        ax = axes[1, 1]
        ax.hist(ranked_df['Total Residues'], bins=8, edgecolor='black', color='lightcoral', alpha=0.7)
        ax.axvline(ranked_df['Total Residues'].mean(), color='red', linestyle='--', linewidth=2, label=f"Mean: {ranked_df['Total Residues'].mean():.0f}")
        ax.set_xlabel('Total Residues', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Distribution of Fragment Sizes')
        ax.legend()
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'disorder_distributions.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'disorder_distributions.png'}")
        plt.close()
    
    def plot_multi_tool_comparison(self, multi_tool_df: pd.DataFrame) -> None:
        """
        Create comparison plot showing all 4 tools for each protein.
        
        Parameters
        ----------
        multi_tool_df : pd.DataFrame
            DataFrame with multi-tool metrics
        """
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        # Create bar plot comparing mean disorder across tools
        fig, ax = plt.subplots(figsize=(16, 8))
        
        x = np.arange(len(multi_tool_df))
        width = 0.2
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        
        for i, tool in enumerate(tools):
            tool_key = f'{tool}_Mean'
            if tool_key in multi_tool_df.columns:
                ax.bar(x + i*width, multi_tool_df[tool_key], width, 
                      label=tool, color=colors[i], edgecolor='black', linewidth=1)
        
        ax.set_xlabel('Protein Fragment', fontsize=12, fontweight='bold')
        ax.set_ylabel('Mean Disorder Score', fontsize=12, fontweight='bold')
        ax.set_title('Multi-Tool Disorder Prediction Comparison\nAcross All Protein Fragments',
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x + width * 1.5)
        ax.set_xticklabels(multi_tool_df['Protein'], rotation=45, ha='right', fontweight='bold')
        ax.legend(loc='upper right', fontsize=11, framealpha=0.95)
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'multi_tool_comparison_bars.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'multi_tool_comparison_bars.png'}")
        plt.close()
    
    def plot_tool_agreement_heatmap(self, multi_tool_df: pd.DataFrame) -> None:
        """
        Create heatmap showing all tool predictions for each protein.
        
        Parameters
        ----------
        multi_tool_df : pd.DataFrame
            DataFrame with multi-tool metrics
        """
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        # Prepare data for heatmap
        heatmap_data = []
        for _, row in multi_tool_df.iterrows():
            protein_data = [row[f'{tool}_Mean'] for tool in tools]
            heatmap_data.append(protein_data)
        
        heatmap_array = np.array(heatmap_data)
        
        fig, ax = plt.subplots(figsize=(10, 12))
        
        # Create heatmap
        im = ax.imshow(heatmap_array.T, cmap='RdYlGn_r', aspect='auto')
        
        # Set ticks and labels
        ax.set_xticks(np.arange(len(multi_tool_df)))
        ax.set_yticks(np.arange(len(tools)))
        ax.set_xticklabels(multi_tool_df['Protein'], rotation=45, ha='right', fontweight='bold')
        ax.set_yticklabels(tools, fontweight='bold')
        
        # Add text annotations
        for i in range(len(tools)):
            for j in range(len(multi_tool_df)):
                text = ax.text(j, i, f'{heatmap_array[j, i]:.3f}',
                             ha="center", va="center", color="black", fontweight='bold', fontsize=9)
        
        ax.set_xlabel('Protein Fragment', fontsize=12, fontweight='bold')
        ax.set_ylabel('Prediction Tool', fontsize=12, fontweight='bold')
        ax.set_title('Multi-Tool Disorder Prediction Heatmap\n(Mean Disorder Scores)',
                    fontsize=14, fontweight='bold', pad=20)
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Mean Disorder Score', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'multi_tool_heatmap.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'multi_tool_heatmap.png'}")
        plt.close()
    
    def plot_tool_agreement_scores(self, multi_tool_df: pd.DataFrame) -> None:
        """
        Create plot showing tool agreement scores (consensus) for each protein.
        
        Parameters
        ----------
        multi_tool_df : pd.DataFrame
            DataFrame with multi-tool metrics
        """
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Sort by agreement score
        sorted_df = multi_tool_df.sort_values('Tool_Agreement', ascending=True)
        
        # Create color gradient based on agreement
        colors = plt.cm.RdYlGn(sorted_df['Tool_Agreement'] / sorted_df['Tool_Agreement'].max())
        
        bars = ax.barh(range(len(sorted_df)), sorted_df['Tool_Agreement'], 
                      color=colors, edgecolor='black', linewidth=1.5)
        
        ax.set_yticks(range(len(sorted_df)))
        ax.set_yticklabels(sorted_df['Protein'], fontsize=11, fontweight='bold')
        ax.set_xlabel('Tool Agreement Score (1 - Tool Std Dev)\nHigher = Better Consensus',
                     fontsize=12, fontweight='bold')
        ax.set_title('Inter-Tool Agreement: Consensus Quality Across 4 Prediction Methods\n(Measures how well AlphaFold2, IUPred3, DISOPRED3, and SPOT-Disorder agree)',
                    fontsize=13, fontweight='bold', pad=20)
        
        # Add value labels
        for i, (idx, row) in enumerate(sorted_df.iterrows()):
            ax.text(row['Tool_Agreement'] + 0.002, i, f"{row['Tool_Agreement']:.4f}", 
                   va='center', fontweight='bold', fontsize=10)
        
        ax.grid(axis='x', alpha=0.3)
        ax.set_xlim(min(0.85, sorted_df['Tool_Agreement'].min() - 0.01), 1.0)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'tool_agreement_scores.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {self.plots_dir / 'tool_agreement_scores.png'}")
        plt.close()
    
    def save_ranking_table(self, ranked_df: pd.DataFrame) -> None:
        """Save ranking table as CSV."""
        csv_path = self.data_dir / 'protein_disorder_ranking.csv'
        ranked_df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")
        
        # Also save as detailed table with all metrics
        detailed_path = self.data_dir / 'protein_disorder_detailed.csv'
        ranked_df.to_csv(detailed_path, index=False)
        print(f"Saved: {detailed_path}")
    
    def save_multi_tool_comparison(self, multi_tool_df: pd.DataFrame) -> None:
        """Save multi-tool comparison data as CSV."""
        csv_path = self.data_dir / 'protein_multi_tool_comparison.csv'
        multi_tool_df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")
        
        # Create and save summary comparison table
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        summary_list = []
        
        for _, row in multi_tool_df.iterrows():
            summary_list.append({
                'Protein': row['Protein'],
                'AF2_Mean': row['AlphaFold2 (pLDDT)_Mean'],
                'IUPred3_Mean': row['IUPred3_Mean'],
                'DISOPRED3_Mean': row['DISOPRED3_Mean'],
                'SPOT_Mean': row['SPOT-Disorder_Mean'],
                'Consensus_Mean': row['Consensus_Mean'],
                'Tool_Agreement': row['Tool_Agreement'],
            })
        
        summary_df = pd.DataFrame(summary_list)
        summary_path = self.data_dir / 'protein_tool_consensus.csv'
        summary_df.to_csv(summary_path, index=False)
        print(f"Saved: {summary_path}")
    
    def print_multi_tool_summary(self, multi_tool_df: pd.DataFrame) -> None:
        """Print multi-tool analysis summary to console."""
        print("\n" + "="*120)
        print("MULTI-TOOL DISORDER PREDICTION COMPARISON")
        print("="*120)
        
        tools = ['AlphaFold2 (pLDDT)', 'IUPred3', 'DISOPRED3', 'SPOT-Disorder']
        
        print("\n📊 Mean Disorder Scores by Tool (Across All Proteins):\n")
        
        tool_stats = {}
        for tool in tools:
            tool_key = f'{tool}_Mean'
            if tool_key in multi_tool_df.columns:
                values = multi_tool_df[tool_key]
                tool_stats[tool] = {
                    'mean': values.mean(),
                    'std': values.std(),
                    'min': values.min(),
                    'max': values.max(),
                }
        
        for tool, stats in tool_stats.items():
            print(f"  {tool:25s} | Mean: {stats['mean']:.4f} | Std: {stats['std']:.4f} | "
                  f"Range: [{stats['min']:.4f}, {stats['max']:.4f}]")
        
        print("\n📈 Tool Agreement Across Proteins:\n")
        
        agreement_sorted = multi_tool_df.sort_values('Tool_Agreement', ascending=False)
        print(f"  Highest agreement: {agreement_sorted.iloc[0]['Protein']:20s} "
              f"(Agreement: {agreement_sorted.iloc[0]['Tool_Agreement']:.4f})")
        print(f"  Lowest agreement:  {agreement_sorted.iloc[-1]['Protein']:20s} "
              f"(Agreement: {agreement_sorted.iloc[-1]['Tool_Agreement']:.4f})")
        print(f"  Mean agreement across all proteins: {multi_tool_df['Tool_Agreement'].mean():.4f}")
        
        print("\n" + "="*120)
    
    def run_comparison(self, disorder_analysis_dir: str = "disorder_analysis") -> None:
        """Run complete comparison analysis."""
        
        print("\n" + "="*100)
        print("MULTI-PROTEIN DISORDER COMPARISON")
        print("="*100)
        
        # Find protein data
        print("\nFinding protein analysis data...")
        proteins = self.find_protein_data(disorder_analysis_dir)
        
        if not proteins:
            print("\n✗ No protein data found")
            return
        
        print(f"✓ Found {len(proteins)} proteins")
        
        # Extract metrics
        print("\nExtracting disorder metrics...")
        metrics_df = self.extract_disorder_metrics(proteins)
        
        # Extract multi-tool metrics
        print("Extracting multi-tool metrics...")
        multi_tool_df = self.extract_multi_tool_metrics(proteins)
        
        # Rank proteins
        print("Ranking proteins by disorder...")
        ranked_df = self.rank_proteins(metrics_df)
        
        # Print ranking table
        self.print_ranking_table(ranked_df)
        
        # Print multi-tool summary
        self.print_multi_tool_summary(multi_tool_df)
        
        # Generate plots
        print("\n" + "="*100)
        print("GENERATING COMPARISON FIGURES")
        print("="*100)
        
        self.plot_disorder_ranking_bars(ranked_df)
        self.plot_disorder_vs_residues(ranked_df)
        self.plot_comparison_heatmap(ranked_df)
        self.plot_distribution_comparison(ranked_df)
        
        # Multi-tool plots
        print("\nGenerating multi-tool comparison plots...")
        self.plot_multi_tool_comparison(multi_tool_df)
        self.plot_tool_agreement_heatmap(multi_tool_df)
        self.plot_tool_agreement_scores(multi_tool_df)
        
        # Save ranking table
        print("\n" + "="*100)
        print("SAVING DATA")
        print("="*100)
        
        self.save_ranking_table(ranked_df)
        self.save_multi_tool_comparison(multi_tool_df)
        
        # Final summary
        print("\n" + "="*100)
        print("COMPARISON ANALYSIS COMPLETE")
        print("="*100)
        print(f"\n✓ Results saved in: {self.output_dir}")
        print(f"  - Plots: {self.plots_dir}")
        print(f"  - Data: {self.data_dir}")
        
        print("\n📊 Generated Figures:")
        print("   - protein_disorder_ranking.png - Bar chart ranking by disorder")
        print("   - disorder_vs_residues.png - Scatter plot: disorder vs size")
        print("   - disorder_metrics_heatmap.png - Heatmap of normalized metrics")
        print("   - disorder_distributions.png - Distribution of all metrics")
        print("   - multi_tool_comparison_bars.png - Bar chart: all 4 tools vs proteins")
        print("   - multi_tool_heatmap.png - Heatmap of all tool predictions")
        print("   - tool_agreement_scores.png - Inter-tool consensus quality")
        
        print("\n📋 Data Files:")
        print("   - protein_disorder_ranking.csv - Ranking table (AlphaFold2 focus)")
        print("   - protein_disorder_detailed.csv - Detailed metrics")
        print("   - protein_multi_tool_comparison.csv - All tool data for all proteins")
        print("   - protein_tool_consensus.csv - Consensus scores and agreement metrics")


def main():
    parser = argparse.ArgumentParser(
        description='Multi-Protein Disorder Comparison and Ranking',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: analyze disorder_analysis/ and save to protein_comparison/
  python3 03_compare_proteins.py

  # Custom output directory
  python3 03_compare_proteins.py --output_dir my_comparison

  # Custom analysis directory
  python3 03_compare_proteins.py --analysis_dir /path/to/disorder_analysis
        """
    )
    
    parser.add_argument('--analysis_dir', default='disorder_analysis',
                       help='Path to disorder_analysis directory (default: disorder_analysis)')
    parser.add_argument('--output_dir', default='protein_comparison',
                       help='Output directory for comparison results (default: protein_comparison)')
    
    args = parser.parse_args()
    
    # Run comparison
    comparison = ProteinDisorderComparison(output_dir=args.output_dir)
    comparison.run_comparison(disorder_analysis_dir=args.analysis_dir)


if __name__ == '__main__':
    main()
