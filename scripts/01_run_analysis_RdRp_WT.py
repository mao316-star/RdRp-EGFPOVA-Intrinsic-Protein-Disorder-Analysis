#!/usr/bin/env python3
"""
Interactive wrapper for AlphaFold Disorder Analysis - WT Proteins
Analyzes wild-type protein variants from WT_AF_info
Provides user-friendly interface with installation assistance
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def check_dependencies():
    """Check if required dependencies are installed."""
    required = ['numpy', 'pandas', 'matplotlib', 'seaborn']
    missing = []
    
    for package in required:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    return missing


def install_dependencies(include_optional=True):
    """Install required and optional dependencies."""
    print("\n" + "="*70)
    print("Installing Dependencies")
    print("="*70)
    
    # Core dependencies
    print("\n✓ Installing core dependencies...")
    core_packages = ['numpy', 'pandas', 'matplotlib', 'seaborn']
    subprocess.run([sys.executable, '-m', 'pip', 'install'] + core_packages, 
                   check=False, capture_output=True)
    
    # Optional dependencies
    if include_optional:
        print("✓ Installing optional tools (BioPython, IUPred3)...")
        optional_packages = ['biopython', 'iupred-fast']
        subprocess.run([sys.executable, '-m', 'pip', 'install'] + optional_packages,
                       check=False, capture_output=True)
        print("  Note: DISOPRED3 and SPOT-Disorder require manual installation")
    
    print("\n✓ Dependencies installed successfully!")


def find_af_data():
    """Find available AlphaFold data files for wild-type proteins."""
    # Try multiple paths to find RdRp_WT structures
    possible_paths = [
        Path('data/alphafold_raw_data/RdRp_WT_structures'),
        Path('../data/alphafold_raw_data/RdRp_WT_structures'),
        Path('../../RdRp-Disorder-Analysis/data/alphafold_raw_data/RdRp_WT_structures')
    ]
    
    af_dir = None
    for path in possible_paths:
        if path.exists() and any(p.name.startswith('wt_') for p in path.iterdir() if p.is_dir()):
            af_dir = path
            break
    
    if af_dir is None:
        print(f"✗ RdRp_WT structures directory not found")
        print(f"  Expected location: data/alphafold_raw_data/RdRp_WT_structures")
        return None
    
    print(f"✓ Found RdRp_WT AlphaFold data at: {af_dir}")
    
    # Find all protein folders with AlphaFold data
    proteins = {}
    
    # Subdirectories with protein variants (wt_*)
    subdirs = sorted([d for d in af_dir.iterdir() 
                     if d.is_dir() and d.name.startswith('wt_') and d.name not in ['msas', 'templates', '.DS_Store']])
    
    for subdir in subdirs:
        # Find full_data JSON files in each protein directory
        data_files = sorted(subdir.glob('fold_*_full_data_*.json'))
        if data_files:
            proteins[subdir.name] = data_files
    
    if not proteins:
        print(f"✗ No AlphaFold data files found in {af_dir}")
        return None
    
    return proteins


def select_proteins(proteins_dict):
    """Interactive protein selection for batch analysis."""
    print("\n" + "="*70)
    print("Available WT Protein Variants")
    print("="*70)
    
    protein_list = list(proteins_dict.items())
    
    for i, (protein_name, data_files) in enumerate(protein_list):
        print(f"  [{i}] {protein_name}: {len(data_files)} models")
    
    print(f"  [A] Analyze ALL proteins")
    print(f"  [Q] Quit")
    
    print()
    while True:
        choice = input("Select protein to analyze (enter number, A for all, or Q to quit): ").strip().upper()
        
        if choice == 'Q':
            return None
        elif choice == 'A':
            return protein_list  # Return all proteins
        else:
            try:
                idx = int(choice)
                if 0 <= idx < len(protein_list):
                    return [protein_list[idx]]  # Return as list for consistency
            except ValueError:
                pass
        
        print("Invalid selection. Please try again.")


def select_model(data_files):
    """Interactive model selection."""
    print("\n" + "="*70)
    print("Available AlphaFold Models")
    print("="*70)
    
    for i, f in enumerate(data_files):
        model_num = f.stem.split('_')[-1]
        print(f"  [{i}] Model {model_num}: {f.name}")
    
    print()
    while True:
        try:
            choice = input("Select model to analyze (enter number): ").strip()
            idx = int(choice)
            if 0 <= idx < len(data_files):
                return data_files[idx]
        except ValueError:
            pass
        print("Invalid selection. Please try again.")


def find_reference_pdb():
    """Find reference PDB file (optional for WT analysis)."""
    print("\n" + "="*70)
    print("Reference PDB File (Optional)")
    print("="*70)
    
    # WT proteins may not have reference PDB files
    print("\n⚠ Reference PDB not required for WT analysis (optional - disorder tools will be estimated)")
    return None


def run_analysis(af_data_file, confidence_file, protein_name, reference_pdb=None):
    """Run the AlphaFold disorder analysis for a single protein."""
    print("\n" + "="*70)
    print(f"Running AlphaFold Disorder Analysis for {protein_name}")
    print("="*70)
    
    # Create protein-specific output directory
    output_dir = f'disorder_analysis_WT/{protein_name}'
    
    # Get the path to the analysis script (in parent AF_Disorder directory)
    script_path = Path(__file__).parent / '02_analyze_alphafold_disorder.py'
    
    if not script_path.exists():
        print(f"✗ Analysis script not found at {script_path}")
        return False
    
    cmd = [
        sys.executable,
        str(script_path),
        '--af_data', str(af_data_file),
        '--confidence', str(confidence_file),
        '--output_dir', output_dir
    ]
    
    if reference_pdb:
        cmd.extend(['--reference', str(reference_pdb)])
    
    print(f"\nCommand: {' '.join(cmd)}\n")
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Analysis failed with error: {e}")
        return False


def print_results_summary(results_dict=None):
    """Print summary of results."""
    output_dir = Path('protein_analysis')
    
    if not output_dir.exists():
        return
    
    print("\n" + "="*70)
    print("Analysis Results Summary")
    print("="*70)
    
    if results_dict:
        print("\nWT Protein Analysis Status:")
        for protein_name, success in sorted(results_dict.items()):
            status = "✓ SUCCESS" if success else "✗ FAILED"
            print(f"  {status}: {protein_name}")
    
    # List all results by protein
    protein_dirs = sorted([d for d in output_dir.iterdir() if d.is_dir()])
    
    if protein_dirs:
        print("\nGenerated Results by Protein:")
        for pdir in protein_dirs:
            print(f"\n  📁 {pdir.name}/")
            
            plots_dir = pdir / 'plots'
            if plots_dir.exists():
                plots = list(plots_dir.glob('*.png'))
                if plots:
                    print(f"     📊 Figures ({len(plots)}):")
                    for plot in sorted(plots):
                        print(f"        - {plot.name}")
            
            data_dir = pdir / 'data'
            if data_dir.exists():
                csv_files = list(data_dir.glob('*.csv'))
                if csv_files:
                    print(f"     📋 Data Files ({len(csv_files)}):")
                    for csv in sorted(csv_files):
                        print(f"        - {csv.name}")
    
    print(f"\n✓ All results saved in: protein_analysis/")
    print(f"  Total WT proteins analyzed: {len(protein_dirs)}")


def main():
    parser = argparse.ArgumentParser(
        description='Interactive AlphaFold Disorder Analysis - WT Proteins',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (recommended)
  python3 01_run_analysis_WT.py

  # Analyze specific WT protein
  python3 01_run_analysis_WT.py --protein wt_366_581

  # Analyze all WT proteins
  python3 01_run_analysis_WT.py --all

  # Skip dependency check
  python3 01_run_analysis_WT.py --skip-check
        """
    )
    
    parser.add_argument('--protein', type=str, default=None,
                       help='Specific WT protein to analyze (e.g., wt_366_581)')
    parser.add_argument('--model', type=int, default=None,
                       help='Model number to analyze (0-4)')
    parser.add_argument('--all', action='store_true',
                       help='Analyze all WT proteins')
    parser.add_argument('--skip-check', action='store_true',
                       help='Skip dependency check')
    parser.add_argument('--skip-install', action='store_true',
                       help='Skip dependency installation')
    
    args = parser.parse_args()
    
    # Change to repository root directory
    script_dir = Path(__file__).parent.resolve()
    repo_root = script_dir.parent
    os.chdir(repo_root)
    
    print("\n" + "="*70)
    print("AlphaFold Disorder Analysis - RdRp_WT Proteins Setup")
    print("="*70)
    print(f"Repository root: {os.getcwd()}")
    
    # Check dependencies
    if not args.skip_check:
        missing = check_dependencies()
        if missing:
            print(f"\n⚠ Missing dependencies: {', '.join(missing)}")
            if not args.skip_install:
                response = input("\nInstall dependencies now? (y/n): ").strip().lower()
                if response == 'y':
                    install_dependencies(include_optional=True)
                else:
                    print("Skipping dependency installation")
    
    # Find AF data for all WT proteins
    print("\n" + "="*70)
    print("Locating AlphaFold Data")
    print("="*70)
    
    proteins_dict = find_af_data()
    if not proteins_dict:
        print("\n✗ Cannot find AlphaFold data files")
        sys.exit(1)
    
    print(f"\n✓ Found {len(proteins_dict)} WT proteins:")
    for protein_name in sorted(proteins_dict.keys()):
        print(f"   - {protein_name}")
    
    # Determine which proteins to analyze
    if args.all:
        selected_proteins = list(proteins_dict.items())
    elif args.protein:
        if args.protein in proteins_dict:
            selected_proteins = [(args.protein, proteins_dict[args.protein])]
        else:
            print(f"✗ WT Protein '{args.protein}' not found")
            sys.exit(1)
    else:
        selected = select_proteins(proteins_dict)
        if selected is None:
            print("Analysis cancelled")
            sys.exit(0)
        selected_proteins = selected
    
    # Find reference PDB (optional for WT)
    reference_pdb = find_reference_pdb()
    
    # Create main output directory
    output_dir = Path('disorder_analysis_WT')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run analysis for each selected protein
    results = {}
    total_proteins = len(selected_proteins)
    
    for idx, (protein_name, data_files) in enumerate(selected_proteins, 1):
        print(f"\n\n{'='*70}")
        print(f"WT Protein {idx}/{total_proteins}: {protein_name}")
        print(f"{'='*70}")
        
        # Select model for this protein
        if args.model is not None:
            if 0 <= args.model < len(data_files):
                selected_file = data_files[args.model]
            else:
                print(f"✗ Invalid model number: {args.model}")
                results[protein_name] = False
                continue
        elif args.all:
            # When analyzing all, use model 0 by default
            selected_file = data_files[0]
            print(f"  Using model 0 (default for batch analysis)")
        else:
            selected_file = select_model(data_files)
        
        model_num = selected_file.stem.split('_')[-1]
        
        # Find matching confidence file
        confidence_pattern = selected_file.stem.replace('full_data', 'summary_confidences')
        confidence_file = selected_file.parent / f'{confidence_pattern}.json'
        
        if not confidence_file.exists():
            print(f"✗ Confidence file not found: {confidence_file}")
            results[protein_name] = False
            continue
        
        # Run analysis
        success = run_analysis(selected_file, confidence_file, protein_name, reference_pdb)
        results[protein_name] = success
        
        if success:
            print(f"\n✓ Analysis for {protein_name} completed successfully!")
    
    # Print summary
    print_results_summary(results)
    
    # Check if all analyses succeeded
    all_success = all(results.values())
    if all_success:
        print("\n✓ All WT protein analyses completed successfully!")
    else:
        failed = [name for name, success in results.items() if not success]
        print(f"\n⚠ {len(failed)} WT protein(s) failed: {', '.join(failed)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
