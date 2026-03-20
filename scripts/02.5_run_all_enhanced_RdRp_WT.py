#!/usr/bin/env python3
"""
Run enhanced disorder analysis on all WT (wild-type) proteins sequentially.
This ensures clean execution without hanging.

WT Proteins:
- Automatically discovers all wt_* folders in the current directory
- Analyzes each one with AlphaFold disorder prediction tools
"""

import subprocess
import sys
from pathlib import Path

# Dynamically find all WT protein folders
current_dir = Path.cwd()
proteins = sorted([d.name for d in af_info_dir.iterdir() 
                   if d.is_dir() and d.name.startswith('wt_') and d.name != 'WT-disorder analysis'])

print("="*80)
print("ENHANCED ALPHAFOLD DISORDER ANALYSIS - ALL WT PROTEINS")
print("="*80)

successful = []
failed = []

for i, protein in enumerate(proteins, 1):
    print(f"\n[{i}/{len(proteins)}] Processing {protein}...")
    
    # Find data files in current directory
    protein_dir = Path(protein)
    
    if not protein_dir.exists():
        print(f"  ✗ Directory not found for {protein}")
        failed.append(protein)
        continue
    
    # Find full_data JSON files
    af_data_file = list(protein_dir.glob('fold_*_full_data_0.json'))
    conf_file = list(protein_dir.glob('fold_*_summary_confidences_0.json'))
    
    if not af_data_file or not conf_file:
        print(f"  ✗ Missing data files for {protein}")
        failed.append(protein)
        continue
    
    af_data_file = str(af_data_file[0])
    conf_file = str(conf_file[0])
    output_dir = f"disorder_analysis_WT/{protein}"
    
    # Get the path to the analysis script (in parent AF_Disorder directory)
    script_path = Path(__file__).parent.parent / 'AF_Disorder' / '02_analyze_alphafold_disorder.py'
    
    if not script_path.exists():
        print(f"  ✗ Analysis script not found at {script_path}")
        failed.append(protein)
        continue
    
    # Run analysis
    cmd = [
        "python3",
        str(script_path),
        "--af_data", af_data_file,
        "--confidence", conf_file,
        "--output_dir", output_dir
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print(f"  ✓ Successfully analyzed {protein}")
            successful.append(protein)
        else:
            print(f"  ✗ Error analyzing {protein}")
            print(f"    Error: {result.stderr[:200]}")
            failed.append(protein)
    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout analyzing {protein}")
        failed.append(protein)
    except Exception as e:
        print(f"  ✗ Exception: {e}")
        failed.append(protein)

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Successfully analyzed: {len(successful)}/{len(proteins)}")
print(f"  {', '.join(successful)}")
if failed:
    print(f"\nFailed: {len(failed)}")
    print(f"  {', '.join(failed)}")

sys.exit(0 if len(failed) == 0 else 1)
