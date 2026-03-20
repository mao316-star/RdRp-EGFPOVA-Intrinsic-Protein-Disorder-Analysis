#!/usr/bin/env python3
"""
Run enhanced disorder analysis on all 12 proteins sequentially.
This ensures clean execution without hanging.
"""

import subprocess
import sys
from pathlib import Path

# List of all 13 proteins
proteins = [
    "OG_Protein",
    "del366_581", 
    "1_130",
    "1_276",
    "1_365",
    "1_581",
    "131_249",
    "250_932",
    "251_365",
    "277_932",
    "366_581",
    "582_932",
    "del251_365",
]

af_info_dir = Path("data/alphafold_raw_data/RdRp_EGFPOVA_structures")

print("="*80)
print("ENHANCED ALPHAFOLD DISORDER ANALYSIS - ALL PROTEINS")
print("="*80)

successful = []
failed = []

for i, protein in enumerate(proteins, 1):
    print(f"\n[{i}/{len(proteins)}] Processing {protein}...")
    
    # Find data files
    af_data_file = list(af_info_dir.glob(f"{protein}/*full_data_0.json"))
    conf_file = list(af_info_dir.glob(f"{protein}/*summary_confidences_0.json"))
    
    if not af_data_file or not conf_file:
        print(f"  ✗ Missing data files for {protein}")
        failed.append(protein)
        continue
    
    af_data_file = str(af_data_file[0])
    conf_file = str(conf_file[0])
    output_dir = f"disorder_analysis/{protein}"
    
    # Run analysis
    cmd = [
        "python3",
        "AF_Disorder/02_analyze_alphafold_disorder.py",
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
print(f"Successfully analyzed: {len(successful)}/{ len(proteins)}")
print(f"  {', '.join(successful)}")
if failed:
    print(f"\nFailed: {len(failed)}")
    print(f"  {', '.join(failed)}")

sys.exit(0 if len(failed) == 0 else 1)
