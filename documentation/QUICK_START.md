# Quick Start Guide

## For Users Who Want to Run the Analysis

### Prerequisites

1. **Python 3.7+**
2. **Dependencies**:
   ```bash
   pip install numpy pandas matplotlib seaborn biopython openpyxl
   ```
3. **AlphaFold3 Output Files** (JSON format):
   - `fold_[range]_full_data_[0-4].json` (5 models)
   - `fold_[range]_summary_confidences_[0-4].json` (confidence scores)

### Directory Structure Setup

Organize your AlphaFold data like this:

```
your_project/
├── protein_folder_1/
│   ├── fold_[range]_full_data_0.json
│   ├── fold_[range]_full_data_1.json
│   ├── fold_[range]_full_data_2.json
│   ├── fold_[range]_full_data_3.json
│   ├── fold_[range]_full_data_4.json
│   ├── fold_[range]_summary_confidences_0.json
│   ├── fold_[range]_summary_confidences_1.json
│   ├── fold_[range]_summary_confidences_2.json
│   ├── fold_[range]_summary_confidences_3.json
│   └── fold_[range]_summary_confidences_4.json
├── protein_folder_2/
│   ├── (same structure)
└── scripts/
    ├── 01_run_analysis.py
    ├── 02_analyze_alphafold_disorder.py
    ├── 02.5_run_all_enhanced.py
    └── 03_compare_proteins.py
```

### Step 1: Install Dependencies

```bash
pip install numpy pandas matplotlib seaborn biopython openpyxl
```

### Step 2: Run Analysis on Single Protein (Optional)

```bash
python3 01_run_analysis.py
```

Then select your protein and model. This is good for testing.

### Step 3: Run Batch Analysis (Recommended)

This analyzes all protein folders automatically:

```bash
python3 02.5_run_all_enhanced.py
```

**Output**: Creates `disorder_analysis/` folder with results for each protein.

### Step 4: Generate Comparison Report

```bash
python3 03_compare_proteins.py
```

**Output**:
- `protein_comparison/plots/` - Ranking charts, scatter plots, multi-tool comparisons
- `protein_comparison/data/` - CSV files with statistics and rankings

### Step 5: Create Excel Files (From Output Data)

The analysis scripts generate CSV files. To convert to Excel publication-ready format:

```python
# Run this Python snippet
import pandas as pd
from pathlib import Path

# Read comparison data
csv_file = Path('protein_comparison/data/protein_disorder_ranking.csv')
df = pd.read_csv(csv_file)

# Create Excel with formatting
output_file = 'Protein_Fragment_Ranking.xlsx'
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='Data', index=False)
```

### Workflow Summary

```
Step 1: Prepare AlphaFold3 JSON files
        ↓
Step 2: Run 02.5_run_all_enhanced.py
        ↓
Step 3: Run 03_compare_proteins.py
        ↓
Step 4: View results in protein_comparison/
        ↓
Step 5: Create publication figures from Excel/CSV data
```

---

## For Users Who Just Want to Explore the Data

### Read-Only Access

You don't need to run anything! The Excel files are ready to use:

#### Location 1: RdRp_EGFPOVA Proteins
```
data/RdRp_EGFPOVA_excel_files/
├── 366_581_plddt.xlsx          # Per-residue pLDDT data
├── [Other fragments]_plddt.xlsx # (12 more files)
├── Protein_Fragment_Ranking_For_Chart.xlsx
├── Multi_Tool_Disorder_Comparison.xlsx
└── Inter_Tool_Agreement_Comparison.xlsx
```

#### Location 2: RdRp_WT Proteins
```
data/RdRp_WT_excel_files/
├── 1_130_plddt.xlsx            # Per-residue pLDDT data
├── [Other fragments]_plddt.xlsx # (12 more files)
├── RdRp_WT_Protein_Fragment_Ranking_For_Chart.xlsx
├── RdRp_WT_Multi_Tool_Disorder_Comparison.xlsx
└── RdRp_WT_Inter_Tool_Agreement_Comparison.xlsx
```

### Exploring the Data

1. **Per-Residue Files** ([Fragment]_plddt.xlsx):
   - Open in Excel/Google Sheets
   - Columns: Residue, pLDDT_Score, High_Confidence_Over_70, Disorder_Threshold_Below_50, Confidence_Classification
   - Color-coded by confidence level
   - Plot: Residue Position (X-axis) vs pLDDT_Score (Y-axis)

2. **Ranking Files** (Protein_Fragment_Ranking_For_Chart.xlsx):
   - Shows all 13 fragments ranked by disorder percentage
   - Includes pLDDT standard deviation
   - Ready for bar chart creation

3. **Comparison Files** (Multi_Tool_Disorder_Comparison.xlsx):
   - Comparison of 4 prediction tools (AlphaFold3, IUPred3, DISOPRED3, SPOT-Disorder)
   - Ready for clustered column chart

4. **Agreement Files** (Inter_Tool_Agreement_Comparison.xlsx):
   - Shows how well prediction tools agree
   - Higher = Better consensus
   - Range: 0.96-0.98 (excellent agreement)

### Creating Publication Figures

#### Figure 1: Fragment Ranking
```
File: Protein_Fragment_Ranking_For_Chart.xlsx
Steps:
  1. Select columns: Protein_Fragment, Intrinsic_Disorder_Percent
  2. Insert → Chart → Bar Chart (Horizontal)
  3. Format: Red for high disorder, green for low
  4. Title: "Protein Fragment Ranking by Intrinsic Disorder"
```

#### Figure 2: Multi-Tool Comparison
```
File: Multi_Tool_Disorder_Comparison.xlsx
Steps:
  1. Select columns: Protein_Fragment, AlphaFold3, IUPred3, DISOPRED3, SPOT-Disorder
  2. Insert → Chart → Column Chart (Clustered)
  3. Format: Different colors for each tool
  4. Title: "Multi-Tool Disorder Prediction Comparison"
```

#### Figure 3: Inter-Tool Agreement
```
File: Inter_Tool_Agreement_Comparison.xlsx
Steps:
  1. Select columns: Protein_Fragment, Tool_Agreement_Score
  2. Insert → Chart → Bar Chart (Horizontal)
  3. Format: Green color
  4. Title: "Inter-Tool Agreement: Consensus Quality"
```

#### Figure 4: Per-Residue Detail
```
File: [Fragment]_plddt.xlsx (e.g., 366_581_plddt.xlsx)
Steps:
  1. Select columns: Residue, pLDDT_Score
  2. Insert → Chart → Line Chart
  3. Add horizontal threshold lines at 70 (green) and 50 (orange)
  4. Title: "AlphaFold3 Confidence Across 366-581 Fragment"
  5. Annotation: Highlight disordered regions (below 50)
```

---

## Understanding the Results

### Key Metrics Explained

**Intrinsic Disorder (%)**
- Percentage of residues with pLDDT < 50
- Higher = More disordered
- 366-581: 31.5% (RdRp_EGFPOVA) vs 43.1% (RdRp_WT)

**pLDDT Score (0-100)**
- AlphaFold3 confidence metric
- ≥70: High confidence (ordered)
- <50: Low confidence (disordered)

**Standard Deviation**
- Reflects variability of pLDDT across residues
- High = Mixed ordered/disordered regions
- Low = Consistent predictions
- HIGH is NORMAL for disordered proteins!

**Inter-Tool Agreement (0-1)**
- How well 4 tools agree on disorder prediction
- 0.97-1.00: Excellent agreement
- 0.95-0.97: Good agreement
- Lower values = More uncertain predictions

### Key Findings at a Glance

| Metric | RdRp_EGFPOVA | RdRp_WT | Significance |
|--------|--------|-----------|-------------|
| **366-581 Disorder** | 31.5% | 43.1% | RdRp_WT is MORE disordered |
| **Overall Mean Disorder** | 3-31% | 0-43% | High variability across fragments |
| **Tool Agreement** | 0.957-0.975 | 0.957-0.975 | Excellent consensus |
| **Most Ordered** | OG_Protein (3.2%) | RdRp_WT_1_581 (0%) | Full-length is most structured |

---

## Troubleshooting

### "Module not found" errors
```bash
pip install --upgrade numpy pandas matplotlib seaborn biopython openpyxl
```

### Excel files won't open
- Use Excel 2016+ or LibreOffice Calc
- Check file permissions: `chmod +r *.xlsx`

### Analysis running too slowly
- Run on subset of fragments first
- Check system RAM (recommended: 8GB+)
- Reduce number of AlphaFold models analyzed

### Need more help?
- Check METHODS.md for detailed methodology
- Review README.md for complete documentation
- Examine Excel files directly - data is self-documented

---

## Next Steps

1. **Explore the Data**: Open Excel files and browse per-residue data
2. **Understand the Methods**: Read METHODS.md for detailed explanation
3. **Create Figures**: Use ranking/comparison Excel files for publication plots
4. **Run Full Pipeline**: Execute all scripts on your own data
5. **Integrate**: Use scripts for your own protein disorder analysis

Enjoy exploring your protein disorder data! 🔬
