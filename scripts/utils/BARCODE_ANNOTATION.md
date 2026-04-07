# Barcode Annotation System

This directory contains tools for annotating cells by experimental conditions encoded in CapWGS combinatorial barcodes.

## Architecture

The system uses a **layered design** to support both interactive analysis and pipeline automation:

```
┌─────────────────────────────────────────────────────────┐
│  HIGH-LEVEL APIs (Use-case specific)                    │
├─────────────────────────────────────────────────────────┤
│  barcode_annotation.py   │  assign_conditions.py        │
│  (for notebooks)         │  (for pipeline/CLI)          │
├─────────────────────────────────────────────────────────┤
│  CORE LOGIC (Shared utilities)                          │
├─────────────────────────────────────────────────────────┤
│  barcode_core.py - Constants, plate loading, lookup     │
│                    building, condition assignment        │
└─────────────────────────────────────────────────────────┘
```

## Files

### `barcode_core.py` (Core Utilities)
Shared logic used by both high-level APIs:
- **Constants**: Barcode position slices, column mappings
- **Plate loading**: CSV → dict or DataFrame
- **Well parsing**: Flexible well specification (rows, columns, specific wells)
- **Lookup building**: Conditions → barcode subsequence mapping
- **Assignment**: Cell barcodes → condition labels

### `barcode_annotation.py` (Interactive Analysis)
DataFrame-based API for Jupyter notebooks:
- Input: pandas DataFrame with cell barcodes as index
- Output: DataFrame with new condition column
- Flexible well specification: `'A'` (rows), `1` (columns), `'A1'` (wells)
- Works with any barcode position (A/B/C/D)

### `../CapWGS/assign_conditions.py` (Pipeline Automation)
CLI tool for dashboard generation and automation:
- Input: JSON config + barcode map CSV + cell list TXT
- Output: CSV with barcode-condition assignments
- Integrates with HTML dashboard generation
- Warnings instead of exceptions for robustness

## Usage

### Interactive Analysis (Jupyter Notebooks)

```python
import pandas as pd
import sys
sys.path.append('../scripts/utils')
from barcode_annotation import annotate_barcodes

# Load barcode plate map
bcs = pd.read_csv('../barcodes/map.csv', index_col=0)
bcs.columns = bcs.columns.astype(int)

# Load cell data
df = pd.read_csv('../results/sample/compiled_qc_metrics.csv', index_col=0)

# Annotate by rows in position A
df = annotate_barcodes(
    df, bcs,
    column_name='treatment',
    conditions={
        'Control': ['A', 'B'],      # Rows A-B (all position A columns: 1-3)
        'Drug_X': ['C', 'D'],       # Rows C-D
        'Drug_Y': ['E', 'F', 'G', 'H']
    },
    barcode_set='A'
)

# Annotate by specific wells in position D
df = annotate_barcodes(
    df, bcs,
    column_name='replicate',
    conditions={
        'Rep1': ['A10', 'B10'],     # Specific wells
        'Rep2': ['A11', 'B11'],
        'Rep3': ['A12', 'B12']
    },
    barcode_set='D'
)

# Annotate by columns
df = annotate_barcodes(
    df, bcs,
    column_name='batch',
    conditions={
        'Batch1': [1],              # Column 1 (all rows)
        'Batch2': [2],
        'Batch3': [3]
    },
    barcode_set='A'
)
```

### Pipeline/CLI Usage

```bash
# Create conditions config JSON
cat > conditions.json << 'EOF'
{
    "Control": ["A1", "A2", "B1", "B2"],
    "Treatment": ["C1", "C2", "D1", "D2"]
}
EOF

# Run assignment
python3 scripts/CapWGS/assign_conditions.py \
    conditions.json \
    barcodes/map.csv \
    results/sample/real_cells.txt \
    --output results/sample/condition_assignments.csv

# Use from Python (for dashboard generation)
from assign_conditions import assign_conditions_to_barcodes
assignments = assign_conditions_to_barcodes(config, 'barcodes/map.csv', cell_list)
```

## Barcode Structure

Cell barcodes are 45 bases with 4 positions separated by constant sequences:

```
{bcD}AGGA{bcC}ACTC{bcB}AAGG{bcA}A
0-7  8-11 12-19 20-23 24-31 32-35 36-43 44
```

Each 8bp position corresponds to different columns in the 8×12 plate:
- **Position A**: columns 1-3 (typically used for experimental conditions)
- **Position B**: columns 4-6 (added in first split-pool)
- **Position C**: columns 7-9 (added in second split-pool)
- **Position D**: columns 10-12 (added in third split-pool)

## Well Specification Formats

Both APIs support flexible well specifications:

| Format | Example | Meaning |
|--------|---------|---------|
| Row letter | `'A'`, `'B'` | All valid columns in that row |
| Column number | `1`, `2`, `3` | All rows in that column |
| Specific well | `'A1'`, `'B3'` | Single well |

**Note**: Column numbers must match the barcode position:
- Position A: columns 1-3 only
- Position B: columns 4-6 only
- Position C: columns 7-9 only
- Position D: columns 10-12 only

## Design Rationale

**Why separate APIs?**
- **Notebooks** need DataFrame integration and flexible well specs
- **Pipelines** need file I/O, CLI interface, and robust error handling

**Why shared core?**
- Eliminates code duplication (~150 lines → ~50 lines per API)
- Single source of truth for barcode logic
- Easier testing and maintenance
- Consistent behavior across use cases

**Error handling differences:**
- `barcode_annotation.py`: Raises exceptions (notebooks should fail fast)
- `assign_conditions.py`: Prints warnings (pipelines should be robust)

## Testing

```bash
# Test core utilities
python3 -c "
import sys; sys.path.append('scripts/utils')
from barcode_core import BC_SLICES, POSITION_COLUMNS
print('Positions:', list(BC_SLICES.keys()))
print('Position A cols:', POSITION_COLUMNS['A'])
"

# Test notebook API
python3 << 'EOF'
import sys, pandas as pd
sys.path.append('scripts/utils')
from barcode_annotation import annotate_barcodes

bcs = pd.read_csv('barcodes/map.csv', index_col=0)
bcs.columns = bcs.columns.astype(int)
df = pd.read_csv('results/sample/compiled_qc_metrics.csv', index_col=0)

result = annotate_barcodes(df, bcs, 'condition', {'PTA': ['A', 'B']})
print(result['condition'].value_counts())
EOF

# Test CLI
python3 scripts/CapWGS/assign_conditions.py \
    test_conditions.json barcodes/map.csv \
    results/sample/real_cells.txt -o /tmp/test.csv
```

## Migration Notes

If you have existing code using the old `barcode_annotation.py`, it should work without changes. If you were using the old `assign_conditions.py` from `scripts/CapWGS`, the only change is:

```python
# Old: only accepted specific wells
conditions = {"Control": ["A1", "A2", "B1", "B2"]}

# New: also accepts rows and columns
conditions = {"Control": ["A", "B"]}  # Simpler!
```

Both formats still work, so existing configs are compatible.
