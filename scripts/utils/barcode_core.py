#!/usr/bin/env python3
"""
Core utilities for CapWGS combinatorial barcode annotation.

Cell barcodes have the format (45 bases):
    {bcD}AGGA{bcC}ACTC{bcB}AAGG{bcA}A
    0-7  8-11 12-19 20-23 24-31 32-35 36-43 44

Each position (A/B/C/D) corresponds to a different set of columns in the
8x12 barcode plate:
    - Position A: columns 1-3
    - Position B: columns 4-6
    - Position C: columns 7-9
    - Position D: columns 10-12

This module provides shared utilities for both:
    - Interactive DataFrame annotation (barcode_annotation.py)
    - Pipeline/CLI tools (assign_conditions.py)
"""

import csv
import sys
from pathlib import Path

# Barcode position slices within the 45-base combinatorial string
BC_SLICES = {
    'D': slice(0, 8),
    'C': slice(12, 20),
    'B': slice(24, 32),
    'A': slice(36, 44),
}

# Column ranges for each barcode position in the 8x12 plate
POSITION_COLUMNS = {
    'A': [1, 2, 3],
    'B': [4, 5, 6],
    'C': [7, 8, 9],
    'D': [10, 11, 12],
}

# Valid row letters
VALID_ROWS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']


def load_barcode_map_dict(map_csv):
    """Load the 8x12 barcode plate map into a dict.

    Args:
        map_csv: Path to barcode map CSV file

    Returns:
        dict: (row_letter, col_int) -> barcode_sequence
    """
    plate = {}
    with open(map_csv) as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = [int(c) for c in header[1:]]
        for row in reader:
            row_letter = row[0]
            for i, col in enumerate(cols):
                plate[(row_letter, col)] = row[i + 1]
    return plate


def load_barcode_map_dataframe(map_csv):
    """Load the 8x12 barcode plate map into a pandas DataFrame.

    Args:
        map_csv: Path to barcode map CSV file

    Returns:
        pd.DataFrame: 8x12 DataFrame (rows A-H, columns 1-12)
    """
    import pandas as pd
    df = pd.read_csv(map_csv, index_col=0)
    df.columns = df.columns.astype(int)
    return df


def parse_well_spec(well_spec, position='A'):
    """Parse a well specification into row/column pair(s).

    Args:
        well_spec: Can be:
            - Row letter: 'A', 'B', etc. (all valid columns for position)
            - Column number: 1, 2, etc. (all rows)
            - Specific well: 'A1', 'B3', etc.
        position: Barcode position ('A', 'B', 'C', 'D') for column validation

    Returns:
        list of (row_letter, col_int) tuples

    Raises:
        ValueError: If well_spec is invalid or column not valid for position
    """
    valid_cols = POSITION_COLUMNS[position]

    # Integer: column number (all rows)
    if isinstance(well_spec, int):
        if well_spec not in valid_cols:
            raise ValueError(
                f"Column {well_spec} is not valid for barcode position {position} "
                f"(valid columns: {valid_cols})"
            )
        return [(row, well_spec) for row in VALID_ROWS]

    # String: row letter OR specific well
    if isinstance(well_spec, str):
        # Single letter: row (all valid columns)
        if len(well_spec) == 1:
            row = well_spec.upper()
            if row not in VALID_ROWS:
                raise ValueError(f"Invalid row letter: {row} (valid: {VALID_ROWS})")
            return [(row, col) for col in valid_cols]

        # Specific well: 'A1', 'B12', etc.
        row = well_spec[0].upper()
        try:
            col = int(well_spec[1:])
        except ValueError:
            raise ValueError(f"Cannot parse column from well '{well_spec}'")

        if row not in VALID_ROWS:
            raise ValueError(f"Invalid row in well '{well_spec}': {row}")
        if col not in valid_cols:
            raise ValueError(
                f"Well {well_spec} column {col} is not valid for barcode position {position} "
                f"(valid columns: {valid_cols})"
            )
        return [(row, col)]

    raise ValueError(f"Unrecognized well spec type: {well_spec} ({type(well_spec)})")


def build_barcode_lookup(conditions, plate, position='A', warn=False):
    """Build a lookup from barcode subsequence to condition name.

    Args:
        conditions: dict mapping condition names to lists of well specs
            Example: {'Control': ['A', 'B'], 'Treatment': ['C1', 'D2']}
        plate: dict from load_barcode_map_dict() OR pandas DataFrame
        position: Barcode position to use ('A', 'B', 'C', 'D')
        warn: If True, print warnings instead of raising exceptions

    Returns:
        dict: barcode_subsequence -> condition_name

    Raises:
        ValueError: If duplicate barcode assignments (unless warn=True)
    """
    # Handle both dict and DataFrame plate formats
    if hasattr(plate, 'loc'):  # DataFrame
        def get_barcode(row, col):
            return plate.loc[row, col]
    else:  # dict
        def get_barcode(row, col):
            return plate.get((row, col))

    lookup = {}

    for condition, wells in conditions.items():
        for well_spec in wells:
            try:
                well_pairs = parse_well_spec(well_spec, position)
            except ValueError as e:
                if warn:
                    print(f"Warning: {e}", file=sys.stderr)
                    continue
                else:
                    raise

            for row, col in well_pairs:
                seq = get_barcode(row, col)
                if seq is None:
                    msg = f"Well {row}{col} not found in plate map"
                    if warn:
                        print(f"Warning: {msg}", file=sys.stderr)
                        continue
                    else:
                        raise ValueError(msg)

                if seq in lookup:
                    msg = (
                        f"Barcode {seq} (well {row}{col}) assigned to both "
                        f"'{lookup[seq]}' and '{condition}'"
                    )
                    if warn:
                        print(f"Warning: {msg}, overwriting with '{condition}'",
                              file=sys.stderr)
                    else:
                        raise ValueError(msg)

                lookup[seq] = condition

    return lookup


def extract_barcode_position(cell_barcode, position='A'):
    """Extract barcode subsequence for a specific position.

    Args:
        cell_barcode: 45-base cell barcode string
        position: Position to extract ('A', 'B', 'C', 'D')

    Returns:
        8-base barcode subsequence
    """
    return cell_barcode[BC_SLICES[position]]


def assign_conditions(cell_barcodes, lookup, position='A', unassigned_label='Unassigned'):
    """Assign conditions to a list of cell barcodes.

    Args:
        cell_barcodes: List of 45-base cell barcode strings
        lookup: Dict from build_barcode_lookup()
        position: Barcode position used in lookup
        unassigned_label: Label for cells not in lookup

    Returns:
        dict: cell_barcode -> condition_name
    """
    assignments = {}
    for barcode in cell_barcodes:
        subseq = extract_barcode_position(barcode, position)
        assignments[barcode] = lookup.get(subseq, unassigned_label)
    return assignments
