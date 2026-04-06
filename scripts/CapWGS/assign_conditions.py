#!/usr/bin/env python3
"""
Assign experimental conditions to cell barcodes based on a conditions config
and the barcode plate map.

Conditions are encoded on barcode position A (columns 1-3 of the plate map,
the last 8bp of the 45-base combinatorial barcode). The other positions
(B/C/D) are added through split-pooling and are not condition-specific.

The conditions config is a simple JSON file mapping condition names to wells:
{
    "Control": ["A1", "A2", "B1", "B2"],
    "Drug_X": ["C1", "C2", "D1", "D2"]
}

Wells are specified as row letter + column number (1-3 only, since position A
uses columns 1-3 of the 8x12 plate). Row letters A-H, columns 1-3.

Usage:
    python3 assign_conditions.py <conditions.json> <barcode_map.csv> <real_cells.txt> \
        --output conditions_assignment.csv

    # Or use from Python:
    from assign_conditions import assign_conditions_to_barcodes
    assignments = assign_conditions_to_barcodes(config, barcode_map_path, cell_barcodes)
"""

import argparse
import csv
import json
import sys


# Barcode position A within the 45-base combinatorial barcode string:
# {bcD}AGGA{bcC}ACTC{bcB}AAGG{bcA}A
# 0-7  8-11 12-19 20-23 24-31 32-35 36-43 44
BC_A_SLICE = slice(36, 44)

# Position A uses plate columns 1-3
VALID_COLS = [1, 2, 3]


def load_barcode_map(map_csv):
    """Load the 8x12 barcode plate map.

    Returns dict: (row_letter, col_int) -> barcode_sequence
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


def build_barcode_to_condition(config, plate):
    """Build a lookup from barcode-A subsequence -> condition name.

    Args:
        config: dict mapping condition names to lists of well strings (e.g. "A1")
        plate: dict from load_barcode_map

    Returns:
        dict: barcode_subsequence -> condition_name
    """
    lookup = {}
    for condition, wells in config.items():
        for w in wells:
            if not isinstance(w, str) or len(w) < 2:
                print("Warning: unrecognized well spec: {}".format(w), file=sys.stderr)
                continue

            row_letter = w[0].upper()
            try:
                col = int(w[1:])
            except ValueError:
                print("Warning: can't parse column from well '{}'"
                      .format(w), file=sys.stderr)
                continue

            if col not in VALID_COLS:
                print("Warning: well {} column {} not valid for barcode position A "
                      "(valid: {})".format(w, col, VALID_COLS), file=sys.stderr)
                continue

            key = (row_letter, col)
            if key not in plate:
                print("Warning: well {} not found in plate map".format(w),
                      file=sys.stderr)
                continue

            seq = plate[key]
            if seq in lookup:
                print("Warning: barcode {} (well {}) already assigned to '{}', "
                      "overwriting with '{}'".format(seq, w, lookup[seq], condition),
                      file=sys.stderr)
            lookup[seq] = condition

    return lookup


def assign_conditions_to_barcodes(config, barcode_map_path, cell_barcodes):
    """Assign conditions to a list of cell barcodes.

    Args:
        config: dict mapping condition names to lists of well strings
        barcode_map_path: path to barcode map CSV (barcodes/map.csv)
        cell_barcodes: list of 45-base cell barcode strings

    Returns:
        dict: cell_barcode -> condition_name (or 'Unassigned')
    """
    plate = load_barcode_map(barcode_map_path)
    lookup = build_barcode_to_condition(config, plate)

    assignments = {}
    for barcode in cell_barcodes:
        subseq = barcode[BC_A_SLICE]
        assignments[barcode] = lookup.get(subseq, 'Unassigned')

    return assignments


def main():
    parser = argparse.ArgumentParser(
        description='Assign experimental conditions to cell barcodes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('conditions_json',
                        help='Conditions config JSON (maps condition names to wells)')
    parser.add_argument('barcode_map',
                        help='Barcode plate map CSV (barcodes/map.csv)')
    parser.add_argument('real_cells',
                        help='real_cells.txt (one barcode per line)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output CSV file (barcode, condition)')
    args = parser.parse_args()

    with open(args.conditions_json) as f:
        config = json.load(f)

    with open(args.real_cells) as f:
        barcodes = [line.strip() for line in f if line.strip()]

    assignments = assign_conditions_to_barcodes(config, args.barcode_map, barcodes)

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['barcode', 'condition'])
        for bc in barcodes:
            writer.writerow([bc, assignments[bc]])

    # Summary
    from collections import Counter
    counts = Counter(assignments.values())
    print("Condition assignments:")
    for cond, count in sorted(counts.items()):
        print("  {}: {} cells".format(cond, count))
    print("Written to: {}".format(args.output))


if __name__ == '__main__':
    main()
