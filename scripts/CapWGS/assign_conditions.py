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

This is a CLI tool designed for pipeline automation and dashboard generation.
For interactive DataFrame annotation in notebooks, see: scripts/utils/barcode_annotation.py

Usage:
    python3 assign_conditions.py <conditions.json> <barcode_map.csv> <real_cells.txt> \\
        --output conditions_assignment.csv

    # Or use from Python:
    from assign_conditions import assign_conditions_to_barcodes
    assignments = assign_conditions_to_barcodes(config, barcode_map_path, cell_barcodes)
"""

import argparse
import csv
import json
import os
import sys

# Add utils to path for barcode_core import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'utils'))

from barcode_core import (
    load_barcode_map_dict, build_barcode_lookup, assign_conditions
)


def assign_conditions_to_barcodes(config, barcode_map_path, cell_barcodes, position='A'):
    """Assign conditions to a list of cell barcodes.

    Args:
        config: dict mapping condition names to lists of well strings
        barcode_map_path: path to barcode map CSV (barcodes/map.csv)
        cell_barcodes: list of 45-base cell barcode strings
        position: barcode position to use (default: 'A')

    Returns:
        dict: cell_barcode -> condition_name (or 'Unassigned')
    """
    plate = load_barcode_map_dict(barcode_map_path)
    lookup = build_barcode_lookup(config, plate, position, warn=True)
    assignments = assign_conditions(cell_barcodes, lookup, position)
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
    parser.add_argument('--position', default='A',
                        choices=['A', 'B', 'C', 'D'],
                        help='Barcode position to use for conditions (default: A)')
    args = parser.parse_args()

    with open(args.conditions_json) as f:
        config = json.load(f)

    with open(args.real_cells) as f:
        barcodes = [line.strip() for line in f if line.strip()]

    assignments = assign_conditions_to_barcodes(
        config, args.barcode_map, barcodes, args.position
    )

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
