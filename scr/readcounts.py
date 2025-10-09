#!/usr/bin/env python

import sys
import pandas as pd
import argparse
from tqdm import tqdm

def count_barcodes(output_file=None):
    # Initialize dictionary to hold barcode counts
    barcode_counts = {}

    # Process reads with tqdm progress bar, without pre-counting lines
    with tqdm(desc="Computing reads per barcode", unit="reads") as pbar:
        for line in sys.stdin:
            # Split the line by tabs to parse the fields
            fields = line.strip().split()

            # Find the barcode (BC:Z:<barcode>)
            barcode = None
            for field in fields[11:]:  # SAM optional fields start at the 12th column
                if field.startswith("CB:Z:"):
                    barcode = field.split(":", 2)[2]  # Extract the barcode part
                    break

            if barcode:
                # Increment count for this barcode
                if barcode in barcode_counts:
                    barcode_counts[barcode] += 1
                else:
                    barcode_counts[barcode] = 1

            # Update progress bar after each read
            pbar.update(1)

    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(list(barcode_counts.items()), columns=["barcode", "read_count"])

    # Output to CSV or standard output
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"Saved barcode read counts to {output_file}")
    else:
        print(df.to_csv(index=False))

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Count reads per barcode in a BAM file.")
    parser.add_argument("-o", "--output", help="Name of the output CSV file to save the barcode read counts. If not specified, prints to stdout.")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the count_barcodes function with the output file argument
    count_barcodes(args.output)

if __name__ == "__main__":
    main()

