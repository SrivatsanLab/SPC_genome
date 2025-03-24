#!/usr/bin/env python

import sys
import pandas as pd
import argparse
from tqdm import tqdm
import pysam
from collections import defaultdict

def compute(bam_file, output_file=None):
    # Dictionaries for barcode counts and unique reads
    barcode_counts = defaultdict(int)
    barcode_reads = defaultdict(list)

    # Dictionary to store total insert length per barcode
    insert_lengths = defaultdict(list)

    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Process reads with a progress bar
        with tqdm(desc="Processing BAM file", unit="reads") as pbar:
            for read in bam:
                if not read.is_read1:
                    # Process only read1 for insert length calculation
                    continue

                # Parse the barcode
                cell_barcode = None
                for tag, value in read.tags:
                    if tag == "BC":
                        # Handle concatenated BC and CB
                        if "CB:Z:" in value:
                            parts = value.split(" CB:Z:")
                            cell_barcode = parts[1]
                        else:
                            continue
                    elif tag == "CB":
                        cell_barcode = value
                    
                if not cell_barcode:
                    # Skip if no cell barcode is found
                    continue
                
                if read.is_unmapped or read.mate_is_unmapped:
                    # Skip reads where either mate is unmapped, but increment unmapped count
                    unmapped_counts[cell_barcode] += 1
                    continue

                # Increment total read count for this barcode
                barcode_counts[cell_barcode] += 1

                # Compute insert length (requires paired-end reads)
                if read.is_proper_pair:
                    # Calculate 5' and 3' coordinates for insert length
                    five_prime = read.reference_start + 1  # 5' end (1-based)
                    three_prime = read.next_reference_start + read.template_length

                    barcode_reads[cell_barcode].append((five_prime, three_prime))

                    insert_length = abs(three_prime - five_prime)

                    # Store insert length for the barcode, clip alignment artifacts (inserts longer thank 1kb)
                    if insert_length < 1000:
                        insert_lengths[cell_barcode].append(insert_length)

                # Update progress bar
                pbar.update(1)

    # Compute average insert length per barcode
    avg_insert_lengths = {
        barcode: sum(lengths) / len(lengths) if lengths else 0
        for barcode, lengths in insert_lengths.items()
    }

    unique_read_counts = {barcode:len(set(barcode_reads[barcode])) for barcode in barcode_reads.keys()}


    # Create DataFrames for outputs
    barcode_counts_df = pd.DataFrame(list(barcode_counts.items()), columns=["barcode", "read_count"])
    unmapped_counts_df = pd.DataFrame(list(unmapped_counts.items()), columns=["barcode", "unmapped_read_count"])

    insert_length_df = pd.DataFrame(
        list(avg_insert_lengths.items()), columns=["barcode", "average_insert_length"]
    )

    unique_read_counts_df = pd.DataFrame(list(unique_read_counts.items()), columns=["barcode", "unique_read_count"])

    # Merge both DataFrames on barcode
    result_df = pd.merge(pd.merge(pd.merge(barcode_counts_df, unique_read_counts_df, on="barcode", how="outer"), insert_length_df, on="barcode", how="outer"), unmapped_counts_df,on="barcode", how="outer").fillna(0)
    
    # Output to CSV or standard output
    if output_file:
        result_df.to_csv(output_file, index=False)
        print(f"Saved barcode read counts and insert lengths to {output_file}")
    else:
        print(result_df.to_csv(index=False))

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Count unique reads per barcode and gDNA insert length in a paired end BAM file.")
    parser.add_argument("-i", "--input", help="path to input bam")
    parser.add_argument("-o", "--output", help="desired output file path, if not provided prints to stdout")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the count_barcodes function with the output file argument
    compute(args.input, args.output)

if __name__ == "__main__":
    main()