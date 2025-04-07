import os
import pysam
from collections import defaultdict
import sys


def filter_and_split_sam_by_read_group(input_sam, output_dir, threshold):
    """
    Filter a SAM file by read group counts and split into multiple SAM files.

    Args:
        input_sam (str): Path to the input SAM file.
        output_dir (str): Directory to store the filtered and split SAM files.
        threshold (int): Minimum number of reads required for a read group to be kept.
    """
    # Open the input SAM file
    with pysam.AlignmentFile(input_sam, "r") as infile:
        header = infile.header.to_dict()

        # Count the occurrences of each read group
        read_group_counts = defaultdict(int)
        for read in infile:
            if read.has_tag("RG"):
                rg = read.get_tag("RG")
                read_group_counts[rg] += 1

    # Identify read groups that meet the threshold
    valid_read_groups = {rg for rg, count in read_group_counts.items() if count > threshold}

    # Reopen the input SAM to filter reads and split by read group
    with pysam.AlignmentFile(input_sam, "r") as infile:
        os.makedirs(output_dir, exist_ok=True)

        # Dictionary to store file handles for each valid read group
        rg_files = {}

        for read in infile:
            if read.has_tag("RG"):
                rg = read.get_tag("RG")
                if rg in valid_read_groups:
                    # Open a new SAM file for this RG if not already open
                    if rg not in rg_files:
                        output_file = os.path.join(output_dir, f"{rg}.sam")
                        rg_files[rg] = pysam.AlignmentFile(output_file, "w", header=header)

                    # Write the read to the appropriate SAM file
                    rg_files[rg].write(read)

        # Close all open SAM files
        for file in rg_files.values():
            file.close()

    # Write summary of filtered read groups and counts
    counts_file = os.path.join(output_dir, "filtered_read_group_counts.txt")
    with open(counts_file, "w") as f:
        f.write("ReadGroup\tReadCount\n")
        for rg, count in read_group_counts.items():
            if rg in valid_read_groups:
                f.write(f"{rg}\t{count}\n")

    print(f"Filtering and splitting complete. Output written to {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_sam> <output_dir> <threshold>")
        sys.exit(1)

    # Get input SAM file, output directory, and threshold from command-line arguments
    input_sam = sys.argv[1]
    output_dir = sys.argv[2]
    threshold = int(sys.argv[3])

    # Run the filter and split function
    filter_and_split_sam_by_read_group(input_sam, output_dir, threshold)
