import gzip
from collections import defaultdict

# Function to parse the barcode and assign read groups
def assign_read_groups(fastq_file, output_fastq, count_output_file):
    barcode_to_group = {}
    group_counter = 1
    cell_counts = defaultdict(int)  # Dictionary to track counts for each cell

    # Open the input FASTQ file
    with gzip.open(fastq_file, 'rt') as infile, gzip.open(output_fastq, 'wt') as outfile:
        while True:
            # Read four lines for each FASTQ record
            header = infile.readline().strip()
            if not header:
                break
            seq = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Extract the barcode (CB:Z:) from the header
            if 'CB:Z:' in header:
                barcode = header.split('CB:Z:')[1].split()[0]  # Extract barcode

                # Assign a new read group if barcode is new
                if barcode not in barcode_to_group:
                    barcode_to_group[barcode] = f"cell{group_counter}"
                    group_counter += 1

                # Increment the count for the cell
                read_group = barcode_to_group[barcode]
                cell_counts[read_group] += 1

                # Remove CB:Z: from the header and add RG:Z:
                new_header = header.split('CB:Z:')[0].strip() + f" RG:Z:{read_group}"
            else:
                # Skip reads without a barcode
                new_header = header

            # Write the updated record to the output FASTQ file
            outfile.write(f"{new_header}\n{seq}\n{plus}\n{quality}\n")

    # Save cell counts to a file
    with open(count_output_file, 'w') as count_file:
        count_file.write("Cell\tReadCount\n")
        for cell, count in cell_counts.items():
            count_file.write(f"{cell}\t{count}\n")

    print(f"Finished processing. Read groups assigned: {len(barcode_to_group)}")
    print(f"Cell counts saved to {count_output_file}")

# Input and output file paths
input_fastq = "./read1_10_90_K562_PBMC.fastq.gz"  # Replace with your input FASTQ file path
output_fastq = "read1_10_90_K562_PBMC_with_read_groups.fastq.gz"  # Replace with your desired output path
count_output_file = "read1_read1_10_90_K562_PBMC_cell_counts.txt"  # Replace with your desired cell count file path

# Run the function
assign_read_groups(input_fastq, output_fastq, count_output_file)

# Input and output file paths
input_fastq = "./read2_10_90_K562_PBMC.fastq.gz"  # Replace with your input FASTQ file path
output_fastq = "read2_10_90_K562_PBMC_with_read_groups.fastq.gz"  # Replace with your desired output path
count_output_file = "read2_read1_10_90_K562_PBMC_cell_counts.txt"  # Replace with your desired cell count file path

# Run the function
assign_read_groups(input_fastq, output_fastq, count_output_file)