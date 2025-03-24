#!/usr/bin/env python

import gzip
import numpy as np
import sys
from Levenshtein import distance
from tqdm import tqdm
import argparse

def get_bcs(fileA,fileB,fileC,fileD):
    true_BCs={}
    true_BCs['A'] = []
    true_BCs['B'] = []
    true_BCs['C'] = []
    true_BCs['D'] = []
    with open(fileA, 'r') as fin:
        for line in fin:
            line = line.strip()
            true_BCs['A'].append(line)
    with open(fileB, 'r') as fin:
        for line in fin:
            line = line.strip()
            true_BCs['B'].append(line)
    with open(fileC, 'r') as fin:
        for line in fin:
            line = line.strip()
            true_BCs['C'].append(line)
    with open(fileD, 'r') as fin:
        for line in fin:
            line = line.strip()
            true_BCs['D'].append(line)
            
    return true_BCs

def pairwise_levenshtein_distances(string, vector_of_strings):
    """
    Calculate pairwise levenshtein distances between a string and a vector of strings.
    """
    return [distance(string, seq) for seq in vector_of_strings]
    
def find_max_distance(true_BCs):
    BC_set_distances = {}
    for key, strings in true_BCs.items():
        distances = []
        for i, string1 in enumerate(strings):
            for j, string2 in enumerate(strings):
                if i < j:  # Avoid duplicate pairs
                    dist = distance(string1, string2)
                    distances.append(dist)
        BC_set_distances[key] = distances

    max_distances = {}
    for key, distances in BC_set_distances.items():
        max_distances[key] =  min(distances)-2
    return max_distances

def correction(sequence, true_overhangs, true_BCs, max_distances):
    """
    sequence       : The DNA sequence from read2 (string)
    true_overhangs : The overhang sequences (dictionary)
    true_BCs       : The barcodes (dictionary from get_bcs())
    max_distances  : The maximum permissible levenshtein distances (dictionary from find_max_distance())
    """
    # Parse index sequence
    parts = [sequence[0:8], sequence[8:12], sequence[12:20], sequence[20:24], sequence[24:32], sequence[32:36], sequence[36:44], sequence[44:45]]
    seq_overhangs = {"D":parts[1],
                     "C":parts[3],
                     "B":parts[5],
                     "A":parts[7]}

    seq_bcs = {"D":parts[0],
               "C":parts[2],
               "B":parts[4],
               "A":parts[6]}

    # compute levenshtein distances, perform error correction
    corrected_index = ""
    
    bc_distances = {}
    for l in ["D","C","B","A"]:
        
        # overhang sequences:
        overhang_dist = distance(seq_overhangs[l], true_overhangs[l])
        if overhang_dist > 1:
            return None
        
        # barcode distances
        distances = pairwise_levenshtein_distances(seq_bcs[l], true_BCs[l])
        bc_dist = min(distances)
        bc_distances[l] = bc_dist
        corected_bc = true_BCs[l][np.argmin(distances)]
        corrected_index += corected_bc+true_overhangs[l]
    
    # Return corrected index, genomic sequence
    genomic_seq = sequence[46:]
    
    if bc_distances['A'] <= max_distances['A'] and bc_distances['B'] <= max_distances['B'] and bc_distances['C'] <= max_distances['C'] and bc_distances['D'] <= max_distances['D']:
        return corrected_index,genomic_seq
    else:
        return None
        

def correct_fastqs_gz(R1_file, R2_file, R1_output, R2_output, true_BCs, true_overhangs):
    """
    Read in R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file        : input read1 fastq, gzipped
    R2_file        : input read2 fastq, gzipped
    R1_output      : path for read1 output, should end with .fastq.gz
    R2_output      : path for read2 output, should end with .fastq.gz
    true_BCs       : dictionary of true A,B,C,D barcodes from get_bcs()
    true_overhangs : dictionary of true A,B,C,D overhangs
    
    """
    
    max_distances = find_max_distance(true_BCs)
    
    with gzip.open(R1_file, "rt") as R1in, gzip.open(R1_output, "wt") as R1out, gzip.open(R2_file, "rt") as R2in, gzip.open(R2_output, "wt") as R2out:
        
        mapped_reads = 0
        total_reads = 0
        
        total_reads_estimated = total_reads or sum(1 for _ in R1in) // 4
        R1in.seek(0)  # Reset file pointer after estimating total reads
        
        with tqdm(total=total_reads_estimated, desc="Processing reads", unit="read pair") as pbar:
            
            while True:
                # Read 4 lines from each file (1 FASTQ record)
                read1_lines = [R1in.readline() for _ in range(4)]
                read2_lines = [R2in.readline() for _ in range(4)]

                # Break if we've reached the end of either file
                if not read1_lines[0] or not read2_lines[0]:
                    break

                total_reads += 1

                read2_seq = read2_lines[1].strip()

                corrected = correction(read2_seq, true_overhangs, true_BCs, max_distances) # if barcode maps, corrected = (corrected_index,genomic_seq). else it will = None
                if corrected == None:
                    # print("barcode not mapped")
                    pass
                else:
                    mapped_reads += 1
                    # add barcode to read1 header
                    read1_lines[0] = read1_lines[0].strip().split()
                    # check for illumina demultiplexing index
                    if len(read1_lines[0][1]) > 7:
                        read1_lines[0] = f"{read1_lines[0][0]} BC:Z:{read1_lines[0][1][7:]} CB:Z:{corrected[0]}\n"  # Add barcode to header as optional flag, also preserve illumina multiplexing as optional flag
                    elif len(read1_lines[0][1]) == 7:
                        read1_lines[0] = f"{read1_lines[0][0]} CB:Z:{corrected[0]}\n"                               # Add barcode to header as optional flag, no illumina multiplexing index
                    R1out.writelines(read1_lines)

                    # add barcode to read2 header, delete index from read2 sequence (and shorten quality scores to match)
                    read2_lines[0] = read2_lines[0].strip().split()
                    # check for illumina demultiplexing index
                    if len(read2_lines[0][1]) > 7:
                        read2_lines[0] = f"{read2_lines[0][0]} BC:Z:{read2_lines[0][1][7:]} CB:Z:{corrected[0]}\n"  # Add barcode to header as optional flag, also preserve illumina multiplexing as optional flag
                    elif len(read2_lines[0][1]) == 7:
                        read2_lines[0] = f"{read2_lines[0][0]} CB:Z:{corrected[0]}\n"                               # Add barcode to header as optional flag, no illumina multiplexing index
                    read2_lines[1] = f"{corrected[1]}\n"                                                            # Remove barcode from genomic seq
                    read2_lines[3] = f"{read2_lines[3][len(corrected[0])+1:]}"                                      # Remove barcode from quality scores
                    R2out.writelines(read2_lines)
                
                pbar.update(1)  # Update progress bar
        
        print(f"proportion of reads with barcode mapped: {mapped_reads/total_reads}")
        
def correct_fastqs(R1_file, R2_file, R1_output, R2_output, true_BCs, true_overhangs):
    """
    Read in UNZIPPED R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file        : input read1 fastq
    R2_file        : input read2 fastq
    R1_output      : path for read1 output
    R2_output      : path for read2 output
    true_BCs       : dictionary of true A,B,C,D barcodes from get_bcs()
    true_overhangs : dictionary of true A,B,C,D overhangs
    
    """
    
    max_distances = find_max_distance(true_BCs)
    
    with open(R1_file, "r") as R1in, open(R1_output, "w") as R1out, open(R2_file, "r") as R2in, open(R2_output, "w") as R2out:
        
        mapped_reads = 0
        total_reads = 0
        
        total_reads_estimated = total_reads or sum(1 for _ in R1in) // 4
        R1in.seek(0)  # Reset file pointer after estimating total reads
        
        with tqdm(total=total_reads_estimated, desc="Processing reads", unit="read pair") as pbar:
            while True:
                # Read 4 lines from each file (1 FASTQ record)
                read1_lines = [R1in.readline() for _ in range(4)]
                read2_lines = [R2in.readline() for _ in range(4)]

                # Break if we've reached the end of either file
                if not read1_lines[0] or not read2_lines[0]:
                    break

                total_reads += 1

                read2_seq = read2_lines[1].strip()

                corrected = correction(read2_seq, true_overhangs, true_BCs, max_distances) # if barcode maps, corrected = (corrected_index,genomic_seq). else it will = None
                if corrected == None:
                    # print("barcode not mapped")
                    pass
                else:
                    mapped_reads += 1
                    # add barcode to read1 header
                    read1_lines[0] = read1_lines[0].strip().split()
                    # check for illumina demultiplexing index
                    if len(read1_lines[0][1]) > 7:
                        read1_lines[0] = f"{read1_lines[0][0]} BC:Z:{read1_lines[0][1][7:]} CB:Z:{corrected[0]}\n"  # Add barcode to header as optional flag, also preserve illumina multiplexing as optional flag
                    elif len(read1_lines[0][1]) == 7:
                        read1_lines[0] = f"{read1_lines[0][0]} CB:Z:{corrected[0]}\n"                               # Add barcode to header as optional flag, no illumina multiplexing index
                    R1out.writelines(read1_lines)

                    # add barcode to read2 header, delete index from read2 sequence (and shorten quality scores to match)
                    read2_lines[0] = read2_lines[0].strip().split()
                    # check for illumina demultiplexing index
                    if len(read2_lines[0][1]) > 7:
                        read2_lines[0] = f"{read2_lines[0][0]} BC:Z:{read2_lines[0][1][7:]} CB:Z:{corrected[0]}\n"  # Add barcode to header as optional flag, also preserve illumina multiplexing as optional flag
                    elif len(read2_lines[0][1]) == 7:
                        read2_lines[0] = f"{read2_lines[0][0]} CB:Z:{corrected[0]}\n"                               # Add barcode to header as optional flag, no illumina multiplexing index
                    read2_lines[1] = f"{corrected[1]}\n"                                                            # Remove barcode from genomic seq
                    read2_lines[3] = f"{read2_lines[3][len(corrected[0])+1:]}"                                      # Remove barcode from quality scores
                    R2out.writelines(read2_lines)
                
                pbar.update(1)  # Update progress bar
        
        print(f"proportion of reads with barcode mapped: {mapped_reads/total_reads}")

def rejects(R1_file, R2_file, R1_output, R2_output, true_BCs, true_overhangs):
    """
    Read in R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file        : input read1 fastq, gzipped
    R2_file        : input read2 fastq, gzipped
    true_BCs       : dictionary of true A,B,C,D barcodes from get_bcs()
    true_overhangs : dictionary of true A,B,C,D overhangs
    
    """
    
    max_distances = find_max_distance(true_BCs)
    with open(R1_file, "r") as R1in, open(R1_output, "w") as R1out, open(R2_file, "r") as R2in, open(R2_output, "w") as R2out:
        
        mapped_reads = 0
        total_reads = 0
        
        total_reads_estimated = total_reads or sum(1 for _ in R1in) // 4
        R1in.seek(0)  # Reset file pointer after estimating total reads
        
        with tqdm(total=total_reads_estimated, desc="Processing reads", unit="read pair") as pbar:
            
            while True:
                # Read 4 lines from each file (1 FASTQ record)
                read1_lines = [R1in.readline() for _ in range(4)]
                read2_lines = [R2in.readline() for _ in range(4)]

                # Break if we've reached the end of either file
                if not read1_lines[0] or not read2_lines[0]:
                    break

                total_reads += 1

                read2_seq = read2_lines[1].strip()

                corrected = correction(read2_seq, true_overhangs, true_BCs, max_distances) # if barcode maps, corrected = (corrected_index,genomic_seq). else it will = None
                if corrected == None:
                    R1out.writelines(read1_lines)
                    R2out.writelines(read2_lines)
                    
                else:
                    pass
                
                pbar.update(1)  # Update progress bar

def rejects_gz(R1_file, R2_file, R1_output, R2_output, true_BCs, true_overhangs):
    """
    Read in R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file        : input read1 fastq, gzipped
    R2_file        : input read2 fastq, gzipped
    true_BCs       : dictionary of true A,B,C,D barcodes from get_bcs()
    true_overhangs : dictionary of true A,B,C,D overhangs
    
    """
    
    max_distances = find_max_distance(true_BCs)
    
    with gzip.open(R1_file, "rt") as R1in, gzip.open(R1_output, "wt") as R1out, gzip.open(R2_file, "rt") as R2in, gzip.open(R2_output, "wt") as R2out:
        
        mapped_reads = 0
        total_reads = 0
        
        total_reads_estimated = total_reads or sum(1 for _ in R1in) // 4
        R1in.seek(0)  # Reset file pointer after estimating total reads
        
        with tqdm(total=total_reads_estimated, desc="Processing reads", unit="read pair") as pbar:
            
            while True:
                # Read 4 lines from each file (1 FASTQ record)
                read1_lines = [R1in.readline() for _ in range(4)]
                read2_lines = [R2in.readline() for _ in range(4)]

                # Break if we've reached the end of either file
                if not read1_lines[0] or not read2_lines[0]:
                    break

                total_reads += 1

                read2_seq = read2_lines[1].strip()

                corrected = correction(read2_seq, true_overhangs, true_BCs, max_distances) # if barcode maps, corrected = (corrected_index,genomic_seq). else it will = None
                if corrected == None:
                    R1out.writelines(read1_lines)
                    R2out.writelines(read2_lines)
                    
                else:
                    pass
                
                pbar.update(1)  # Update progress bar


def main():

    parser = argparse.ArgumentParser(description="Index mapping and correction + fastq modification for paired end fastqs from Atrandi style combinatorial indexing")
    
    # Positional arguments
    parser.add_argument("R1_file", type=str, help="Path to R1 FASTQ file (gzipped)")
    parser.add_argument("R2_file", type=str, help="Path to R2 FASTQ file (gzipped)")
    parser.add_argument("barcode_dir", type=str, help="Path to directory containing barcode files.")
    
    # Optional arguments with default values
    parser.add_argument("--R1_output", type=str, default="R1_corrected.fastq.gz", 
                        help="Output file name for corrected R1 (default: R1_corrected.fastq.gz).")
    parser.add_argument("--R2_output", type=str, default="R2_corrected.fastq.gz", 
                        help="Output file name for corrected R2 (default: R2_corrected.fastq.gz).")
    parser.add_argument("--gzip", type=str, default="True",
                        help="Are the input fastqs gzipped (True/False, default=True)")
    parser.add_argument("--rejects", action='store_true',
                        help="add this flag to run in 'rejects mode' which writes reads whose barcodes don't map to separate files: R[12]_unmapped.fastq.gz")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Extract arguments
    R1_file = args.R1_file
    R2_file = args.R2_file
    barcode_dir = args.barcode_dir
    R1_output = args.R1_output
    R2_output = args.R2_output
    
    if args.gzip == "True":
        gzip = True
    elif args.gzip == "False":
        gzip = False

    # Get barcodes (supplied in 4 files in barcode_dir)
    true_barcodes = get_bcs(f"{barcode_dir}/bcA_24_exD.txt", f"{barcode_dir}/bcB_24_exC.txt", f"{barcode_dir}/bcC_24_exB.txt", f"{barcode_dir}/bcD_24_exA.txt")
    
    # Overhangs defined in atrandi library structure
    true_overhangs = {'A':'A',
                      'B':'AAGG',
                      'C':'ACTC',
                      'D':'AGGA'
                      }
    # if running in rejects mode, only print rejects
    if args.rejects:
        if gzip:
            rejects_gz(R1_file, R2_file, R1_output, R2_output, true_barcodes, true_overhangs)
        else: 
            rejects_gz(R1_file, R2_file, R1_output, R2_output, true_barcodes, true_overhangs)
    else:
        if gzip:
            correct_fastqs_gz(R1_file, R2_file, R1_output, R2_output, true_barcodes, true_overhangs)

        else:
            correct_fastqs(R1_file, R2_file, R1_output, R2_output, true_barcodes, true_overhangs)

if __name__ == "__main__":
    main()