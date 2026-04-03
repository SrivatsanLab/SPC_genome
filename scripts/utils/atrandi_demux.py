#!/usr/bin/env python

import gzip
import json
import numpy as np
import os
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

def load_correction_dict(dict_path=None):
    """Load BFS-based barcode correction dictionary from JSON.
    Keys A/B/C/D map observed 8bp sequences to their corrected true sequences.
    """
    if dict_path is None:
        dict_path = os.path.join(os.path.dirname(__file__), '..', '..', 'barcodes', 'barcode_correction_dicts.json')
    with open(dict_path) as f:
        return json.load(f)

def correction(sequence, true_overhangs, correction_dict, ignore_overhangs=False):
    """
    sequence         : The DNA sequence from read2 (string)
    true_overhangs   : The overhang sequences (dictionary)
    correction_dict  : BFS correction dictionary from load_correction_dict()
    ignore_overhangs : If True, skip overhang sequence validation (boolean, default False)
    """
    # Parse index sequence
    seq_overhangs = {"D": sequence[8:12],
                     "C": sequence[20:24],
                     "B": sequence[32:36],
                     "A": sequence[44:45]}

    seq_bcs = {"D": sequence[0:8],
               "C": sequence[12:20],
               "B": sequence[24:32],
               "A": sequence[36:44]}

    corrected_index = ""

    for l in ["D","C","B","A"]:

        # overhang sequences:
        if not ignore_overhangs:
            overhang_dist = distance(seq_overhangs[l], true_overhangs[l])
            if overhang_dist > 1:
                return None

        # barcode lookup
        corrected_bc = correction_dict[l].get(seq_bcs[l])
        if corrected_bc is None:
            return None
        corrected_index += corrected_bc + true_overhangs[l]

    # Return corrected index, genomic sequence
    genomic_seq = sequence[46:]
    return corrected_index, genomic_seq
        

def correct_fastqs_gz(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, ignore_overhangs=False):
    """
    Read in R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file          : input read1 fastq, gzipped
    R2_file          : input read2 fastq, gzipped
    R1_output        : path for read1 output, should end with .fastq.gz
    R2_output        : path for read2 output, should end with .fastq.gz
    correction_dict  : BFS correction dictionary from load_correction_dict()
    true_overhangs   : dictionary of true A,B,C,D overhangs
    ignore_overhangs : if True, skip overhang validation (boolean, default False)

    """

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

                corrected = correction(read2_seq, true_overhangs, correction_dict, ignore_overhangs)
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
        
def correct_fastqs(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, ignore_overhangs=False):
    """
    Read in UNZIPPED R1 and R2 files, correct barcodes, and write corrected sequences.
    R1_file          : input read1 fastq
    R2_file          : input read2 fastq
    R1_output        : path for read1 output
    R2_output        : path for read2 output
    correction_dict  : BFS correction dictionary from load_correction_dict()
    true_overhangs   : dictionary of true A,B,C,D overhangs
    ignore_overhangs : if True, skip overhang validation (boolean, default False)

    """

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

                corrected = correction(read2_seq, true_overhangs, correction_dict, ignore_overhangs)
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

def rejects(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, ignore_overhangs=False):
    """
    Read in R1 and R2 files, write reads whose barcodes don't map.
    R1_file          : input read1 fastq
    R2_file          : input read2 fastq
    correction_dict  : BFS correction dictionary from load_correction_dict()
    true_overhangs   : dictionary of true A,B,C,D overhangs
    ignore_overhangs : if True, skip overhang validation (boolean, default False)

    """

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

                corrected = correction(read2_seq, true_overhangs, correction_dict, ignore_overhangs)
                if corrected == None:
                    R1out.writelines(read1_lines)
                    R2out.writelines(read2_lines)
                    
                else:
                    pass
                
                pbar.update(1)  # Update progress bar

def rejects_gz(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, ignore_overhangs=False):
    """
    Read in R1 and R2 files, write reads whose barcodes don't map (gzipped).
    R1_file          : input read1 fastq, gzipped
    R2_file          : input read2 fastq, gzipped
    correction_dict  : BFS correction dictionary from load_correction_dict()
    true_overhangs   : dictionary of true A,B,C,D overhangs
    ignore_overhangs : if True, skip overhang validation (boolean, default False)

    """

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

                corrected = correction(read2_seq, true_overhangs, correction_dict, ignore_overhangs)
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
    
    # Optional arguments with default values
    parser.add_argument("--R1_output", type=str, default="R1_corrected.fastq.gz", 
                        help="Output file name for corrected R1 (default: R1_corrected.fastq.gz).")
    parser.add_argument("--R2_output", type=str, default="R2_corrected.fastq.gz", 
                        help="Output file name for corrected R2 (default: R2_corrected.fastq.gz).")
    parser.add_argument("--gzip", type=str, default="True",
                        help="Are the input fastqs gzipped (True/False, default=True)")
    parser.add_argument("--rejects", action='store_true',
                        help="add this flag to run in 'rejects mode' which writes reads whose barcodes don't map to separate files: R[12]_unmapped.fastq.gz")
    parser.add_argument("--ignore-overhangs", action='store_true',
                        help="Skip overhang sequence validation (use when constant sequences have systematic errors)")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Extract arguments
    R1_file = args.R1_file
    R2_file = args.R2_file
    R1_output = args.R1_output
    R2_output = args.R2_output
    
    if args.gzip == "True":
        gzip = True
    elif args.gzip == "False":
        gzip = False

    # Load BFS correction dictionary
    correction_dict = load_correction_dict()

    # Overhangs defined in atrandi library structure
    true_overhangs = {'A':'A',
                      'B':'AAGG',
                      'C':'ACTC',
                      'D':'AGGA'
                      }
    # if running in rejects mode, only print rejects
    if args.rejects:
        if gzip:
            rejects_gz(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, args.ignore_overhangs)
        else:
            rejects(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, args.ignore_overhangs)
    else:
        if gzip:
            correct_fastqs_gz(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, args.ignore_overhangs)

        else:
            correct_fastqs(R1_file, R2_file, R1_output, R2_output, correction_dict, true_overhangs, args.ignore_overhangs)

if __name__ == "__main__":
    main()