Installation:

for now, simply clone the repo, and ensure scripts are executable:

	git clone https://github.com/SrivatsanLab/SPC_genome
	
	chmod -R u+x  SPC_genome

The main script for preprocessing SPC WGS data will split your input fastq's into a desired number of chunks, and process them in parallel using SLURM job arrays: 

	SPC_genome/WGS_PP.sh -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome.fa> -r <read_count>

	Required arguments:
	  -o    <output_name>           Desired sample name (prefix for outputs)
	  -1    <read1.fastq.gz>        Read 1 FASTQ file
	  -2    <read2.fastq.gz>        Read 2 FASTQ file
	  -g    <reference_genome>      Path to directory containing and genome fasta, fasta index, and BWA index folder
	  -r    <READ_COUNT>            Number of reads (from sequencing run info)

	Optional arguments:
	  -s    <scripts_DIR>           path to the SPC_genome directory (default: ./SPC_genome)
	  -O    <output_dir>            desired output directory (default: .)
	  -n    <N_CHUNKS>              Number of subjobs for SLURM arrays. Default: 500
	  -t    <TMP_DIR>               Temp directory for fastq chunks. Provide in order to override mktemp (e.g. when scratch space is limited)
	  -h                            Show this help message and exit
	  
The subdirectory `SPC_genome/scr` contains scripts used for each sub job, and other useful utilities.