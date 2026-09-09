project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library(ggtree)
library(tidyverse)
library(Biostrings)  # for reverseComplement
library(data.table)
set.seed(50)

setwd(file.path(project_root, "paper_figures/data/Single_cell_bottlenecking_summary_statistics/"))

input_file <- if (file.exists("filtered_bulk_vaf.csv")) {
  "filtered_bulk_vaf.csv"
} else if (file.exists("filtered_bulk_vaf.csv.gz")) {
  "filtered_bulk_vaf.csv.gz"
} else {
  stop("Missing filtered_bulk_vaf.csv(.gz) in paper_figures/data/Single_cell_bottlenecking_summary_statistics")
}
variants = fread(input_file, sep = ",")
head(variants)

variants =
  variants %>%
  mutate(mutation = stringr::str_split_fixed(V1,pattern = "-",3)[,3])


# Output for 
# Project	Sample	ID	Genome	mut_type	chrom	pos_start	pos_end	ref	alt	Type
# BRCA	example1	.	GRCh37	SNP	X	3332848	3332848	C	T	SOMATIC
# BRCA	example1	.	GRCh37	SNP	X	3358662	3358662	G	A	SOMATIC
# BRCA	example1	.	GRCh37	SNP	X	3894402	3894402	G	A	SOMATIC
# BRCA	example1	.	GRCh37	SNP	X	4071484	4071484	C	T	SOMATIC

variants_parsed <- 
  variants %>%
  separate(V1, into = c("chrom", "pos", "change"), sep = "-", remove = FALSE) %>%
  separate(change, into = c("ref", "alt"), sep = ">", remove = FALSE) %>%
  mutate(
    Project = "CapWGS",
    Sample = "scTree",
    ID = ".",
    Genome = "GRCh38",
    mut_type = "SNP",
    pos_start = as.integer(pos),
    pos_end = as.integer(pos),
    Type = "SOMATIC",
    chrom_short = stringr::str_replace(chrom,"chr","")
  ) %>%
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)


write.table(variants_parsed, "test/input/filtered_variants.txt", sep = "\t", quote = FALSE, row.names = FALSE)








