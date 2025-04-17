#!/bin/usr/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from itertools import product, combinations
from cyvcf2 import VCF
from tqdm import tqdm
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
import sys

vcf_file = sys.argv[1]
output_name = sys.argv[2]

depth = np.zeros((variants.shape[0],cells.shape[0]))
alt = np.zeros((variants.shape[0],cells.shape[0]))
binary = np.zeros((variants.shape[0],cells.shape[0]))

vcf_file = "vcf_chunks/chr10:100000001-105000000.vcf.gz"
vcf = VCF(vcf_file)
with tqdm(total=len(variants),desc="Processing Variants", unit="var") as pbar:
    for count,variant in enumerate(vcf):
        # print(variant)
        alt_dp = variant.gt_alt_depths
        ref_dp = variant.gt_ref_depths
        tot_dp = np.add(alt_dp,ref_dp)
        # store alt depth
        alt[count,]+=alt_dp
        # store total depth
        depth[count,]+=tot_dp
        
        # mask = (depth[count,] > 0) & (alt[count,] > 0)
        mask = (variant.gt_types != 0)
        binary[count,mask] = 1
        
        pbar.update(1)

vcf = VCF(vcf_file)
cells = np.array(vcf.samples)
variants = [f"{variant.CHROM}-{variant.start}-{variant.REF}>{variant.ALT[0]}" for variant in vcf]
variants = np.array(variants)

binary = csr_matrix(binary.T, dtype=np.float32)
adata = ad.AnnData(binary)
adata.obs_names = cells
adata.var_names = variants

adata.write_h5ad(output_name)