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

vcf = VCF(vcf_file)

# First: load all variants into a list
variants_list = list(vcf)

# Extract variant names from that
variants = [f"{v.CHROM}-{v.start}-{v.REF}>{v.ALT[0]}" for v in variants_list]
variants = np.array(variants)

cells = np.array(vcf.samples)

depth = np.zeros((variants.shape[0],cells.shape[0]))
alt = np.zeros((variants.shape[0],cells.shape[0]))
binary = np.zeros((variants.shape[0],cells.shape[0]))

print('setup complete, counting reads now:')
print(f"there are {len(cells)} cells and {len(variants)} variants")

with tqdm(total=len(variants),desc="Processing Variants", unit="var") as pbar:
    for count,variant in enumerate(variants_list):
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

binary = csr_matrix(binary.T, dtype=np.float32)
adata = ad.AnnData(binary)
adata.obs_names = cells
adata.var_names = variants

adata.layers['DP']=csr_matrix(depth.T, dtype=np.float32)
adata.layers['AD']=csr_matrix(alt.T, dtype=np.float32)

adata.write_h5ad(output_name)