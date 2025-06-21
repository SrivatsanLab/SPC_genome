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

input_dir = sys.argv[1] # Directory containing your .h5ad files

# Load all AnnData objects
adatas = []
for fname in sorted(os.listdir(input_dir)):
    if fname.endswith(".h5ad"):
        adata = ad.read_h5ad(os.path.join(input_dir, fname))
        adatas.append(adata)

# Concatenate along variables (features) — axis=1 means concatenate columns
adata = ad.concat(adatas, axis=1, join="outer", merge="unique")