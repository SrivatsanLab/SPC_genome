#!/usr/bin/env python3
"""Build expanded L4-appropriate marker gene sets for C. elegans cell types.

Sources:
  - Cao et al. 2017 (whole-animal L2 sci-RNA-seq, 117 cell types, WormBase-wrangled
    h5ad from https://data.caltech.edu/records/j4ycn-dpv05). Used for all tissue
    markers, including "neuron" (as the aggregate of every cell type not in the
    tissue map).

Method:
    1. Aggregate Cao's 117 fine cell types into 14 broad tissues (TISSUE_MAP);
       anything not mapped becomes 'neuron'.
    2. Log-normalize + Wilcoxon rank_genes_groups per tissue vs rest, pulling
       top CANDIDATES=200 per tissue.
    3. Specificity filter: keep only genes that appear in <= MAX_SHARED tissues'
       top-CANDIDATES lists. This removes cross-tissue bleed (e.g. collagens
       appearing in both seam_cell and hypodermis).
    4. Emit top TOP_N=40 filtered markers per tissue as WBGene IDs.

Output: <output.json>  {tissue: [WBGeneID, ...]}

Usage:
    scripts/CapGTA/build_expanded_markers.py <cao2017.h5ad> <output.json>
"""

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc


TISSUE_MAP = {
    'body_muscle':           ['BWM_middle_and_posterior','BWM_neck','BWM_posterior','BWM_anterior','BWM_far_posterior','BWM_head'],
    'pharyngeal_muscle':     ['pm1_pm2','pm3_pm4_pm5','pm6_and_possibly_pm7'],
    'intestine':             ['Intestine_anterior','Intestine_middle_and_posterior','Intestine_far_posterior'],
    'hypodermis':            ['hyp_4_to_7_bin_1_mid_L2','hyp_4_to_7_bin_2_late_L2','hyp_4_to_7_bin_3_around_L2_molt','Head_hypodermis','Tail_hypodermis_bin_1_late_L2','Tail_hypodermis_bin_2_around_L2_molt'],
    'seam_cell':             ['Seam_cell_bin_1_mid_L2','Seam_cell_bin_2_late_L2','Seam_cell_bin_3_around_L2_molt','Seam_cell_bin_4_ecdysis_post_L2_molt','Seam_cell_bin_5_early_L3'],
    'germline':              ['Germline'],
    'coelomocyte':           ['Coelomocyte'],
    'glia':                  ['AMsh','AMso_PHso','Includes_CEPso_ILso','CEPsh','PHsh','GLR'],
    'excretory':             ['Excretory_gland','Excretory_duct','Excretory_cell'],
    'somatic_gonad':         ['Somatic_gonad_precursor','Distal_tip_cell','Vulval_precursor_bin_1_late_L2','Vulval_precursor_bin_2_early_L3','Sex_myoblast'],
    'pharyngeal_epithelium': ['Anterior_arcade_bin_1_mid_to_late_L2','Anterior_arcade_bin_2_around_L2_molt','Anterior_arcade_bin_3_early_L3','Posterior_arcade_bin_1_mid_to_late_L2','Posterior_arcade_bin_2_around_L2_molt','Posterior_arcade_bin_3_early_L3','Pharyngeal_epithelial_cell_tentative','Marginal_cell_bin_1_mid_to_late_L2','Marginal_cell_bin_2_around_L2_molt_to_early_L3'],
    'pharyngeal_gland':      ['g1A','g1P','g2'],
    'rectal_mesoderm':       ['mu_int_mu_anal','mu_sph','XXX','hmc'],
}
TOP_N = 40
CANDIDATES = 200
MAX_SHARED = 2


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('cao_h5ad', type=Path, help='Cao 2017 wrangled h5ad')
    p.add_argument('output_json', type=Path, help='output JSON path')
    args = p.parse_args()

    if not args.cao_h5ad.is_file():
        print(f'Error: {args.cao_h5ad} not found', file=sys.stderr)
        return 1

    a = ad.read_h5ad(args.cao_h5ad)
    a = a[a.obs['cell_type'].notna() & (a.obs['cell_type'].astype(str) != 'nan')].copy()

    reverse = {ct: tissue for tissue, cts in TISSUE_MAP.items() for ct in cts}
    a.obs['tissue'] = a.obs['cell_type'].map(reverse).fillna('neuron').astype('category')

    print('Tissue cell counts:')
    print(a.obs['tissue'].value_counts())

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, 'tissue', method='wilcoxon', n_genes=CANDIDATES)

    tissues = list(a.obs['tissue'].cat.categories)
    cand = {t: list(a.uns['rank_genes_groups']['names'][t][:CANDIDATES]) for t in tissues}

    appear_in: dict[str, list[str]] = {}
    for t, genes in cand.items():
        for g in genes:
            appear_in.setdefault(g, []).append(t)

    markers = {t: [g for g in genes if len(appear_in[g]) <= MAX_SHARED][:TOP_N] for t, genes in cand.items()}

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(markers, indent=2))
    print(f'\nWrote {args.output_json} ({len(markers)} tissues)')

    name_map = dict(zip(a.var.index, a.var['gene_name']))
    print('\nTop 10 markers per tissue:')
    for t, wb in markers.items():
        named = [name_map.get(w, w) for w in wb[:10]]
        print(f"  {t:24s} {', '.join(named)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
