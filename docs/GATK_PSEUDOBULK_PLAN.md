# Implementation plan: per-worm pseudobulk variant calling with GATK

**For an agent executing on a SLURM cluster (Fred Hutch, this repo).** Read section 0
before starting — this task has a specific failure mode that will silently corrupt
results if not handled.

---

## 0. Context and the one thing you must not get wrong

### Project specifics (this repo)

| item | value |
|---|---|
| assay | **CapGTA** (genome + transcriptome coassay), not CapWGS — BAMs are STAR-aligned and contain **both DNA and RNA reads** (see section 3) |
| reference | `/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa` (same one STAR aligned against; `.fai` and `.dict` already present) |
| h5ad with assignments | `results/worm6_final/DNA_analysis/joint_variants_process.h5ad` (6 GB — open with `backed='r'`) |
| donor column to use | `adata.obs["donor"]` from the K=20, ncomp=16, all_variants sweep run — this is what got written to the top-level `obs.csv` and the processed h5ad |
| per-cell BAM paths | already stored in `adata.obs["bam_path"]` — use it directly, don't rebuild it; layout is `/fh/working/srivatsan_s/Dustin/worm6-dubnova-100k-align/runs/worm6_UDI{1..24}/results/worm6_UDI{N}/bams/{cell}.bam` |
| existing joint callset | `results/worm6_final/joint_variants/joint_variants.vcf.gz` (bcftools mpileup, 933,159 variants) and `results/worm6_final/joint_variants/regions.txt` (104 × 1 Mb regions — **reuse for scatter**) |
| output root | **`results/worm6_final/pseudobulk_gatk/`** — all deliverables live under here |
| scripts location | **`scripts/CapGTA/pseudobulk_gatk/`** (project convention: SLURM scripts under `scripts/{pipeline}/`) |
| SLURM logs | `SLURM_outs/array_outs/%x_%A_%a.out` (project convention) |
| environment | `micromamba activate spc_genome` for Python; `module load GATK / SAMtools / BCFtools` for tools |

### What this is

~2,400 single-cell whole-genome amplification libraries from a pool of ~16–20
*C. elegans* individuals. The individuals are selfing descendants of one founder
carrying a **pole-1 P270R hypermutator** allele, so they are effectively
**homozygous** at germline sites. Cells have already been demultiplexed to
inferred worms of origin by an SVD/haplotype-module method; assignments are in
`adata.obs["donor"]` (values `worm01`…`wormNN`, plus `background`).

The goal here is to merge each inferred worm's BAMs into a pseudobulk, call
variants with GATK HaplotypeCaller, and joint-genotype across worms.

### Why: three specific objectives, in priority order

1. **Resolve germline vs somatic by zygosity.** Individual cells are ~1×, so
   zygosity is unmeasurable. Pooled, a worm reaches 10–150×, where hom and het
   separate cleanly. Germline in a selfing line is **homozygous**; somatic
   mutations are **heterozygous**. This directly tests an open question: whether
   inferred modules beyond ~20 are additional worms or sub-worm somatic clones.
2. **Discover variants pileup calling misses.** HaplotypeCaller performs local
   de novo reassembly and recovers indels and variants in misaligned regions
   that the existing bcftools joint callset does not contain.
3. **Build a per-worm germline blacklist** for downstream somatic mutation
   calling.

**This is NOT to find more sites for demultiplexing.** The demultiplexing is not
site-limited (median ~27,000 jointly-covered informative sites per cell pair,
against ~150 needed). Do not evaluate success by counting sites.

### THE CRITICAL PITFALL — contamination breaks zygosity calls

Every capsule contains ambient ("soup") DNA from the whole pool. Contamination
fraction per cell is estimated as `adata.obs["alpha_hap"]` (= 1 − purity) and is
**large and graded** — commonly 0.2–0.6, not a small nuisance.

In a pooled BAM, a **homozygous** germline site therefore reads at

```
VAF ≈ (1 − α)·1 + α·(1/K)
```

| α | pooled VAF at a hom-alt germline site | GATK's diploid model calls |
|---|---|---|
| 0.0 | 1.00 | 1/1 ✓ |
| 0.2 | 0.81 | **0/1 — WRONG** |
| 0.4 | 0.62 | **0/1 — WRONG** |
| 0.6 | 0.43 | **0/1 — WRONG** |

**Consequence: if you trust GATK's `GT` field, every homozygous germline site in
every worm will be miscalled heterozygous, which inverts the exact discriminator
this analysis exists to measure.**

**Required mitigation — do all three:**

- Use HaplotypeCaller for **allele discovery** (which sites, which alleles).
- **Re-genotype yourself from the AD/DP fields**, not from `GT`, using a
  contamination-aware threshold (section 5).
- Pass a per-worm `--contamination-fraction-to-filter` estimated from
  `alpha_hap` so GATK's own filtering is at least approximately right.

---

## 1. Inputs (already located — verify before running)

All paths absolute unless noted. Verify each before dispatch:

```bash
REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
H5AD=$REPO/results/worm6_final/DNA_analysis/joint_variants_process.h5ad
REF=/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa
REGIONS=$REPO/results/worm6_final/joint_variants/regions.txt      # reuse — 104 × 1 Mb
EXISTING_VCF=$REPO/results/worm6_final/joint_variants/joint_variants.vcf.gz
OUT=$REPO/results/worm6_final/pseudobulk_gatk
SCRIPTS=$REPO/scripts/CapGTA/pseudobulk_gatk

test -s "$REF.fai" && test -s "${REF%.fa}.dict"  # both present
test -s "$H5AD" && test -s "$REGIONS" && test -s "$EXISTING_VCF"
```

`adata.obs` provides `donor`, `capsule_type`, `purity`, `alpha_hap`,
`n_core_sites`, `sample` (UDI), `bam_path`. `donor` values are
`worm01..worm20` plus `background`. Do NOT ask the user for per-cell BAM paths
— they are already in `bam_path`. Spot-check a handful (`ls $(head -3 ...)`)
before merging.

---

## 2. Build per-worm cell lists

Output goes to `$OUT/worm_lists/`. Each worm gets **two** files: `<worm>.cells.txt`
(barcodes, for diagnostics) and `<worm>.bampaths` (BAM paths, feeds `samtools merge -b`).

```python
import anndata as ad
import pandas as pd
from pathlib import Path

H5AD = "results/worm6_final/DNA_analysis/joint_variants_process.h5ad"
OUT  = Path("results/worm6_final/pseudobulk_gatk")
(OUT / "worm_lists").mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs

# EXCLUDE background — they carry pool-average genotypes and would contaminate
# every worm's pseudobulk with every other worm's variants.
keep = obs["donor"].astype(str).str.startswith("worm")

# Permissive alpha cut only. Do NOT filter hard: excluding high-alpha cells
# biases toward clean cells and reduces depth.
keep &= obs["alpha_hap"] < 0.9

# Guard: verify bam_path exists on disk for a random sample (10 cells).
sample_paths = obs.loc[keep, "bam_path"].sample(10, random_state=0)
missing = [p for p in sample_paths if not Path(p).is_file()]
assert not missing, f"BAM paths not found on disk: {missing[:3]}"

summary = (obs[keep].groupby("donor", observed=True)
           .agg(n_cells=("purity", "size"),
                med_alpha=("alpha_hap", "median"),
                med_purity=("purity", "median"),
                med_sites=("n_core_sites", "median")))
summary["est_depth"] = summary.n_cells * 1.15 * 0.25   # rough; measure in §3
print(summary.sort_values("n_cells", ascending=False).to_string())
summary.to_csv(OUT / "worm_summary.csv")

for worm, grp in obs[keep].groupby("donor", observed=True):
    (OUT / "worm_lists" / f"{worm}.cells.txt").write_text(
        "\n".join(grp.index.astype(str)) + "\n")
    (OUT / "worm_lists" / f"{worm}.bampaths").write_text(
        "\n".join(grp["bam_path"].astype(str)) + "\n")
```

**Known composition at K=20 (from `results/worm6_final/DNA_analysis/obs.csv`):**
17 worms with ≥25 cells (worm07=186 down to worm14=25), 3 tiny modules
(worm04=6, worm02=5, worm19=2), 591 background. Expect the 3 tiny modules to
fall out at the `<50 cells` threshold — report them, then skip.

---

## 3. Merge BAMs per worm

**This is Case A** (one BAM per cell, at the paths in `obs["bam_path"]`). The
merged pseudobulks go to `$OUT/worm_bams/`.

### CapGTA-specific step: strip RNA-artifact reads before HC

The per-cell BAMs are STAR alignments of a **combined DNA + RNA** library.
Three read categories are RNA (or RNA-derived artifacts) that confound
HaplotypeCaller's local reassembly and should be removed before calling:

1. **Spliced reads** (CIGAR contains `N`, an intron skip) — bona fide RNA
   alignments. The existing `scripts/CapGTA/filter_spliced_reads_array.sh`
   **keeps** these (writes RNA-only BAMs for GEX counting) — here we want the
   complement.
2. **Poly-A reads** — reads dominated by a long `A` run, typically from mRNA
   poly(A) tails or internal-priming artifacts. These map somewhere in the
   genome (usually poorly) with no `N` in CIGAR and would seed spurious HC
   reassemblies.
3. **Poly-T reads** — same class, seen on the opposite strand.

Chain all three into one `samtools view -e` pass after merge:

```bash
samtools merge -@ 8 -f "$OUT/worm_bams/${WORM}.raw.bam" \
    -b "$OUT/worm_lists/${WORM}.bampaths"

# Filter: no splice-junction reads, no poly-A/T reads (≥20 consecutive A or T).
# samtools filter expressions use POSIX ERE; `seq` and `cigar` are addressable.
samtools view -@ 8 -b -h \
    -e '!(cigar =~ "N") && !(seq =~ "A{20,}") && !(seq =~ "T{20,}")' \
    "$OUT/worm_bams/${WORM}.raw.bam" \
    -o "$OUT/worm_bams/${WORM}.bam"
samtools index -@ 8 "$OUT/worm_bams/${WORM}.bam"
rm "$OUT/worm_bams/${WORM}.raw.bam"    # keep only DNA-only, RNA-artifact-free pseudobulk
```

**Threshold choice.** 20 consecutive A/T is conservative for the C. elegans
genome (AT-rich but 20+ homopolymer runs are rare outside intergenic regions).
Drop to 15 if downstream inspection still shows RNA-tail contamination; raise
to 25 if genuine genomic AT-runs are being lost. Log the fraction of reads
excluded by this filter per worm so the choice is auditable:

```bash
n_in=$(samtools view -c "$OUT/worm_bams/${WORM}.raw.bam")
n_out=$(samtools view -c "$OUT/worm_bams/${WORM}.bam")
echo "${WORM} kept ${n_out}/${n_in} reads ($(bc -l <<< "scale=4; ${n_out}/${n_in}"))" \
    >> "$OUT/filter_report.tsv"
```

### Read groups and duplicates — already handled upstream

Confirmed from a per-cell BAM header (`worm6_UDI1/A10_A7_B4_F2_...bam`):

- **`@RG`** — one per cell, `ID=SM=LB=<cell_barcode>`, `PL=ILLUMINA`. `samtools
  merge` preserves all of them. GATK handles many read groups fine; leave them.
- **`MarkDuplicates`** — already run per cell by
  `scripts/CapWGS/sc_extract_preprocess_qc_array.sh` (via Picard, `REMOVE_DUPLICATES=false`,
  duplicates flagged). Do NOT re-run at the worm level: cross-cell "duplicates"
  from different WGA capsules are independent observations, not PCR duplicates,
  and MarkDuplicates would wrongly collapse them.
- **Sample tag for the joint VCF** — the per-cell `SM` tags will make GATK
  emit hundreds of sample columns. To collapse to one column per worm, either
  (a) rewrite RG SM tags to the worm name after merge:
  ```bash
  samtools addreplacerg -@ 8 -r "ID:${WORM}\tSM:${WORM}\tLB:${WORM}\tPL:ILLUMINA" \
      -m overwrite_all -o "$OUT/worm_bams/${WORM}.bam.tmp" "$OUT/worm_bams/${WORM}.bam"
  mv "$OUT/worm_bams/${WORM}.bam.tmp" "$OUT/worm_bams/${WORM}.bam"
  samtools index -@ 8 "$OUT/worm_bams/${WORM}.bam"
  ```
  or (b) merge with `-r` / `-p` flags that rewrite RGs on the fly. Choose (a).

### Measure actual depth per worm before calling

```bash
# Reference chromosomes only, skip flagged duplicates and secondary alignments.
samtools depth -a -Q 20 -G 0x400 "$OUT/worm_bams/${WORM}.bam" \
  | awk '{s+=$3; n++; if($3>0){c++; cs+=$3}} END {
      print "mean depth all sites:", s/n;
      print "breadth (frac covered):", c/n;
      print "mean depth at covered sites:", cs/c}' \
  > "$OUT/depth_report/${WORM}.txt"
```

Aggregate into `$OUT/depth_report.tsv`. This determines which worms can be
genotyped at all (§8, failure mode 2).

---

## 4. HaplotypeCaller in GVCF mode, then joint genotyping

Reuse the existing scatter — `results/worm6_final/joint_variants/regions.txt`
already partitions WBcel235 into 104 × 1 Mb regions and is the same partition
the current bcftools joint callset uses. Convert to GATK `-L` on the fly
(`chrom:start-end` is native GATK syntax; no conversion needed).

**SLURM array shape:** `(n_worms × n_regions)`. With 17 worms passing the
`≥50 cells` bar × 104 regions = 1,768 tasks. Run as one giant array with
`SLURM_ARRAY_TASK_ID` decoded to `(worm_i, region_j)` — mirrors
`scripts/CapWGS/bcftools_joint_calling_parallel_array.sh`. Emit logs to
`SLURM_outs/array_outs/` per project convention.

### Per-worm × per-region GVCF

```bash
REF=/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa
OUT_ROOT=results/worm6_final/pseudobulk_gatk
REGION=$(sed -n "${REGION_IDX}p" results/worm6_final/joint_variants/regions.txt)
ALPHA=$(awk -F, -v w="$WORM" '$1==w {print $3}' "$OUT_ROOT/worm_summary.csv")  # med_alpha

module load GATK

gatk --java-options "-Xmx14g" HaplotypeCaller \
    -R "$REF" \
    -I "$OUT_ROOT/worm_bams/${WORM}.bam" \
    -O "$OUT_ROOT/gvcf/${WORM}/${REGION//[:\-]/_}.g.vcf.gz" \
    -L "$REGION" \
    -ERC GVCF \
    --sample-ploidy 2 \
    --contamination-fraction-to-filter "$ALPHA" \
    --min-base-quality-score 20
```

- `--sample-ploidy 2` — *C. elegans* hermaphrodite is XX diploid. Do not use
  pooled/high ploidy; the point is to resolve zygosity in a diploid individual.
- `--contamination-fraction-to-filter ${ALPHA}` — the worm's **median `alpha_hap`**
  from `worm_summary.csv`. Per-worm, not global. This is only a coarse GATK-side
  filter; the real contamination correction is in §5.
- Region names contain `:` and `-` which are illegal in some filesystems; substitute
  to `_` for the output filename.

### Combine + joint-genotype per region

```bash
gatk CombineGVCFs -R "$REF" \
    $(printf -- '-V %s ' "$OUT_ROOT"/gvcf/*/"${REGION//[:\-]/_}".g.vcf.gz) \
    -O "$OUT_ROOT/combined_gvcf/${REGION//[:\-]/_}.g.vcf.gz"

gatk GenotypeGVCFs -R "$REF" \
    -V "$OUT_ROOT/combined_gvcf/${REGION//[:\-]/_}.g.vcf.gz" \
    -O "$OUT_ROOT/joint/${REGION//[:\-]/_}.vcf.gz" \
    --sample-ploidy 2
```

At ~17 worms `CombineGVCFs` is fine (`GenomicsDBImport` only wins above ~30
samples). Joint genotyping across worms is important: it forces every worm to be
genotyped at every site, giving a complete sites × worms matrix rather than a
ragged union of per-worm calls.

### Concat per-region VCFs into one final VCF

```bash
module load BCFtools
ls "$OUT_ROOT"/joint/region_*.vcf.gz | sort > "$OUT_ROOT/joint/vcf_list.txt"
bcftools concat -a -O z -f "$OUT_ROOT/joint/vcf_list.txt" \
    -o "$OUT_ROOT/joint.vcf.gz"
bcftools index "$OUT_ROOT/joint.vcf.gz"
```

---

## 5. Re-genotype from AD/DP — do not use GATK's GT

This is the step that makes the results correct. Extract allele depths, not
genotypes:

```bash
gatk VariantsToTable -V joint.vcf.gz -O joint_table.tsv \
  -F CHROM -F POS -F REF -F ALT -F TYPE -GF AD -GF DP -GF GQ
```

```python
import numpy as np, pandas as pd

tab = pd.read_csv("joint_table.tsv", sep="\t")
alpha = pd.read_csv("worm_summary.csv", index_col=0)["med_alpha"]
worms = [c[:-3] for c in tab.columns if c.endswith(".AD")]

MIN_DP = 10          # per-worm depth required to genotype at all
rows = {}
for w in worms:
    adw = tab[f"{w}.AD"].astype(str).str.split(",", expand=True)
    ref = pd.to_numeric(adw[0], errors="coerce")
    alt = pd.to_numeric(adw[1], errors="coerce")
    dp = ref + alt
    vaf = alt / dp.replace(0, np.nan)
    a = float(alpha.get(w, 0.4))
    # a HOMOZYGOUS site in a capsule pool with contamination a reads at ~(1-a);
    # a HETEROZYGOUS site reads at ~(1-a)/2. Put the cut between them.
    hom_lo = 0.75 * (1 - a)
    het_lo = 0.25 * (1 - a)
    gt = np.where(dp < MIN_DP, "./.",
         np.where(vaf >= hom_lo, "1/1",
         np.where(vaf >= het_lo, "0/1", "0/0")))
    rows[f"{w}_vaf"] = vaf
    rows[f"{w}_dp"] = dp
    rows[f"{w}_gt"] = gt
geno = pd.concat([tab[["CHROM","POS","REF","ALT","TYPE"]], pd.DataFrame(rows)], axis=1)
geno.to_csv("regenotyped.tsv", sep="\t", index=False)
```

**Sanity check that must pass before interpreting anything:** plot the VAF
histogram per worm at sites where that worm is non-reference. It should be
**bimodal**, with modes near `(1−α)` and `(1−α)/2`. If it is unimodal, either
the contamination estimate is wrong or the pseudobulk is mixing worms — stop and
report rather than proceeding.

---

## 6. Analyses to run and report

### 6.1 Germline vs somatic per worm (objective 1)

```
per worm: n hom sites, n het sites, hom:het ratio
```

Germline in a selfing line is homozygous. A large het fraction is either somatic
load (expected under a hypermutator) or residual heterozygosity from recent
mutations not yet fixed by selfing. Report both counts.

### 6.2 Are the extra modules worms or sub-clones? (the open question)

For each inferred worm, compute the fraction of its **private** variants (called
non-reference in that worm and reference in all others) that are **homozygous**.

- A genuine worm: private variants are largely **hom** → it is a distinct
  germline lineage.
- A somatic sub-clone wrongly split from a worm: private variants are largely
  **het**, and the module will share most of its *homozygous* sites with its
  parent worm.

Also compute the pairwise sharing matrix of **homozygous** calls between worms.
Two modules from the same worm will share nearly all hom sites. Report any pair
sharing >90% — those should be merged, and their count reduces the effective K.

### 6.3 New sites vs the existing callset (objective 2)

Intersect with the original bcftools VCF. Report counts of sites gained, broken
down by `TYPE` (SNP vs INDEL). Indels are where reassembly is expected to help
most.

### 6.4 Per-worm germline blacklist (objective 3)

Write, per worm, a BED/VCF of positions called non-reference in **any** worm.
This is the exclusion list for downstream somatic calling. Note it is
**necessary but not sufficient** — the ambient pool also carries other cells'
somatic mutations, so read-count thresholds are still required at the somatic
calling stage.

---

## 7. Scripts to create

Under `scripts/CapGTA/pseudobulk_gatk/` (project convention):

| script | purpose | invocation |
|---|---|---|
| `build_worm_lists.py` | reads h5ad, writes `worm_summary.csv` + `worm_lists/*.{cells,bampaths}` | `python … build_worm_lists.py $H5AD $OUT` |
| `merge_pseudobulk_array.sh` | SLURM array (one task per worm): merge → strip spliced → rewrite RG SM → index → depth | `sbatch --array=1-N …` |
| `haplotype_caller_array.sh` | SLURM array (worms × regions): per-worm GVCF over one region, `--contamination-fraction-to-filter ${ALPHA}` | `sbatch --array=1-1768 …` |
| `combine_gvcfs_array.sh` | SLURM array (one task per region): `CombineGVCFs` → `GenotypeGVCFs` | `sbatch --array=1-104 …` |
| `concat_and_final.sh` | one-shot: `bcftools concat` per-region VCFs → `joint.vcf.gz` | `sbatch …` |
| `regenotype_from_ad.py` | `VariantsToTable` → contamination-aware GT from AD/DP → `regenotyped.tsv` + VAF histograms | `python … regenotype_from_ad.py` |
| `analyze_worms_vs_subclones.py` | hom:het ratios, pairwise hom-sharing matrix, new-sites-vs-existing, blacklist BED, REPORT.md | `python … analyze_worms_vs_subclones.py` |

Wire them into `scripts/CapGTA/CapGTA_pseudobulk_gatk_pipeline.sh` — dispatches
the arrays with the right `--dependency=afterok:...` chain.

## 8. Deliverables (all under `results/worm6_final/pseudobulk_gatk/`)

```
results/worm6_final/pseudobulk_gatk/
├── worm_summary.csv          cells, median alpha, estimated + MEASURED depth per worm
├── worm_lists/
│   ├── {worm}.cells.txt      cell barcodes
│   └── {worm}.bampaths       BAM paths (input to samtools merge)
├── worm_bams/
│   ├── {worm}.bam            merged, spliced-filtered, RG-rewritten pseudobulk
│   └── {worm}.bam.bai
├── filter_report.tsv         per-worm reads-kept fractions after splice + polyA/T filter
├── depth_report/{worm}.txt   per-worm samtools depth summary
├── depth_report.tsv          aggregated
├── gvcf/{worm}/{region}.g.vcf.gz     per-worm × per-region GVCFs
├── combined_gvcf/{region}.g.vcf.gz   combined across worms, per region
├── joint/{region}.vcf.gz             joint-genotyped per region
├── joint.vcf.gz                      final concat (+ .csi)
├── regenotyped.tsv                   contamination-aware genotypes from AD/DP
├── analyses/
│   ├── vaf_histograms.png            per-worm VAF distribution — correctness check
│   ├── hom_sharing_matrix.csv        pairwise hom-site sharing between worms
│   ├── new_sites_summary.txt         sites gained vs existing callset, by TYPE
│   └── germline_blacklist.bed        union of non-reference positions
└── REPORT.md                         findings, with explicit caveats
```

SLURM logs go to the repo-standard `SLURM_outs/array_outs/` (do NOT nest inside
the output tree).

---

## 9. Failure modes to watch for and report rather than work around

1. **Unimodal VAF histogram** — pseudobulk is mixing worms, or `alpha_hap` is
   wrong. Do not proceed to zygosity conclusions.
2. **A worm with <10× measured depth** — cannot be genotyped; report as
   inconclusive rather than calling it.
3. **Two worms sharing >90% of homozygous sites** — they are one individual
   split in two, or siblings too closely related to resolve. This is a *finding*
   about K, not an error to fix.
4. **Duplicate rate very high after MarkDuplicates** — WGA libraries can be
   dominated by duplicates; if effective depth after deduplication is far below
   the estimate, report it before running the full callset.
5. **Barcode-to-BAM mapping mismatch** — verify on a subset first. This is
   silent and catastrophic.

## 10. Interpretive caveats to carry into REPORT.md

- Cell assignments are **inferred**, so every pseudobulk inherits assignment
  error. Do not present per-worm callsets as ground truth.
- This analysis **cannot validate the demultiplexing** that produced its input —
  the pseudobulks are defined by the assignment. It can flag internal
  inconsistency (worms sharing hom sites) but not confirm correctness.
- Contamination is **graded, not discrete**: there is no clean cell/empty split
  in this data. A single per-worm α is a summary of a distribution; sites near
  the hom/het boundary are genuinely ambiguous.
- Worms differ greatly in cell count (observed range ~539 to ~7), so per-worm
  depth and therefore sensitivity differ by more than an order of magnitude.
  Never compare raw variant counts between worms without depth normalization.
