# CapWGS benchmarking: false-positive rate and dropout vs. public WGA methods

## Goal

Recreate Figures 3 and 4 of Gonzalez-Pena et al. (PNAS 2021,
[10.1073/pnas.2024176118](https://www.pnas.org/doi/10.1073/pnas.2024176118)) for
CapWGS HSCs, comparing them to 100M-read downsampled public single-cell WGA
data so depth is matched.

* **Fig 3** — per-cell sensitivity, false-positive rate, allele-dropout, GIV.
* **Fig 4** — mutation accumulation in CD34+ HSCs (treatment series).

## Donor groups

Built by `scripts/benchmarking/build_donor_manifest.py` →
`data/benchmarking/manifests/donor_manifest.tsv` (149 samples).

| donor_id | source | cells | bulk(s) | matched to HSC4? |
|---|---|---|---|---|
| `HSC4` | this study (CapWGS) | 50 (CapWGS) | `bulkHSC4.bqsr.bam` | self |
| `CD34_cord_blood` | PNAS PTA paper, P0 male, St. Louis Cord Blood Bank | 29 PTA (Mannitol/ENU/VHC) | SRR8438271 | ✅ best match — both HSC populations |
| `B_lymphocyte` | PTA paper, female, NIGMS | 10 PTA + 10 SCMDA | SRR8438299 (**missing on disk**) | different cell type |
| `BJ_fibroblast` | LIANTI paper, newborn male, ATCC BJ | 30 LIANTI + 3 each {MDA-Qiagen, MDA-GE, MALBAC-Yikon, MALBAC-Rubicon, DOP-PCR} | SRR5365377, SRR5365378 | different cell type |

Notes on the metadata fix: in `LIANTI_meta.tsv` the human-readable label
("LIANTI Single Cell BJ1", "Bulk 1", etc.) is in the **`replicate`** column,
not `isolate` (which is literally "cell line" for every row). The manifest
builder relies on `replicate`. For PTA, blank `amplification_method` is PTA
(cord-blood cells), `'None'` flags bulk/no-amp controls, `'PTA'`/`'SCMDA'` map
straight through.

## Variant-calling pipeline

```
                            +-----------------------------+
                            |  Per-cell BAMs (100M for    |
                            |  public; native for HSC4)   |
                            +--------------+--------------+
                                           |
                       (cells)             |             (bulks)
        ----------------------------       |       --------------------------
        | preprocess_gvcf_array.sh |       |       | preprocess_bulk_call_  |
        | AddRG -> BQSR -> HC GVCF |       |       | array.sh: MarkDup ->   |
        | per cell, per donor      |       |       | BQSR -> HC variants    |
        ----------------------------       |       --------------------------
                       \                                      /
                        \                                    /
                         v                                  v
                     +------------------+        +-------------------+
                     | per-donor GVCFs  |        | per-donor bulk    |
                     | data/benchmarking|        | truth VCF         |
                     | /gvcfs/<donor>/  |        | data/benchmarking |
                     +--------+---------+        | /bulks/<donor>/   |
                              |                   +-------------------+
                              v
                  +-----------------------+
                  | joint_call_donor.sh   |
                  | (array, 1 task/donor) |
                  | GenomicsDBImport ->   |
                  | GenotypeGVCFs         |
                  +-----------+-----------+
                              |
                              v
                +----------------------------+
                | joint VCFs (per donor):    |
                | data/benchmarking/joint/   |
                | <donor>/<donor>_joint.vcf  |
                +----------------------------+
```

* Bulks are **not** joint-called with cells — they're the independent truth
  set, so joint-calling them with cells would inflate cell sensitivity and
  bias the comparison.
* The cell pipeline takes the existing markdup'd 100M BAMs as input (markdup
  was applied by `markdup_and_metrics_array.sh` in an earlier step). The bulk
  pipeline starts from the raw full-depth BAMs and runs MarkDuplicates first.
* HSC4 bulk already has a single-sample VCF (`bulkHSC4.vcf.gz`) and HSC4 cells
  already have GVCFs in `data/HSC4/sc_outputs/*.g.vcf.gz`. The HSC4 cell GVCFs
  are symlinked into `data/benchmarking/gvcfs/HSC4/` so `joint_call_donor.sh`
  picks them up uniformly.

## Scripts

* `scripts/benchmarking/build_donor_manifest.py` — generates the donor manifest.
* `scripts/benchmarking/preprocess_gvcf_array.sh` — per-cell BQSR + HC GVCF.
  Adapted from `sc_extract_preprocess_qc_array.sh` + `sc_var_array_gatk.sh`,
  with the read-extraction and MarkDup steps removed (already done upstream).
* `scripts/benchmarking/preprocess_bulk_call_array.sh` — bulk MarkDup + BQSR +
  HC (single-sample, variants-only) for the public bulks.
* `scripts/benchmarking/joint_call_donor.sh` — per-donor joint calling.
  Adapted from `gatk_genomicsdb_import_array.sh` + `gatk_joint_calling_parallel_array.sh`,
  collapsed to a single whole-genome job (cohorts are 10–50 cells; the 100-way
  scatter isn't worth the overhead).

## Open prerequisite: HSC4 depth matching

HSC4 cells have native depths of 170M – 1.7B reads (mean 488M); public cells
are downsampled to 100M. For a strict 100M-matched comparison, HSC4 cells need
to be downsampled to 100M too. The infrastructure exists
(`scripts/benchmarking/downsample_array.sh`,
`scripts/benchmarking/markdup_and_metrics_array.sh`,
`scripts/benchmarking/create_downsample_tasks.py`) but the array hasn't been
submitted for HSC4. Until that's done we'll have to either:

1. Compare HSC4 (native) vs. public (100M) and report a depth caveat, or
2. Add an HSC4 downsample step before this pipeline. With that path, we'd
   re-run `preprocess_gvcf_array.sh` on the HSC4 100M BAMs and replace the
   symlinks in `data/benchmarking/gvcfs/HSC4/`.

Initial joint calling uses option 1.

## Open prerequisite: missing B_lymphocyte bulk

`SRR8438299` (the only B_lymphocyte bulk in the PTA paper) was never
downloaded — its `.bam` is absent from `data/public_WGA/`. Joint calling for
the B_lymphocyte donor still runs (PTA vs. SCMDA cell-only metrics are
informative), but bulk-vs-cell sensitivity / FPR / dropout is not computable
for that donor unless the BAM is fetched from SRA later.

## Open prerequisite: SRR5365332 below 100M

`SRR5365332` has only 95M mapped reads in the original BAM, so it can't be
downsampled to 100M. It is dropped from the cell task list (94 cells instead
of 95).

## Downstream analysis (post joint call)

1. **Hard-filter cell joint VCFs** (per GATK best practices):
   * SNV: `QD<2 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8`
   * Indel: `QD<2 || FS>200 || ReadPosRankSum<-20`
2. **Hard-filter bulk truth VCFs** with the same filters; restrict to
   bulk-callable regions (DP ≥ 10).
3. **Per-cell metrics vs bulk truth** (per GonzalezPena Fig 3):
   * **Sensitivity** = `(cell PASS ∩ bulk PASS in callable region)` / `bulk PASS`
   * **FPR** = `(cell PASS ∉ bulk and bulk has DP≥10 + REF call)` / `cell PASS`
   * **Allele dropout** = `(bulk het called as hom-ref or hom-alt or no-call in cell)` / `bulk het`
   * **GIV** = strand-specific transition ratios per `bcftools stats`
4. **Mutation accumulation (Fig 4)** for CD34_cord_blood:
   * Restrict to private SNVs (singletons across the donor), exclude bulk
     variants, exclude blacklisted regions, intersect callable masks across
     cells. Tally per-cell SNV burden by treatment level.
5. **Comparison plots** at 100M:
   * HSC4 (CapWGS) vs CD34_cord_blood (PTA) per-cell sensitivity / FPR /
     dropout.
   * Cross-method panel using `BJ_fibroblast` (LIANTI/MDA/MALBAC/DOP-PCR) and
     `B_lymphocyte` (PTA vs SCMDA) for context.

## Job execution

Submission script: `scripts/benchmarking/submit_variant_calling.sh`.
Submits in this order with SLURM dependencies:

1. `preprocess_bulk_call_array.sh` (4 tasks — public bulks)
2. `preprocess_gvcf_array.sh` (96 tasks — public cells at 100M)
3. Symlink HSC4 cell GVCFs → `data/benchmarking/gvcfs/HSC4/`
4. `joint_call_donor.sh` (4 tasks, one per donor; depends on step 2)

Outputs:

* `data/benchmarking/gvcfs/<donor>/<sample>.g.vcf.gz`
* `data/benchmarking/bulks/<donor>/<sample>.vcf.gz`
* `data/benchmarking/joint/<donor>/<donor>_joint.vcf.gz`
