#!/bin/bash
#SBATCH -J sc_align_qc_test
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 0:30:00
#SBATCH -p short

##########################################################################################################################
# Per-cell STAR alignment + QC for the new single-cell-fastq demux test.
# One array task = one cell (input fastqs already trimmed, CB:Z tag already in headers).
# Outputs one BAM per cell and a one-row CSV of QC metrics.
##########################################################################################################################

set -uo pipefail

CELLS_FILE="$1"            # text file, one cell ID per line
INPUT_DIR="$2"             # dir containing {cell}_{SAMPLE_SUFFIX}_R{1,2}.fastq.gz
SAMPLE_SUFFIX="$3"         # e.g. "L4-4_test3" — the bit between the cell ID and "_R1.fastq.gz"
BAM_DIR="$4"               # data/L4-4_test3/sc_outputs
QC_DIR="$5"                # data/L4-4_test3/qc_metrics
PER_CELL_CSV_DIR="$6"      # results/L4-4_test3/per_cell_qc
STAR_INDEX_DIR="$7"        # dir containing the STAR index (e.g. .../STARIndex)
REFERENCE_FASTA="$8"       # FASTA whose contig names match the STAR index (needs sibling .dict for Picard)
GTF_FILE="$9"
SCRIPTS_DIR="${10}"

mkdir -p SLURM_outs/array_outs "${BAM_DIR}" "${QC_DIR}" "${PER_CELL_CSV_DIR}"

CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CELLS_FILE}")
if [ -z "${CELL}" ]; then
    echo "ERROR: empty cell ID at line ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

R1="${INPUT_DIR}/${CELL}_${SAMPLE_SUFFIX}_R1.fastq.gz"
R2="${INPUT_DIR}/${CELL}_${SAMPLE_SUFFIX}_R2.fastq.gz"
if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "ERROR: missing fastq for ${CELL}: ${R1} / ${R2}" >&2
    exit 1
fi

if [ ! -d "${STAR_INDEX_DIR}" ]; then
    echo "ERROR: STAR index dir not found: ${STAR_INDEX_DIR}" >&2; exit 1
fi
if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "ERROR: reference FASTA not found: ${REFERENCE_FASTA}" >&2; exit 1
fi
# Picard CollectWgsMetrics needs a .dict next to the FASTA.
DICT="${REFERENCE_FASTA%.*}.dict"
if [ ! -f "${DICT}" ]; then
    echo "ERROR: Reference dict not found at ${DICT} — build it with Picard CreateSequenceDictionary." >&2
    exit 1
fi

# Per-task scratch dir (node-local /tmp — never on NFS, deleted at exit).
SCRATCH="${TMPDIR:-/tmp}/sc_align_qc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${CELL}"
mkdir -p "${SCRATCH}"
trap 'rm -rf "${SCRATCH}"' EXIT

BAM="${BAM_DIR}/${CELL}.bam"
BIGWIG="${SCRATCH}/${CELL}.bw"
DUP_METRICS="${QC_DIR}/${CELL}_duplication_metrics.txt"
WGS_METRICS="${QC_DIR}/${CELL}_wgs_metrics.txt"
ALIGN_METRICS="${QC_DIR}/${CELL}_alignment_metrics.txt"
LORENZ_CSV="${QC_DIR}/${CELL}_lorenz.csv"
GINI_TXT="${QC_DIR}/${CELL}_gini.txt"
EXONIC_TXT="${QC_DIR}/${CELL}_exonic_enrichment.txt"
PER_CELL_CSV="${PER_CELL_CSV_DIR}/${CELL}.csv"

echo "=========================================="
echo "Cell: ${CELL}  (array task ${SLURM_ARRAY_TASK_ID})"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Scratch: ${SCRATCH}"
echo "=========================================="

##########################################################################################################################
# 1. STAR alignment (splice-aware, same params as CapGTA)
##########################################################################################################################
module load STAR/2.7.9a-GCC-11.2.0

echo ">>> STAR alignment"
# --outSAMunmapped Within keeps unmapped reads in the BAM so we can compute a true
# total-vs-aligned read count from the BAM directly.
STAR --runThreadN 8 \
     --genomeDir "${STAR_INDEX_DIR}" \
     --readFilesIn "${R1}" "${R2}" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${SCRATCH}/star_" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM NM MD \
     --outSAMunmapped Within \
     --outFilterMultimapNmax 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 15000 \
     --alignMatesGapMax 15000

STAR_BAM="${SCRATCH}/star_Aligned.sortedByCoord.out.bam"
module unload STAR

module load SAMtools

# CB:Z is in the FASTQ header but STAR's default fastq path drops it.
# Add a samtools step to re-inject CB:Z from the read name — but since headers
# look like "@READID CB:Z:XXX", STAR doesn't propagate the comment. We accept this
# for the test (cell ID is already in the filename / BAM path).
cp "${STAR_BAM}" "${BAM}"
samtools index -@ 4 "${BAM}"

##########################################################################################################################
# 2. MarkDuplicates (QC only — write markdup BAM to scratch, keep original BAM)
##########################################################################################################################
module load picard

MARKDUP_BAM="${SCRATCH}/${CELL}.markdup.bam"
echo ">>> MarkDuplicates"
java -Xmx12g -jar "${EBROOTPICARD}/picard.jar" MarkDuplicates \
    INPUT="${BAM}" \
    OUTPUT="${MARKDUP_BAM}" \
    METRICS_FILE="${DUP_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES=false

##########################################################################################################################
# 3. Picard CollectAlignmentSummaryMetrics + CollectWgsMetrics (for PCT_1X)
##########################################################################################################################
echo ">>> Picard alignment + WGS metrics"
java -Xmx12g -jar "${EBROOTPICARD}/picard.jar" CollectAlignmentSummaryMetrics \
    I="${MARKDUP_BAM}" O="${ALIGN_METRICS}" R="${REFERENCE_FASTA}" || true

java -Xmx12g -jar "${EBROOTPICARD}/picard.jar" CollectWgsMetrics \
    I="${MARKDUP_BAM}" O="${WGS_METRICS}" R="${REFERENCE_FASTA}" || true

module unload picard

##########################################################################################################################
# 4. bamCoverage → Lorenz / gini
##########################################################################################################################
echo ">>> bamCoverage + lorenz"
module load deepTools
bamCoverage -b "${BAM}" -o "${BIGWIG}" --binSize 1000 -of bigwig -p 4 || true
module unload deepTools

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

if [ -f "${BIGWIG}" ]; then
    python "${SCRIPTS_DIR}/scripts/utils/lorenz.py" \
        "${BIGWIG}" -o "${LORENZ_CSV}" --gini-output "${GINI_TXT}" || true
fi

##########################################################################################################################
# 5. Splice-junction read count (primary alignments only, CIGAR contains N)
##########################################################################################################################
echo ">>> splice junction count"
SPLICE_READS=$(samtools view -F 0x900 "${BAM}" | awk '$6 ~ /N/' | wc -l)
echo "splice_junction_reads=${SPLICE_READS}"

##########################################################################################################################
# 6. Exonic enrichment
##########################################################################################################################
echo ">>> exonic enrichment"
python "${SCRIPTS_DIR}/scripts/CapGTA/calculate_exonic_enrichment.py" \
    "${BAM}" "${GTF_FILE}" "${REFERENCE_FASTA}" > "${EXONIC_TXT}" || true

##########################################################################################################################
# 7. Assemble one-row CSV for this cell
##########################################################################################################################
echo ">>> writing per-cell CSV: ${PER_CELL_CSV}"

# Total + aligned read counts (primary alignments; total = mapped + unmapped primary)
TOTAL_READS=$(samtools view -c -F 0x900 "${BAM}")
ALIGNED_READS=$(samtools view -c -F 0x904 "${BAM}")   # 0x900 + 0x004 (unmapped)

# Picard metrics files are tab-delimited. Force FS='\t' so a library name like
# "Unknown Library" (with a space) doesn't shift every data column by one.
PCT_DUP=$(awk -F'\t' '
    /^LIBRARY\t/ {
        for (i=1; i<=NF; i++) if ($i == "PERCENT_DUPLICATION") col=i
        getline
        if (col) print $col
        exit
    }' "${DUP_METRICS}" 2>/dev/null)
PCT_DUP=${PCT_DUP:-NA}

PCT_1X=$(awk -F'\t' '
    /^GENOME_TERRITORY\t/ {
        for (i=1; i<=NF; i++) if ($i == "PCT_1X") col=i
        getline
        if (col) print $col
        exit
    }' "${WGS_METRICS}" 2>/dev/null)
PCT_1X=${PCT_1X:-NA}

GINI=$(cat "${GINI_TXT}" 2>/dev/null | tr -d '\n')
GINI=${GINI:-NA}

EXONIC=$(cat "${EXONIC_TXT}" 2>/dev/null | tr -d '\n')
EXONIC=${EXONIC:-NA}

# Write header + row
{
    echo "barcode,total_reads,aligned_reads,pct_duplicates,pct_1x_coverage,gini,splice_junction_reads,exonic_enrichment"
    echo "${CELL},${TOTAL_READS},${ALIGNED_READS},${PCT_DUP},${PCT_1X},${GINI},${SPLICE_READS},${EXONIC}"
} > "${PER_CELL_CSV}"

echo "✓ ${PER_CELL_CSV}"
cat "${PER_CELL_CSV}"
echo "Done."
