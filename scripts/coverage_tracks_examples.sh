#!/bin/bash
###########################################################################################################################
# Example usage of coverage_tracks.py and coverage_tracks.sh
###########################################################################################################################

# Activate Python environment first
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# Example 1: Basic whole genome plot
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_whole_genome.png

# Example 2: Single chromosome plot (chr1)
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_chr1.png \
    -c chr1

# Example 3: Dark mode whole genome plot
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_dark_mode.png \
    --dark-mode

# Example 4: Custom color and line style
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_custom_color.png \
    --color '#ff1a5e' \
    --linewidth 1.5

# Example 5: Linear scale with custom y-limits
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_linear_scale.png \
    --no-log-scale \
    --ylim 0 50

# Example 6: Minimal plot without labels or dividers
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_minimal.png \
    --no-contig-labels \
    --no-dividers

# Example 7: Dark mode chr1 plot for presentation
python scripts/coverage_tracks.py \
    data/HSC_CellGrowth_coverage/sc_outputs/ACCTCGTAAGGATTCGGTAGACTCGCTGGATAAAGGCGACCTATA.bw \
    -o test_chr1_dark.png \
    -c chr1 \
    --dark-mode \
    --color '#80f15d' \
    --linewidth 2.0

# Example 8: Using the SLURM array wrapper for all cells
# Whole genome plots for all cells with default settings
sbatch --array=1-40 scripts/coverage_tracks.sh \
    bin/HSC_CellGrowth_coverage/sc_barcodes.txt \
    data/HSC_CellGrowth_coverage/sc_outputs \
    results/HSC_CellGrowth_coverage/coverage_tracks

# Example 9: Chr1 plots for all cells in dark mode
sbatch --array=1-40 scripts/coverage_tracks.sh \
    bin/HSC_CellGrowth_coverage/sc_barcodes.txt \
    data/HSC_CellGrowth_coverage/sc_outputs \
    results/HSC_CellGrowth_coverage/coverage_tracks_chr1_dark \
    -c chr1 \
    --dark-mode \
    --color '#80f15d'

# Example 10: Custom style for all cells
sbatch --array=1-40 scripts/coverage_tracks.sh \
    bin/HSC_CellGrowth_coverage/sc_barcodes.txt \
    data/HSC_CellGrowth_coverage/sc_outputs \
    results/HSC_CellGrowth_coverage/coverage_tracks_custom \
    --color '#ff1a5e' \
    --linewidth 1.5 \
    --gap-size 5000000 \
    --width 20 \
    --height 3
