#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
# Pool from multiple UDIs for a representative sample
mkdir -p /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/r2_primer_scan
# scan 100 cells per UDI from UDI1,UDI6,UDI11,UDI18 = 400 cells pooled
for UDI in worm6_UDI1 worm6_UDI6 worm6_UDI11 worm6_UDI18; do
    echo "=== $UDI ==="
    python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/scan_r2_primer.py \
        --demux-dir /fh/working/srivatsan_s/Dustin/worm6-full-demux-r2trim/runs/${UDI}/results/${UDI}/cells \
        --n-cells 100 \
        --out-dir /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/r2_primer_scan/${UDI}
done
