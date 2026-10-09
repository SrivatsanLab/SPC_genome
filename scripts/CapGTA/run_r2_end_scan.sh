#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
for UDI in worm6_UDI1 worm6_UDI6 worm6_UDI11 worm6_UDI18; do
    echo "=== $UDI ==="
    python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/scan_r2_primer_end.py \
        --demux-dir /fh/working/srivatsan_s/Dustin/worm6-full-demux-r2trim/runs/${UDI}/results/${UDI}/cells \
        --n-cells 100
done
