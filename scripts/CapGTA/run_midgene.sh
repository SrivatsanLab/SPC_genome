#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/count_midgene.py
