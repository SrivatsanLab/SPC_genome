#!/bin/bash
# Download K562 replication timing (ENCODE UW Repli-seq, hg19) and the UCSC
# hg19 -> hg38 liftOver chain into data/reference/replication_timing/.
#
# The six bigWigs are UW's percentage-normalised signal for the G1b, S1, S2, S3,
# S4 and G2 sort fractions (ENCODE accessions ENCFF001GRX, GSF, GSJ, GSN, GSP,
# GSB). The ENCODE portal no longer serves these ENCODE2 files directly, so they
# come from UCSC's encodeDCC mirror and are saved under their ENCODE names;
# md5 sums are UCSC's and match the copies in
# /fh/fast/srivatsan_s/grp/SrivatsanLab/Sanjay/SPC_paper/res/repli-seq-analysis/data/raw_data.
#
# Next: scripts/utils/replication_timing_to_hg38.R
#
# Usage: bash scripts/utils/download_replication_timing.sh [out_dir]

set -euo pipefail

out_dir=${1:-data/reference/replication_timing}
mkdir -p "$out_dir"
ucsc=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq

# fraction  ENCODE accession  UCSC file  md5
files="
G1 ENCFF001GRX wgEncodeUwRepliSeqK562G1PctSignalRep1.bigWig 28c96061658d207976713e68af2d3b75
S1 ENCFF001GSF wgEncodeUwRepliSeqK562S1PctSignalRep1.bigWig 2851e2d9e504eb5c6658ffe33536d851
S2 ENCFF001GSJ wgEncodeUwRepliSeqK562S2PctSignalRep1.bigWig c58fd69eea44a556de8a94de4fad9056
S3 ENCFF001GSN wgEncodeUwRepliSeqK562S3PctSignalRep1.bigWig a29a665494b342915c57bcb1fefd92a8
S4 ENCFF001GSP wgEncodeUwRepliSeqK562S4PctSignalRep1.bigWig 18e369582e4f26cb4f62bb177231b077
G2 ENCFF001GSB wgEncodeUwRepliSeqK562G2PctSignalRep1.bigWig 36e5cbe48571171ca2175b8d1ded9b4f
"

echo "$files" | while read -r phase acc ucsc_file md5; do
    [[ -z "$phase" ]] && continue
    dest="$out_dir/K562_${phase}_${acc}.bigWig"
    if [[ ! -s "$dest" ]]; then
        echo "downloading $phase ($acc)"
        curl -sSfL --retry 3 -o "$dest" "$ucsc/$ucsc_file"
    fi
    echo "$md5  $dest" | md5sum -c --quiet - || { echo "md5 mismatch for $dest" >&2; exit 1; }
done

chain="$out_dir/hg19ToHg38.over.chain.gz"
if [[ ! -s "$chain" ]]; then
    echo "downloading hg19ToHg38 chain"
    curl -sSfL --retry 3 -o "$chain" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi
gzip -t "$chain"
echo "replication timing inputs in $out_dir"
