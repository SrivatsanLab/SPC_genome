#!/usr/bin/env python3
"""Move UB:Z:<umi> from a FASTQ header comment into the read name suffix.

STAR takes the first whitespace-separated token as the read's query_name and
discards everything after it. To carry the UMI through alignment, we encode
it in the read name as `_UMI_<umi>` so pysam can parse it off the BAM later.

Input header format (from extract_primer_trimmed.py):
    @LH00740:234:... CB:Z:<cell> UB:Z:<umi>

Output header format:
    @LH00740:234:..._UMI_<umi> CB:Z:<cell>

The CB:Z: tag is left in the comment half — STAR will discard it in the BAM
but we don't need it per-read anyway (one cell = one BAM = one cell barcode).

Usage:
    encode_umi_in_readname.py <in.fastq.gz> <out.fastq.gz>
"""
import gzip
import re
import sys


UMI_PAT = re.compile(r'\bUB:Z:([ACGTN]+)')


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    in_path, out_path = sys.argv[1], sys.argv[2]

    opener_in = gzip.open if in_path.endswith('.gz') else open
    opener_out = gzip.open if out_path.endswith('.gz') else open

    n_total = 0
    n_umi = 0
    with opener_in(in_path, 'rt') as fin, opener_out(out_path, 'wt') as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                n_total += 1
                if line.startswith('@'):
                    m = UMI_PAT.search(line)
                    if m:
                        n_umi += 1
                        umi = m.group(1)
                        line_stripped = line.rstrip('\n')
                        # Split into name + comment
                        parts = line_stripped.split(None, 1)
                        name = parts[0]
                        comment = parts[1] if len(parts) == 2 else ''
                        # Remove UB:Z:<umi> from the comment to avoid duplication
                        comment = UMI_PAT.sub('', comment).strip()
                        # Rewrite header
                        new_header = f'{name}_UMI_{umi}'
                        if comment:
                            new_header += f' {comment}'
                        line = new_header + '\n'
            fout.write(line)

    print(f'  {in_path} -> {out_path}: {n_umi}/{n_total} records had UB:Z: tag',
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
