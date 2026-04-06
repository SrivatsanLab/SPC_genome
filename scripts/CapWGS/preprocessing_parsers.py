#!/usr/bin/env python3
"""
Shared parsing functions for extracting preprocessing statistics from
demux, Trim Galore, and BWA-MEM log output.

Used by both save_chunk_stats.py (pipeline integration) and
parse_slurm_logs.py (retroactive parsing of old SLURM logs).
"""

import re


def parse_barcode_mapping(text):
    """Extract barcode mapping stats from demux output text.

    Looks for:
      - "proportion of reads with barcode mapped: 0.84..."
      - tqdm progress bar total (e.g., "0/29680") for input read count
      - "total_reads: X" and "mapped_reads: Y" if present (newer demux versions)

    Returns dict with keys: total_reads, mapped_reads, barcode_mapping_rate
    """
    result = {}

    # Try explicit counts first (newer demux output)
    m = re.search(r'total_reads:\s*(\d+)', text)
    if m:
        result['total_reads'] = int(m.group(1))
    m = re.search(r'mapped_reads:\s*(\d+)', text)
    if m:
        result['mapped_reads'] = int(m.group(1))

    # Mapping rate
    m = re.search(r'proportion of reads with barcode mapped:\s*([\d.]+)', text)
    if m:
        result['barcode_mapping_rate'] = float(m.group(1))

    # Fallback: extract total reads from tqdm progress bar
    if 'total_reads' not in result:
        # tqdm shows "0/29680" at start - find the largest denominator
        totals = re.findall(r'Processing reads:.*?\|.*?(\d+)/(\d+)', text)
        if totals:
            result['total_reads'] = int(totals[0][1])

    # Compute mapped_reads if we have total and rate but not mapped
    if 'mapped_reads' not in result and 'total_reads' in result and 'barcode_mapping_rate' in result:
        result['mapped_reads'] = int(round(result['total_reads'] * result['barcode_mapping_rate']))

    return result


def parse_trimming_sections(text):
    """Parse Trim Galore output for R1 and R2 trimming stats.

    Each read's section contains:
      - "Total reads processed: 24,995"
      - "Reads with adapters: 15,330 (61.3%)"
      - "Total basepairs processed: 3,874,225 bp"
      - "Quality-trimmed: 15,823 bp (0.4%)"
      - "Total written (filtered): 3,393,109 bp (87.6%)"
      - Adapter position histogram (length/count table)

    Returns list of dicts (one for R1, one for R2), each with keys:
      total_reads, reads_with_adapters, pct_reads_with_adapters,
      total_bp_processed, bp_quality_trimmed, bp_written,
      adapter_position_histogram (list of [length, count])
    """
    results = []

    # Split on "=== Summary ===" markers to find each read's section
    summary_sections = re.split(r'=== Summary ===', text)

    for section in summary_sections[1:]:  # skip text before first summary
        stats = {}

        m = re.search(r'Total reads processed:\s+([\d,]+)', section)
        if m:
            stats['total_reads'] = int(m.group(1).replace(',', ''))

        m = re.search(r'Reads with adapters:\s+([\d,]+)\s+\(([\d.]+)%\)', section)
        if m:
            stats['reads_with_adapters'] = int(m.group(1).replace(',', ''))
            stats['pct_reads_with_adapters'] = float(m.group(2))

        m = re.search(r'Total basepairs processed:\s+([\d,]+)\s+bp', section)
        if m:
            stats['total_bp_processed'] = int(m.group(1).replace(',', ''))

        m = re.search(r'Quality-trimmed:\s+([\d,]+)\s+bp', section)
        if m:
            stats['bp_quality_trimmed'] = int(m.group(1).replace(',', ''))

        m = re.search(r'Total written \(filtered\):\s+([\d,]+)\s+bp', section)
        if m:
            stats['bp_written'] = int(m.group(1).replace(',', ''))

        # Parse adapter position histogram
        histogram = []
        hist_match = re.search(r'Overview of removed sequences\nlength\tcount\texpect\tmax\.err\terror counts\n(.*?)(?:\n\n|\nRUN STATISTICS)', section, re.DOTALL)
        if hist_match:
            for line in hist_match.group(1).strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        length = int(parts[0])
                        count = int(parts[1])
                        histogram.append([length, count])
                    except ValueError:
                        continue
        stats['adapter_position_histogram'] = histogram

        if stats.get('total_reads') is not None:
            results.append(stats)

    return results


def parse_pairs_removed(text):
    """Extract number and percentage of pairs removed due to length cutoff.

    Looks for: "Number of sequence pairs removed ... : 688 (2.75%)"

    Returns dict with keys: pairs_removed, pct_pairs_removed
    """
    result = {}
    m = re.search(
        r'Number of sequence pairs removed.*?:\s*([\d,]+)\s+\(([\d.]+)%\)',
        text
    )
    if m:
        result['pairs_removed'] = int(m.group(1).replace(',', ''))
        result['pct_pairs_removed'] = float(m.group(2))
    return result


def parse_bwa_stats(text):
    """Extract alignment stats from BWA-MEM output.

    Parses:
      - [main] CMD line for reference index path
      - [main] Version line
      - [M::mem_process_seqs] for reads processed
      - [M::mem_pestat] for insert size stats (FR orientation)

    Returns dict with keys: bwa_version, reference_index_path,
      reads_processed, insert_size_mean, insert_size_stddev
    """
    result = {}

    m = re.search(r'\[main\] Version:\s*(.+)', text)
    if m:
        result['bwa_version'] = m.group(1).strip()

    # Extract reference from CMD line - it's the arg after -o <samfile>
    m = re.search(r'\[main\] CMD:\s*bwa mem\s+(.+)', text)
    if m:
        cmd = m.group(1)
        # Reference is the first non-option arg after -o <file>
        # Pattern: ... -o /path/to/output.sam /path/to/reference.fa ...
        ref_match = re.search(r'-o\s+\S+\s+(\S+)', cmd)
        if ref_match:
            result['reference_index_path'] = ref_match.group(1)

    # Total reads processed
    processed = re.findall(r'\[M::mem_process_seqs\] Processed (\d+) reads', text)
    if processed:
        result['reads_processed'] = sum(int(x) for x in processed)

    # Insert size stats (FR orientation is the standard one)
    m = re.search(
        r'analyzing insert size distribution for orientation FR\.\.\.\n'
        r'\[M::mem_pestat\] \(25, 50, 75\) percentile: \((\d+), (\d+), (\d+)\)\n'
        r'[^\n]+\n'  # boundaries line
        r'\[M::mem_pestat\] mean and std\.dev: \(([\d.]+), ([\d.]+)\)',
        text
    )
    if m:
        result['insert_size_median'] = int(m.group(2))
        result['insert_size_mean'] = float(m.group(4))
        result['insert_size_stddev'] = float(m.group(5))

    return result


def parse_flagstat(text):
    """Parse samtools flagstat output.

    Example input:
        48614 + 0 in total (QC-passed reads + QC-failed reads)
        0 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        45231 + 0 mapped (93.04% : N/A)
        48614 + 0 paired in sequencing
        24307 + 0 read1
        24307 + 0 read2
        38900 + 0 properly paired (80.02% : N/A)
        42000 + 0 with itself and mate mapped
        3231 + 0 singletons (6.65% : N/A)

    Returns dict with keys: total, secondary, supplementary, duplicates,
      mapped, mapped_pct, paired, properly_paired, properly_paired_pct,
      singletons, singletons_pct
    """
    result = {}

    def parse_line(pattern, text):
        m = re.search(pattern, text)
        if m:
            return m.groups()
        return None

    m = re.search(r'(\d+) \+ \d+ in total', text)
    if m:
        result['total'] = int(m.group(1))

    m = re.search(r'(\d+) \+ \d+ secondary', text)
    if m:
        result['secondary'] = int(m.group(1))

    m = re.search(r'(\d+) \+ \d+ supplementary', text)
    if m:
        result['supplementary'] = int(m.group(1))

    m = re.search(r'(\d+) \+ \d+ duplicates', text)
    if m:
        result['duplicates'] = int(m.group(1))

    m = re.search(r'(\d+) \+ \d+ mapped \(([\d.]+)%', text)
    if m:
        result['mapped'] = int(m.group(1))
        result['mapped_pct'] = float(m.group(2))

    m = re.search(r'(\d+) \+ \d+ paired in sequencing', text)
    if m:
        result['paired'] = int(m.group(1))

    m = re.search(r'(\d+) \+ \d+ properly paired \(([\d.]+)%', text)
    if m:
        result['properly_paired'] = int(m.group(1))
        result['properly_paired_pct'] = float(m.group(2))

    m = re.search(r'(\d+) \+ \d+ singletons \(([\d.]+)%', text)
    if m:
        result['singletons'] = int(m.group(1))
        result['singletons_pct'] = float(m.group(2))

    return result


def parse_full_log(text):
    """Parse a complete PP_array SLURM log or equivalent captured output.

    Returns dict with keys: demux, trimming (r1, r2, pairs_removed, pct_pairs_removed), alignment
    """
    demux = parse_barcode_mapping(text)
    trimming_sections = parse_trimming_sections(text)
    pairs = parse_pairs_removed(text)
    alignment = parse_bwa_stats(text)

    trimming = {}
    if len(trimming_sections) >= 1:
        trimming['r1'] = trimming_sections[0]
    if len(trimming_sections) >= 2:
        trimming['r2'] = trimming_sections[1]
    trimming.update(pairs)

    return {
        'demux': demux,
        'trimming': trimming,
        'alignment': alignment,
    }
