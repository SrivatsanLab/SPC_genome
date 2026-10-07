# K562 sc_PolE copy number: review of the GC-corrected AneuFinder run

Date: 2026-10-02. Run: `scripts/utils/run_aneufinder_K562_sc_PolE_gc.sh` (1000 sc_PolE_novaseq
cells, 1 Mb bins, edivisive, GC correction, ENCODE blacklist, MAPQ >= 10, duplicates removed).
Output: `results/K562_tree/aneufinder_gc/output/`. Analysis: `paper_figures/scripts/K562_offpeak_ploidy.R`
-> `paper_figures/output/K562_offpeak_ploidy/`.

## Run health

- Completed in 133 min, no errors; one warning (Bhattacharyya computed on 2- vs 3-somy for one cell).
- `paper_figures/data/Anneufinder/result.csv` is byte-identical to the run output.
- AneuFinder QC, medians: 1.74 M reads/cell, spikiness 0.15, Bhattacharyya 0.77, 107 segments.

## Consensus profile (median of on-peak cells) matches K562

- Modal autosomal copy number 3 in 910/1000 cells.
- 9p21 (CDKN2A) homozygous deletion: CN0 in 76% of cells, CN1 in the rest.
- BCR (22q11) CN9 in 63% of cells; ABL1 (9q34) CN 3-8; 13q31 (chr13 90-92 Mb) up to CN45 in
  616 cells - consistent with the K562 BCR-ABL1 amplicon (worth confirming against the karyotype
  literature).
- CN2: 9p, 13q, 14q, 17p, 18q, 20p, distal 22q, X. chrY absent.
- Pearson r with the K562 reference profile (`K562_ref_ploidy_hg38.bed`;
  `paper_figures/scripts/K562_cnv_reference_correlation.R`, notebook cells 134-148), after
  dropping 212 unmappable bins (>= 50% N gap or ENCODE blacklist; no reads after blacklisting,
  CN0 in every cell) and chrY: all 1000 cells mean 0.72 (s.d. 0.10, median 0.75); 890 on-peak
  cells 0.74 (s.d. 0.05). Bin-shuffle null (1000 permutations) max 0.06, empirical p < 0.001.
  Without the mask the notebook method gives 0.54 (s.d. 0.07), because those bins count as
  mismatches against the reference (or against the 3 filled in where it has no segment).
  The old uncorrected run gave 0.62 (s.d. 0.09).

## Ploidy

- `ploidy` in `cell_ploidy.csv` is the mean over all 3078 bins, including chrY and 184 bins that
  are CN0 in most cells (centromeres, gaps, acrocentric arms). On-peak median 2.71 by that measure,
  2.93 over informative bins.
- 890 cells (89.0%) at [2.5, 3.0). The other 110 (65 below, 45 above) spread continuously from
  1.0 to 5.6 with no peaks at whole-genome multiples. The 4.3N / 5.4N peaks of the uncorrected run
  were GC-bias scaling errors (`K562_gc_toggle_ploidy.R`).
- Off-peak fraction is even across populations (P_0-P_5: 9-13%).

## Are the 110 off-peak cells artefacts or real?

Tested on GC-corrected read depth before copy-number calling (definitions in the script header):

| Class | below 2.5 | above 3.0 | total |
|---|---|---|---|
| Artefact: noisy or depth not on whole-number levels | 16 | 23 | 39 |
| Artefact: depth matches the consensus, ploidy only mis-scaled | 1 | 4 | 5 |
| Real-looking: large deviations on whole-number levels | 48 | 18 | 66 |

- The depth deviations are genuine signal, not a calling error: CN calls track depth at r = 0.91
  in off-peak cells (0.84 on-peak).
- Not residual GC bias (residual-GC r ~ 0) and not a single shared pattern such as S-phase
  replication timing (deviations are cell-specific).
- Spikiness does not explain them (off-peak median 0.158 vs 0.148); 99-107 of the 110 pass
  spikiness cut-offs between the on-peak 95th percentile and 0.30.
- High-ploidy cells mostly fail the whole-number test (lattice error 0.17 vs 0.105 on-peak);
  low-ploidy cells pass it (0.10).
- Interpretation: ~44 artefacts; ~66 cells (6.6%) with large, cell-specific copy-number changes,
  spread evenly across populations. Depth alone cannot separate one-off aneuploidy
  (missegregation, cells in crisis) from capsule-level amplification bias.
- The per-cell SBS tree is not used as evidence here (not considered reliable); sharing is
  assessed by population instead.

## Other issues

- 79 on-peak cells have Bhattacharyya < 0.5; disagreement with the consensus tracks spikiness
  (rho 0.62) and Bhattacharyya (rho -0.43), not depth. Filter before the CNV event tree.
- Residual artefact bins: chr2 90-92 Mb (pericentromeric, CN up to 87 in 7 cells); chr7 106-127 Mb
  CN ~15 in 5 high-ploidy cells from four populations. The repeat mask in `K562_cnv_event_tree.R`
  is off by default.
- BCR bin: 22% of cells at CN2 vs 63% at CN9 - a ~1 Mb amplicon at 1 Mb bins; not biology.
- Excluding all 110 off-peak cells from CNV event calling remains correct.

## Population-specific events

Populations P_0-P_5 are successive splits of the culture (P_0 first). Scan:
`paper_figures/scripts/K562_cnv_population_scan.R` on the 1 Mb, 500 kb and 250 kb runs; on-peak cells
passing QC (Bhattacharyya >= on-peak 10th percentile, spikiness <= 95th; percentiles because both
metrics shift with bin size - Bhattacharyya median 0.77 / 0.64 / 0.54). Hit = run of >= 3 bins
changed in >= 15% of one population and <= 5% of the rest.

- Population median profiles are identical apart from single pericentromeric or edge bins.
- Scan hits: 1 Mb chr7q11 gain (67-72 Mb; P_4); 500 kb chr9q22 loss (91.5-93.5 Mb; P_4) and chr4q28
  loss (131.5-136.5 Mb; P_5); 250 kb none.
- Residual per-cell GC bias confounds these calls. Tested each on GC-corrected 250 kb depth (region /
  flanks), regressing on per-cell GC bias (`K562_gc_toggle_ploidy.R`, quartic fit) and comparing
  populations on the residual:

| Region | GC (region / flanks) | depth ~ GC bias R^2 | population mean shift (perm. p) | cells beyond one-copy-level residual, in focal pop |
|---|---|---|---|---|
| chr7q11 gain | 0.43 / 0.39 | 0.75 | P_4 +0.017, ~7% of a copy (0.004) | 24 of 64 in P_4 (exp. 13.5; p 1e-3) |
| chr9q22 loss | 0.45 / 0.40 | 0.38 | P_4 -0.009 (0.31) | 22 of 53 in P_4 (exp. 11.2; p 4e-4) |
| chr4q28 loss | 0.35 / 0.37-0.38 | 0.19 | P_5 -0.011 (0.20) | 12 of 55 in P_5 (P_4 15; p 0.16) |

- chr7q11: whether a cell reads high is mostly GC bias (1% of the lowest GC-bias quintile vs 99% of
  the top two); 1 Mb calls caught only the extreme cells and finer bins call ~200 cells across all
  populations. Not a population-wide gain.
- chr4q28: not population-specific once GC bias is accounted for; ~55 cells with loss-level depth in
  every population.
- chr7q11 gain and chr9q22 loss co-occur in 21 cells (3.8 expected), 11 in P_4 (others P_3 4, P_5 4,
  P_2 2). Both regions are GC-rich but change in opposite directions, so not simple GC bias - but 7 of
  the 15 such cells that pass QC fall in the P_4 similarity group below, whose shared signature looks
  technical. Not a convincing subclone.
- Conclusion: no copy-number change marks a single population. Population tests on these calls need
  adjusting for per-cell GC bias.
- Reference correlation (notebook method, masked) is stable across bin sizes: mean r 0.72 / 0.71 /
  0.67 at 1 Mb / 500 kb / 250 kb (all cells).

## Small clones: do highly similar cells share a population?

`paper_figures/scripts/K562_cnv_cell_similarity.R` (1 Mb, 777 on-peak QC-passing cells, 2830 unmasked
bins) -> `paper_figures/output/K562_cnv_cell_similarity/`. Similarity = Pearson r of GC-adjusted depth
profiles (each bin's trend with per-cell GC bias removed). Random pairs share a population 18.5% of
the time.

- Similarity on copy-number deviations tracks GC-bias similarity (rho 0.40) and its top pairs are not
  population-concordant (8-22%). Not used.
- GC-adjusted depth (rho with GC-bias similarity -0.01): top 25 / 250 / 5000 pairs 36% / 29% / 23%
  same population; nearest neighbour 28% (label shuffle p < 0.0005).
- Barcode segment 4 is the round-1 well and encodes population (99% purity; 4 wells per population,
  11 for P_2). Sharing a later-round well (segments 1-3) adds nothing. Sharing a round-1 well does:
  top-0.1% pair rate 0.21% same population + same well, 0.13% same population + different well,
  0.09% different population.
- Shuffling whole round-1 wells between populations: population still adds similarity beyond the well
  on average (p = 0.0005), but for the most similar pairs only 1.5x (p = 0.047).
- Groups of cells joined by top-0.1% pairs (r >= 0.48): 133 cells in 20 groups.
  - Technical-looking: 59 cells from all populations with low GC bias and scattered 3-5 Mb deviations;
    15 cells all P_4 (all four P_4 wells; higher spikiness, 0.195 vs 0.141) plus a 3-cell P_4 group,
    sharing losses at 1q12 heterochromatin, chromosome ends and HLA; 10 cells with GC-rich 1p/16p
    losses; a 5-cell group with dozens of small deviations.
  - Clone-like, with large shared events and normal QC: 6 cells, 1p (2-54 Mb) + 18q (22-43 Mb) loss
    (P_0, P_2 x2, P_4 x2, P_5); 5 cells, 9p21 retained at CN1 + 17p loss (P_2 x4, P_4); 4 cells, 9p21
    retained (P_2 x2, P_4, P_5); 3 cells, 2q loss 95-232 Mb (P_2 x2, P_5); 3 cells, 12q loss
    37-133 Mb + 21q loss (P_0, P_2, P_3).
  - The 10 isolated pairs: 20% same population (random 18.5%).
- Interpretation: the clear small clones mostly span several populations, as expected for clones
  that arose before the culture was split; only the 9p21/17p clone is concentrated in one population
  (P_2). The strongest population-concordant group (P_4) carries a technical-looking signature, and
  some population concordance comes from round-1 well effects.

## Follow-up

Replication (cell cycle) turned out to drive most of the per-cell "GC bias" and noise structure in
the copy-number calls. Fragment-level GC correction, repeat/mappability masking, cell-cycle phase
inference, replication regression, corrected copy number and the tree/population tests are written
up in `docs/K562_cell_cycle_cnv.md`. Where that conflicts with the population and subclone
conclusions above (made on uncorrected calls), the later document supersedes them.

