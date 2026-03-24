# nf-peak-consensus

Nextflow DSL2 module for condition-level consensus peaks and exploratory universe peak sets using `bedtools intersect` + `bedtools merge`.

## What it does

For each configured MACS3 consensus profile and each condition:
- selects two replicates
- sorts both peak BED/narrowPeak files
- finds reciprocal overlaps with `bedtools intersect -f 0.7 -r`
- merges overlapping intervals into one consensus BED
- merges all condition-level consensus beds into one `universe_peaks.bed`

For `consensus_q0.05`, it also writes two explicit universe flavors:
- `consensus_first_universe_peaks.bed`
  - built from per-condition consensus peaks
  - default exploratory universe
- `union_first_universe_peaks.bed`
  - built directly from all sample MACS3 `q<0.05` peak files
  - broader supplementary universe

## Input modes

1. `--consensus_pairs_csv`
2. `--samples_master` (default recommended)

### Auto from `samples_master`

Required columns:
- `sample_id`
- `condition`

Optional columns used:
- `replicate`
- `library_type`
- `is_control`
- `enabled`

Default MACS3 profiles used:
- `strict_q0.01`
- `consensus_q0.05`

Default direct-union profile:
- `consensus_q0.05`

Resolved peak path:
- `${macs3_output}/${profile_name}/${sample_id}_peaks.${consensus_peak_ext}`

## Outputs

Under `${project_folder}/${peak_consensus_output}/${profile_name}`:

- `${condition}_rep1.sorted.bed`
- `${condition}_rep2.sorted.bed`
- `${condition}_overlap_pairs.bed`
- `${condition}_consensus.bed`
- `${condition}_consensus.summary.tsv`
- `universe_peaks.bed`
- `universe_peaks.summary.tsv`

Additional outputs under `${project_folder}/${peak_consensus_output}/consensus_q0.05`:
- `consensus_first_universe_peaks.bed`
- `consensus_first_universe_peaks.summary.tsv`
- `union_first_universe_peaks.bed`
- `union_first_universe_peaks.summary.tsv`

## Universe definitions

### `consensus_first_universe_peaks`
- input: all sample `MACS3 q<0.05` peaks
- within-condition: `bedtools intersect -f 0.7 -r`
- per-condition consensus: `bedtools merge -i -`
- across-condition universe: `cat + sort + bedtools merge -i -`

### `union_first_universe_peaks`
- input: all sample `MACS3 q<0.05` peaks
- direct union: `cat + sort + bedtools merge -i -`
- no extra distance-based stitching

## Run

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv
```

Explicit pairs:
```bash
nextflow run main.nf -profile hpc \
  --consensus_pairs_csv consensus_pairs.csv
```
