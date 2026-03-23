# nf-peak-consensus

Nextflow DSL2 module for condition-level consensus peaks using `bedtools intersect` + `bedtools merge`.

## What it does

For each configured MACS3 consensus profile and each condition:
- selects two replicates
- sorts both peak BED/narrowPeak files
- finds reciprocal overlaps with `bedtools intersect -f 0.7 -r`
- merges overlapping intervals into one consensus BED
- merges all condition-level consensus beds into one `universe_peaks.bed`

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
