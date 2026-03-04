# nf-peak-consensus

Nextflow DSL2 module for condition-level consensus peaks using `bedtools intersect` + `bedtools merge`.

## What it does

For each condition:
- selects two replicates
- sorts both peak BED/narrowPeak files
- finds reciprocal overlaps with `bedtools intersect -f 0.7 -r`
- merges overlapping intervals into one consensus BED

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

Default MACS3 profile used:
- `strict_q0.01`

Resolved peak path:
- `${macs3_output}/${consensus_macs3_profile}/${sample_id}_peaks.${consensus_peak_ext}`

## Outputs

Per condition:
- `${condition}_rep1.sorted.bed`
- `${condition}_rep2.sorted.bed`
- `${condition}_overlap_pairs.bed`
- `${condition}_consensus.bed`
- `${condition}_consensus.summary.tsv`

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
