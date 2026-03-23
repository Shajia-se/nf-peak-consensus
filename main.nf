#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def peak_consensus_output = params.peak_consensus_output ?: 'peak_consensus_output'

def resolveBaseDir = { p ->
  def fp = file(p.toString())
  fp.isAbsolute() ? fp : file("${params.project_folder}/${p}")
}

process peak_consensus_per_condition {
  tag "${profile_name}:${condition}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir { "${params.project_folder}/${peak_consensus_output}/${profile_name}" }, mode: 'copy', overwrite: true

  input:
    tuple val(profile_name), val(condition), path(rep1_peak), path(rep2_peak)

  output:
    tuple val(profile_name), path("${condition}_consensus.bed"), emit: consensus_for_universe
    path("${condition}_rep1.sorted.bed")
    path("${condition}_rep2.sorted.bed")
    path("${condition}_overlap_pairs.bed")
    path("${condition}_consensus.bed")
    path("${condition}_consensus.summary.tsv")

  script:
  def recip = params.reciprocal_overlap ?: 0.7
  """
  set -euo pipefail

  sort -k1,1 -k2,2n ${rep1_peak} > ${condition}_rep1.sorted.bed
  sort -k1,1 -k2,2n ${rep2_peak} > ${condition}_rep2.sorted.bed

  bedtools intersect \
    -a ${condition}_rep1.sorted.bed \
    -b ${condition}_rep2.sorted.bed \
    -f ${recip} \
    -r \
    -wo > ${condition}_overlap_pairs.bed

  cut -f1-3 ${condition}_overlap_pairs.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > ${condition}_consensus.bed

  rep1_n=\$(wc -l < ${condition}_rep1.sorted.bed || echo 0)
  rep2_n=\$(wc -l < ${condition}_rep2.sorted.bed || echo 0)
  overlap_n=\$(wc -l < ${condition}_overlap_pairs.bed || echo 0)
  consensus_n=\$(wc -l < ${condition}_consensus.bed || echo 0)

  cat > ${condition}_consensus.summary.tsv << TSV
profile\tcondition\trep1_peaks\trep2_peaks\toverlap_pairs\tconsensus_peaks\treciprocal_overlap
${profile_name}\t${condition}\t\${rep1_n}\t\${rep2_n}\t\${overlap_n}\t\${consensus_n}\t${recip}
TSV
  """
}

process build_universe_peaks {
  tag "${profile_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir { "${params.project_folder}/${peak_consensus_output}/${profile_name}" }, mode: 'copy', overwrite: true

  input:
    tuple val(profile_name), path(consensus_beds)

  output:
    path("universe_peaks.bed")
    path("universe_peaks.summary.tsv")

  script:
  def bedArgs = consensus_beds.collect { "\"${it}\"" }.join(' ')
  def nInputs = consensus_beds.size()
  """
  set -euo pipefail

  cat ${bedArgs} \
    | awk 'BEGIN{OFS="\\t"} NF >= 3 {print \$1,\$2,\$3}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > universe_peaks.bed

  n_universe=\$(wc -l < universe_peaks.bed || echo 0)

  cat > universe_peaks.summary.tsv << TSV
profile\tn_consensus_inputs\tn_universe_peaks
${profile_name}\t${nInputs}\t\${n_universe}
TSV
  """
}

workflow {
  def peakExt = (params.consensus_peak_ext ?: 'narrowPeak').toString()
  def profiles = (params.consensus_profiles ?: params.consensus_macs3_profile ?: 'strict_q0.01')
    .toString()
    .split(',')
    *.trim()
    .findAll { it }
    .unique()
  def rows

  if (params.consensus_pairs_csv && file(params.consensus_pairs_csv).exists()) {
    rows = Channel
      .fromPath(params.consensus_pairs_csv, checkIfExists: true)
      .splitCsv(header: true)
      .flatMap { row ->
        assert row.condition && row.rep1_peaks && row.rep2_peaks : 'consensus_pairs_csv must contain: condition,rep1_peaks,rep2_peaks'
        def p1 = file(row.rep1_peaks.toString().trim())
        def p2 = file(row.rep2_peaks.toString().trim())
        assert p1.exists() : "rep1_peaks not found for ${row.condition}: ${p1}"
        assert p2.exists() : "rep2_peaks not found for ${row.condition}: ${p2}"
        def rowProfiles = row.profile_name ? [row.profile_name.toString().trim()] : profiles
        rowProfiles.collect { prof ->
          tuple(prof, row.condition.toString().trim(), p1, p2)
        }
      }
  } else if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    def header = null
    def records = []
    master.eachLine { line, n ->
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (n == 1) {
        header = cols
      } else {
        def rec = [:]
        header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
        records << rec
      }
    }

    assert header : "samples_master header not found: ${params.samples_master}"
    assert header.contains('sample_id') : "samples_master missing required column: sample_id"
    assert header.contains('condition') : "samples_master missing required column: condition"

    def isEnabled = { rec ->
      def v = rec.enabled?.toString()?.trim()?.toLowerCase()
      (v == null || v == '' || v == 'true')
    }
    def isControl = { rec ->
      rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
    }
    def isChip = { rec ->
      def lt = rec.library_type?.toString()?.trim()?.toLowerCase()
      !isControl(rec) && (lt == null || lt == '' || lt == 'chip')
    }

    def macsBase = resolveBaseDir(params.macs3_output)
    def pairs = []

    records
      .findAll { rec -> isEnabled(rec) && isChip(rec) }
      .groupBy { rec -> rec.condition?.toString()?.trim() ?: 'NA' }
      .each { cond, list ->
        def ordered = list.sort { a, b ->
          def ra = a.replicate?.toString()?.isInteger() ? a.replicate.toInteger() : Integer.MAX_VALUE
          def rb = b.replicate?.toString()?.isInteger() ? b.replicate.toInteger() : Integer.MAX_VALUE
          ra <=> rb
        }
        if (ordered.size() < 2) return
        assert ordered.size() == 2 : "Consensus currently expects exactly 2 enabled ChIP replicates per condition. Condition '${cond}' has ${ordered.size()}."

        def s1 = ordered[0].sample_id?.toString()?.trim()
        def s2 = ordered[1].sample_id?.toString()?.trim()
        if (!s1 || !s2) return

        profiles.each { prof ->
          def peakDir = file("${macsBase}/${prof}")
          assert peakDir.exists() : "MACS3 profile output not found: ${peakDir}"
          def p1 = file("${peakDir}/${s1}_peaks.${peakExt}")
          def p2 = file("${peakDir}/${s2}_peaks.${peakExt}")
          assert p1.exists() : "Peak file not found for ${s1} in profile ${prof}: ${p1}"
          assert p2.exists() : "Peak file not found for ${s2} in profile ${prof}: ${p2}"
          pairs << tuple(prof, cond, p1, p2)
        }
      }

    rows = Channel
      .fromList(pairs)
      .ifEmpty { exit 1, "ERROR: No valid consensus replicate pairs generated from samples_master: ${params.samples_master}" }
  } else {
    exit 1, 'ERROR: Provide --consensus_pairs_csv or --samples_master.'
  }

  def filtered = rows.filter { profile_name, condition, rep1_peak, rep2_peak ->
    !file("${params.project_folder}/${peak_consensus_output}/${profile_name}/${condition}_consensus.bed").exists()
  }

  def consensus_out = peak_consensus_per_condition(filtered)

  def existing_consensus = rows
    .map { profile_name, condition, rep1_peak, rep2_peak ->
      tuple(profile_name, file("${params.project_folder}/${peak_consensus_output}/${profile_name}/${condition}_consensus.bed"))
    }
    .filter { profile_name, bed -> bed.exists() }

  def universe_in = existing_consensus
    .mix(consensus_out.consensus_for_universe)
    .groupTuple()
    .map { profile_name, beds -> tuple(profile_name, beds) }

  build_universe_peaks(universe_in)
}
