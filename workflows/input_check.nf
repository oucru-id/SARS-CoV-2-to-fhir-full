#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def resolve_file(String path, String sample, String column) {
    if (!path) return null
    def f = file(path)
    if (!f.exists()) {
        f = file("${projectDir}/${path}")
    }
    if (!f.exists()) {
        error "Sample '${sample}': ${column} not found -> ${path}"
    }
    return f
}

def resolve_primer_bed(String scheme, String sample) {
    if (!scheme || scheme.toLowerCase() == 'none') return null
    def bed = file("${params.primer_scheme_dir}/${scheme}/primer.bed")
    if (!bed.exists()) {
        error "Sample '${sample}': primer scheme '${scheme}' has no primer.bed at ${bed}. " +
              "Run scripts/prepare_references.sh, or set primer_scheme to 'none'."
    }
    return bed
}

def resolve_primer_fasta(String scheme) {
    if (!scheme || scheme.toLowerCase() == 'none') return null
    def fa = file("${params.primer_scheme_dir}/${scheme}/primers.fasta")
    return fa.exists() ? fa : null
}

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    def seen = [] as Set

    rows = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def sample = row.sample?.trim()
            if (!sample) {
                error "Samplesheet has a row with an empty 'sample' column"
            }
            if (!seen.add(sample)) {
                error "Samplesheet contains duplicate sample id '${sample}'"
            }

            def platform = (row.platform ?: '').trim().toLowerCase()
            if (!(platform in ['illumina', 'nanopore'])) {
                error "Sample '${sample}': platform must be 'illumina' or 'nanopore', got '${platform}'"
            }

            def protocol = ((row.protocol ?: '').trim() ?: params.protocol).toLowerCase()
            if (!(protocol in ['amplicon', 'metagenomic'])) {
                error "Sample '${sample}': protocol must be 'amplicon' or 'metagenomic', got '${protocol}'"
            }

            def scheme = (row.primer_scheme ?: '').trim() ?: params.primer_scheme
            if (protocol == 'metagenomic') scheme = 'none'

            def lengths = params.scheme_read_lengths[scheme] ?: params.scheme_read_lengths['none']
            def min_len = (row.min_len ?: '').trim() ?: lengths[0]
            def max_len = (row.max_len ?: '').trim() ?: lengths[1]

            def meta = [
                id           : sample,
                platform     : platform,
                protocol     : protocol,
                scheme       : scheme,
                primer_bed   : resolve_primer_bed(scheme, sample),
                primer_fasta : resolve_primer_fasta(scheme),
                medaka_model : (row.medaka_model ?: '').trim() ?: params.medaka_model,
                min_len      : min_len as Integer,
                max_len      : max_len as Integer,
                single_end   : platform == 'nanopore'
            ]

            def r1 = resolve_file(row.fastq_1?.trim(), sample, 'fastq_1')
            def r2 = resolve_file(row.fastq_2?.trim(), sample, 'fastq_2')

            if (!r1) error "Sample '${sample}': fastq_1 is required"
            if (platform == 'illumina' && !r2) {
                error "Sample '${sample}': illumina samples require paired fastq_1 and fastq_2"
            }

            return [meta, r2 ? [r1, r2] : [r1]]
        }

    illumina = rows.filter { meta, reads -> meta.platform == 'illumina' }
    nanopore = rows.filter { meta, reads -> meta.platform == 'nanopore' }

    emit:
    illumina = illumina
    nanopore = nanopore
    all      = rows
}