#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process HOSTILE_SHORT {
    tag "$meta.id"
    publishDir "${params.results_dir}/host_removal", mode: 'copy', pattern: '*.json'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clean/*.fastq.gz"), emit: reads
    path "${meta.id}.hostile.json",            emit: report

    script:
    """
    mkdir -p clean
    hostile clean \\
        --fastq1 ${reads[0]} \\
        --fastq2 ${reads[1]} \\
        --index ${params.hostile_index} \\
        --aligner bowtie2 \\
        --threads ${task.cpus} \\
        --force \\
        --output clean > ${meta.id}.hostile.json
    """

    stub:
    """
    mkdir -p clean
    cp ${reads[0]} clean/${meta.id}_1.clean.fastq.gz
    cp ${reads[1]} clean/${meta.id}_2.clean.fastq.gz
    echo '[]' > ${meta.id}.hostile.json
    """
}

process HOSTILE_LONG {
    tag "$meta.id"
    publishDir "${params.results_dir}/host_removal", mode: 'copy', pattern: '*.json'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("clean/*.fastq.gz"), emit: reads
    path "${meta.id}.hostile.json",            emit: report

    script:
    """
    mkdir -p clean
    hostile clean \\
        --fastq1 ${reads[0]} \\
        --index ${params.hostile_index} \\
        --aligner minimap2 \\
        --threads ${task.cpus} \\
        --force \\
        --output clean > ${meta.id}.hostile.json
    """

    stub:
    """
    mkdir -p clean
    cp ${reads[0]} clean/${meta.id}.clean.fastq.gz
    echo '[]' > ${meta.id}.hostile.json
    """
}

workflow HOST_REMOVAL {
    take:
    illumina_reads
    nanopore_reads

    main:
    if (params.skip_host_removal) {
        illumina_out = illumina_reads
        nanopore_out = nanopore_reads
        reports      = Channel.empty()
    }
    else {
        short_out = HOSTILE_SHORT(illumina_reads)
        long_out  = HOSTILE_LONG(nanopore_reads)

        illumina_out = short_out.reads.map { meta, r -> [meta, r.sort { it.name }] }
        nanopore_out = long_out.reads
        reports      = short_out.report.mix(long_out.report)
    }

    emit:
    illumina = illumina_out
    nanopore = nanopore_out
    reports  = reports
}
