#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_NANOPORE {
    tag "$meta.id"
    publishDir "${params.results_dir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads} --outdir .
    """

    stub:
    "touch ${meta.id}_fastqc.zip"
}

process CHOPPER {
    tag "$meta.id"
    publishDir "${params.results_dir}/qc/chopper", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filt.fastq.gz"), emit: reads
    path "${meta.id}.chopper.log",                     emit: log

    script:
    """
    zcat ${reads} \\
        | chopper \\
            --quality ${params.nanopore_min_q} \\
            --minlength ${meta.min_len} \\
            --maxlength ${meta.max_len} \\
            --threads ${task.cpus} \\
            2> ${meta.id}.chopper.log \\
        | gzip > ${meta.id}.filt.fastq.gz
    """

    stub:
    """
    cp ${reads} ${meta.id}.filt.fastq.gz
    touch ${meta.id}.chopper.log
    """
}

process MINIMAP2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    minimap2 -a -x map-ont -t ${task.cpus} \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \\
        ${reference} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam
    """

    stub:
    "touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai"
}

process AMPLICONCLIP {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.clipped.bam"), path("${meta.id}.clipped.bam.bai"), emit: bam

    script:
    """
    samtools ampliconclip \\
        --both-ends \\
        --strand \\
        --filter-len 50 \\
        --hard-clip \\
        -b ${meta.primer_bed} \\
        -o ${meta.id}.clip.unsorted.bam \\
        ${bam}

    samtools sort -@ ${task.cpus} -o ${meta.id}.clipped.bam ${meta.id}.clip.unsorted.bam
    samtools index ${meta.id}.clipped.bam
    """

    stub:
    """
    cp ${bam} ${meta.id}.clipped.bam
    touch ${meta.id}.clipped.bam.bai
    """
}

process NANOPORE_BAM_STATS {
    tag "$meta.id"
    publishDir "${params.results_dir}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai), emit: bam
    path "${meta.id}.*", emit: stats

    script:
    """
    samtools stats ${bam}    > ${meta.id}.stats.txt
    samtools flagstat ${bam} > ${meta.id}.flagstat.txt
    samtools idxstats ${bam} > ${meta.id}.idxstats.txt
    mosdepth --by 200 --fast-mode ${meta.id} ${bam}
    """

    stub:
    """
    touch ${meta.id}.stats.txt ${meta.id}.flagstat.txt ${meta.id}.idxstats.txt
    """
}

process MEDAKA_INFERENCE {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    script:
    """
    cp -L ${reference} reference.fasta
    samtools faidx reference.fasta

    export CUDA_VISIBLE_DEVICES=""

    medaka inference \\
        --model ${meta.medaka_model} \\
        --threads ${task.cpus} \\
        --batch_size 100 \\
        ${bam} \\
        ${meta.id}.hdf

    medaka vcf ${meta.id}.hdf reference.fasta ${meta.id}.raw.vcf

    medaka tools annotate \\
        --dpsp \\
        ${meta.id}.raw.vcf \\
        reference.fasta \\
        ${bam} \\
        ${meta.id}.ann.vcf

    bcftools filter -i 'INFO/DP >= ${params.min_coverage}' ${meta.id}.ann.vcf -Ou \\
        | bcftools sort -Oz -o ${meta.id}.vcf.gz
    tabix -f -p vcf ${meta.id}.vcf.gz
    """

    stub:
    "echo | gzip > ${meta.id}.vcf.gz"
}

workflow NANOPORE {
    take:
    reads         
    reference

    main:
    qc       = FASTQC_NANOPORE(reads)
    filtered = CHOPPER(reads)
    aligned  = MINIMAP2(filtered.reads, reference)

    aligned.bam
        .branch { meta, bam, bai ->
            clip: meta.primer_bed != null
            keep: true
        }
        .set { by_protocol }

    clipped = AMPLICONCLIP(by_protocol.clip).bam.mix(by_protocol.keep)

    stats    = NANOPORE_BAM_STATS(clipped)
    variants = MEDAKA_INFERENCE(stats.bam, reference)

    emit:
    qc_report   = qc.zip
    chopper_log = filtered.log
    trimmed     = filtered.reads
    bam         = stats.bam
    align_stats = stats.stats
    vcf         = variants.vcf
}