#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CONSENSUS {
    tag "$meta.id"
    publishDir "${params.results_dir}/consensus", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(bam)
    path reference

    output:
    tuple val(meta), path("${meta.id}.consensus.fa"), emit: fasta
    tuple val(meta), path("${meta.id}.consensus_stats.tsv"), emit: stats

    script:
    """
    samtools index -@ ${task.cpus} ${bam} 2>/dev/null || true

    bcftools filter -i 'INFO/DP >= ${params.min_coverage}' ${vcf} -Ou 2>/dev/null \\
        | bcftools norm -m -any -f ${reference} -Oz -o filtered.vcf.gz 2>/dev/null \\
        || bcftools norm -m -any -f ${reference} ${vcf} -Oz -o filtered.vcf.gz
    tabix -f -p vcf filtered.vcf.gz

    bedtools genomecov -bga -ibam ${bam} \\
        | awk -v m=${params.min_coverage} 'BEGIN{OFS="\\t"} \$4 < m {print \$1, \$2, \$3}' \\
        > lowcov.bed

    bcftools consensus \\
        --fasta-ref ${reference} \\
        --mask lowcov.bed \\
        --mask-with N \\
        filtered.vcf.gz \\
        > consensus.raw.fa

    awk -v s="${meta.id}" 'NR==1 {print ">" s; next} {print}' consensus.raw.fa \\
        > ${meta.id}.consensus.fa

    python3 ${projectDir}/scripts/consensus_stats.py \\
        --fasta ${meta.id}.consensus.fa \\
        --bam ${bam} \\
        --sample ${meta.id} \\
        --min-coverage ${params.min_coverage} \\
        --output ${meta.id}.consensus_stats.tsv
    """

    stub:
    """
    echo ">${meta.id}" > ${meta.id}.consensus.fa
    echo "N" >> ${meta.id}.consensus.fa
    touch ${meta.id}.consensus_stats.tsv
    """
}

workflow CONSENSUS_CALLING {
    take:
    vcf_bam      
    reference

    main:
    CONSENSUS(vcf_bam, reference)

    emit:
    fasta = CONSENSUS.out.fasta
    stats = CONSENSUS.out.stats
}