#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_ILLUMINA {
    tag "$meta.id"
    publishDir "${params.results_dir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads.join(' ')} --outdir .
    """

    stub:
    "touch ${meta.id}_fastqc.zip"
}

process FASTP {
    tag "$meta.id"
    publishDir "${params.results_dir}/qc/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json",             emit: json
    path "${meta.id}.fastp.html",             emit: html

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_1.trim.fastq.gz \\
        --out2 ${meta.id}_2.trim.fastq.gz \\
        --cut_front --cut_tail --cut_mean_quality ${params.min_base_quality} \\
        --length_required ${params.illumina_min_len} \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        2> ${meta.id}.fastp.log
    """

    stub:
    """
    cp ${reads[0]} ${meta.id}_1.trim.fastq.gz
    cp ${reads[1]} ${meta.id}_2.trim.fastq.gz
    touch ${meta.id}.fastp.json ${meta.id}.fastp.html
    """
}

process BOWTIE2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path index_dir

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "${meta.id}.bowtie2.log", emit: log

    script:
    """
    bowtie2 \\
        --threads ${task.cpus} \\
        -x ${index_dir}/${params.contig} \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --rg-id ${meta.id} --rg "SM:${meta.id}" --rg "PL:ILLUMINA" \\
        2> ${meta.id}.bowtie2.log \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai ${meta.id}.bowtie2.log
    """
}

process IVAR_TRIM {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.trimmed.bam"), path("${meta.id}.trimmed.bam.bai"), emit: bam
    path "${meta.id}.ivar_trim.log", emit: log

    script:
    """
    ivar trim \\
        -i ${bam} \\
        -b ${meta.primer_bed} \\
        -p ${meta.id}.trim.unsorted \\
        -m 30 -q ${params.min_base_quality} -e \\
        > ${meta.id}.ivar_trim.log 2>&1

    samtools sort -@ ${task.cpus} -o ${meta.id}.trimmed.bam ${meta.id}.trim.unsorted.bam
    samtools index ${meta.id}.trimmed.bam
    """

    stub:
    """
    cp ${bam} ${meta.id}.trimmed.bam
    touch ${meta.id}.trimmed.bam.bai ${meta.id}.ivar_trim.log
    """
}

process MARK_DUPLICATES {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bam.bai"), emit: bam
    path "${meta.id}.markdup.metrics.txt", emit: metrics

    script:
    """
    PicardCommandLine MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.id}.md.bam \\
        --METRICS_FILE ${meta.id}.markdup.metrics.txt \\
        --REMOVE_DUPLICATES false \\
        --VALIDATION_STRINGENCY LENIENT
    samtools index ${meta.id}.md.bam
    """

    stub:
    """
    cp ${bam} ${meta.id}.md.bam
    touch ${meta.id}.md.bam.bai ${meta.id}.markdup.metrics.txt
    """
}

process BAM_STATS {
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

process IVAR_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path gff

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    script:
    """
    samtools mpileup \\
        --ignore-overlaps \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 0 \\
        --min-BQ ${params.min_base_quality} \\
        --reference ${reference} \\
        ${bam} \\
        | ivar variants \\
            -t ${params.min_allele_freq} \\
            -q ${params.min_base_quality} \\
            -m ${params.min_coverage} \\
            -r ${reference} \\
            -p ${meta.id}

    python3 ${projectDir}/scripts/ivar_variants_to_vcf.py \\
        --input ${meta.id}.tsv \\
        --output ${meta.id}.vcf \\
        --contig ${params.contig} \\
        --sample ${meta.id}

    bgzip -f ${meta.id}.vcf
    tabix -f -p vcf ${meta.id}.vcf.gz
    """

    stub:
    "echo | gzip > ${meta.id}.vcf.gz"
}

process BCFTOOLS_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    script:
    """
    bcftools mpileup \\
        --ignore-RG \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 0 \\
        --min-BQ ${params.min_base_quality} \\
        --annotate FORMAT/AD,FORMAT/DP,INFO/AD \\
        --fasta-ref ${reference} \\
        ${bam} \\
        | bcftools call --ploidy 1 --keep-alts --multiallelic-caller --variants-only -Ou \\
        | bcftools norm -m -any -f ${reference} -Ou \\
        | bcftools filter -i 'INFO/DP >= ${params.min_coverage}' -Oz -o ${meta.id}.vcf.gz

    tabix -f -p vcf ${meta.id}.vcf.gz
    """

    stub:
    "echo | gzip > ${meta.id}.vcf.gz"
}

workflow ILLUMINA {
    take:
    reads          
    reference
    gff
    bowtie2_index

    main:
    qc      = FASTQC_ILLUMINA(reads)
    trimmed = FASTP(reads)
    aligned = BOWTIE2(trimmed.reads, bowtie2_index)

    aligned.bam
        .branch { meta, bam, bai ->
            clip: meta.primer_bed != null
            keep: true
        }
        .set { by_protocol }

    trimmed_bam = IVAR_TRIM(by_protocol.clip).bam.mix(by_protocol.keep)

    if (params.skip_markduplicates) {
        final_bam   = trimmed_bam
        dup_metrics = Channel.empty()
    }
    else {
        md          = MARK_DUPLICATES(trimmed_bam)
        final_bam   = md.bam
        dup_metrics = md.metrics
    }

    stats = BAM_STATS(final_bam)

    variants = params.variant_caller == 'bcftools'
        ? BCFTOOLS_VARIANTS(stats.bam, reference).vcf
        : IVAR_VARIANTS(stats.bam, reference, gff).vcf

    emit:
    qc_report   = qc.zip
    fastp_json  = trimmed.json
    trimmed     = trimmed.reads
    bam         = stats.bam
    align_stats = stats.stats
    align_log   = aligned.log
    dup_metrics = dup_metrics
    vcf         = variants
}
