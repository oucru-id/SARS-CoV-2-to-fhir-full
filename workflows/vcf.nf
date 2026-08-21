#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VCF_NORMALIZE {
    tag "$meta.id"

    input:
    tuple val(meta), path(vcf)
    path reference

    output:
    tuple val(meta), path("${meta.id}.norm.vcf.gz"), emit: vcf

    script:
    """
    if [[ "${vcf}" == *.gz ]]; then
        gunzip -c ${vcf} > input.vcf
    else
        cp ${vcf} input.vcf
    fi

    REF_CONTIG=\$(head -1 ${reference} | sed 's/^>//' | awk '{print \$1}')
    VCF_CONTIG=\$(grep -v '^#' input.vcf | head -1 | cut -f1)

    if [ -n "\$VCF_CONTIG" ] && [ "\$VCF_CONTIG" != "\$REF_CONTIG" ]; then
        echo "Renaming contig \$VCF_CONTIG -> \$REF_CONTIG"
        printf '%s\\t%s\\n' "\$VCF_CONTIG" "\$REF_CONTIG" > rename.txt
        bcftools annotate --rename-chrs rename.txt input.vcf -o renamed.vcf
    else
        cp input.vcf renamed.vcf
    fi

    bcftools norm -m -both -f ${reference} renamed.vcf -Oz -o ${meta.id}.norm.vcf.gz
    tabix -f -p vcf ${meta.id}.norm.vcf.gz
    """

    stub:
    "echo | gzip > ${meta.id}.norm.vcf.gz"
}

process VCF_FILTER {
    tag "$meta.id"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf

    script:
    """
    if bcftools view -h ${vcf} | grep -q '##INFO=<ID=DP'; then
        bcftools filter -i 'INFO/DP >= ${params.min_coverage}' ${vcf} -Oz -o filtered.vcf.gz
    elif bcftools view -h ${vcf} | grep -q '##FORMAT=<ID=DP'; then
        bcftools filter -i 'FORMAT/DP >= ${params.min_coverage}' ${vcf} -Oz -o filtered.vcf.gz
    else
        cp ${vcf} filtered.vcf.gz
    fi

    bcftools view -v snps,indels filtered.vcf.gz -Ou \\
        | bcftools sort -Oz -o ${meta.id}.vcf.gz
    tabix -f -p vcf ${meta.id}.vcf.gz
    """

    stub:
    "echo | gzip > ${meta.id}.vcf.gz"
}

workflow VCF_PROCESSING {
    take:
    vcfs
    reference

    main:
    normalized = VCF_NORMALIZE(vcfs, reference)
    filtered   = VCF_FILTER(normalized.vcf)

    emit:
    vcf = filtered.vcf
}
