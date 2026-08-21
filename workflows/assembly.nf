#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CUTADAPT {
    tag "$meta.id"
    publishDir "${params.results_dir}/assembly/cutadapt", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.primerless.fastq.gz"), emit: reads
    path "${meta.id}.cutadapt.log",                 emit: log

    script:
    """
    cutadapt \\
        -g file:${meta.primer_fasta} \\
        -G file:${meta.primer_fasta} \\
        --overlap 5 \\
        --minimum-length 30 \\
        --error-rate 0.1 \\
        --times 2 \\
        --cores ${task.cpus} \\
        -o ${meta.id}_1.primerless.fastq.gz \\
        -p ${meta.id}_2.primerless.fastq.gz \\
        ${reads[0]} ${reads[1]} \\
        > ${meta.id}.cutadapt.log
    """

    stub:
    """
    cp ${reads[0]} ${meta.id}_1.primerless.fastq.gz
    cp ${reads[1]} ${meta.id}_2.primerless.fastq.gz
    touch ${meta.id}.cutadapt.log
    """
}

process MEGAHIT {
    tag "$meta.id"
    publishDir "${params.results_dir}/assembly/megahit", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs
    path "${meta.id}.megahit.log",                  emit: log

    script:
    """
    megahit \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --presets meta-sensitive \\
        --num-cpu-threads ${task.cpus} \\
        --memory 0.4 \\
        --out-dir megahit_out \\
        --out-prefix ${meta.id} \\
        > ${meta.id}.megahit.log 2>&1

    cp megahit_out/${meta.id}.contigs.fa ${meta.id}.contigs.fa
    """

    stub:
    """
    echo ">contig_1" > ${meta.id}.contigs.fa
    echo "ACGT" >> ${meta.id}.contigs.fa
    touch ${meta.id}.megahit.log
    """
}

process BLAST_CONTIGS {
    tag "$meta.id"
    publishDir "${params.results_dir}/assembly/blast", mode: 'copy'

    input:
    tuple val(meta), path(contigs)
    path reference

    output:
    tuple val(meta), path("${meta.id}.blast.tsv"), emit: hits

    script:
    """
    makeblastdb -in ${reference} -dbtype nucl -out ref_db > /dev/null

    printf 'qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n' \\
        > ${meta.id}.blast.tsv

    blastn \\
        -query ${contigs} \\
        -db ref_db \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \\
        -num_threads ${task.cpus} \\
        -evalue 1e-10 \\
        >> ${meta.id}.blast.tsv
    """

    stub:
    "touch ${meta.id}.blast.tsv"
}

process ABACAS {
    tag "$meta.id"
    publishDir "${params.results_dir}/assembly/abacas", mode: 'copy'

    input:
    tuple val(meta), path(contigs)
    path reference

    output:
    path "${meta.id}.abacas.*", emit: results

    script:
    """
    abacas -r ${reference} -q ${contigs} -p nucmer -m -b -o ${meta.id}.abacas \\
        > ${meta.id}.abacas.log 2>&1 || {
        echo "abacas failed for ${meta.id}; see log" >> ${meta.id}.abacas.log
    }
    for f in ${contigs}_*; do
        [ -e "\$f" ] && mv "\$f" "${meta.id}.abacas.\${f##*_}" || true
    done
    touch ${meta.id}.abacas.log
    """

    stub:
    "touch ${meta.id}.abacas.log"
}

process QUAST {
    publishDir "${params.results_dir}/assembly/quast", mode: 'copy'

    input:
    path assemblies
    path reference
    path gff

    output:
    path "quast_results/*", emit: results
    path "quast_results/report.tsv", emit: report

    script:
    """
    quast.py \\
        -r ${reference} \\
        --features ${gff} \\
        --threads ${task.cpus} \\
        --min-contig 200 \\
        -o quast_results \\
        ${assemblies}
    """

    stub:
    """
    mkdir -p quast_results
    touch quast_results/report.tsv
    """
}

workflow ASSEMBLY {
    take:
    reads        
    consensus   
    reference
    gff

    main:
    contigs   = Channel.empty()
    blast     = Channel.empty()
    quast_rep = Channel.empty()

    if (!params.skip_assembly) {
        reads
            .branch { meta, r ->
                trim: meta.primer_fasta != null
                keep: true
            }
            .set { by_scheme }

        prepped = CUTADAPT(by_scheme.trim).reads.mix(by_scheme.keep)

        asm     = MEGAHIT(prepped)
        contigs = asm.contigs
        blast   = BLAST_CONTIGS(contigs, reference).hits

        if (!params.skip_abacas) {
            ABACAS(contigs, reference)
        }
    }

    if (!params.skip_quast) {
        to_assess = contigs.map { meta, f -> f }
            .mix(consensus.map { meta, f -> f })
            .collect()
        quast_rep = QUAST(to_assess, reference, gff).report
    }

    emit:
    contigs = contigs
    blast   = blast
    quast   = quast_rep
}
