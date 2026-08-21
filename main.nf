#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INPUT_CHECK }             from './workflows/input_check.nf'
include { HOST_REMOVAL }            from './workflows/host_removal.nf'
include { ILLUMINA }                from './workflows/illumina.nf'
include { NANOPORE }                from './workflows/nanopore.nf'
include { VCF_PROCESSING }          from './workflows/vcf.nf'
include { ANNOTATE }                from './workflows/annotate.nf'
include { CONSENSUS_CALLING }       from './workflows/consensus.nf'
include { LINEAGE }                 from './workflows/nextclade.nf'
include { ASSEMBLY }                from './workflows/assembly.nf'
include { GENERATE_REPORT }         from './workflows/report.nf'
include { GENERATE_SAMPLE_REPORTS } from './workflows/report.nf'
include { FHIR }                    from './workflows/fhir.nf'
include { MERGE_CLINICAL_DATA }     from './workflows/merge_clinical_data.nf'
include { VALIDATE }                from './workflows/validate_fhir.nf'
include { UPLOAD_FHIR }             from './workflows/upload_fhir.nf'
include { VERSIONS }                from './workflows/utils.nf'

workflow {

log.info """
    SARS-CoV-2 Genomic Surveillance and FHIR Pipeline (v${params.version})
    Reference : ${params.contig} (reported as ${params.fhir_reference_accession})
    Samplesheet: ${params.input}
    Developed by SPHERES OUCRU-ID
"""

    reference     = file(params.reference, checkIfExists: true)
    gff           = file(params.gff,       checkIfExists: true)
    bowtie2_index = file(params.bowtie2_index)

    inputs = INPUT_CHECK(params.input)
    cleaned = HOST_REMOVAL(inputs.illumina, inputs.nanopore)
    illumina_out = ILLUMINA(cleaned.illumina, reference, gff, bowtie2_index)
    nanopore_out = NANOPORE(cleaned.nanopore, reference)
    vcf_ch = Channel
        .fromPath("${params.vcf_dir}/*.vcf{,.gz}", checkIfExists: false)
        .map { f ->
            def id = f.name.replaceFirst(/\.vcf(\.gz)?$/, '')
            [[id: id, platform: 'vcf', protocol: 'metagenomic', scheme: 'none',
              primer_bed: null, primer_fasta: null, single_end: true], f]
        }
    vcf_out = VCF_PROCESSING(vcf_ch, reference)
    all_vcf = illumina_out.vcf
        .mix(nanopore_out.vcf)
        .mix(vcf_out.vcf)
        
    all_bam = illumina_out.bam
        .mix(nanopore_out.bam)

    annotated = ANNOTATE(all_vcf)
    vcf_bam = all_vcf
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(all_bam.map { meta, bam, bai -> [meta.id, bam] })
        .map { id, meta, vcf, bam -> [meta, vcf, bam] }

    consensus = CONSENSUS_CALLING(vcf_bam, reference)

    lineage = LINEAGE(consensus.fasta)

    ASSEMBLY(illumina_out.trimmed, consensus.fasta, reference, gff)

    per_sample = annotated.vcf
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(lineage.lineage.map { meta, j -> [meta.id, j] })
        .join(consensus.stats.map { meta, s -> [meta.id, s] })
        .map { id, meta, vcf, lin, stats -> [meta, vcf, lin, stats] }

    GENERATE_SAMPLE_REPORTS(per_sample)

    multiqc_files = illumina_out.qc_report.map { meta, f -> f }
        .mix(nanopore_out.qc_report.map { meta, f -> f })
        .mix(illumina_out.fastp_json)
        .mix(illumina_out.align_log)
        .mix(illumina_out.dup_metrics)
        .mix(illumina_out.align_stats)
        .mix(nanopore_out.align_stats)
        .mix(annotated.csv)
        .mix(cleaned.reports)
        .flatten()
        .collect()

    GENERATE_REPORT(multiqc_files)

    org_metadata_ch = Channel.fromPath(params.organization_metadata, checkIfExists: false).first()
    patient_metadata_ch = Channel.fromPath(params.patient_metadata, checkIfExists: false).first()
    practitioner_metadata_ch = Channel.fromPath(params.practitioner_metadata, checkIfExists: false).first()

    fhir_out = FHIR(per_sample, org_metadata_ch)

    merged_clinical_out = MERGE_CLINICAL_DATA(
        fhir_out.fhir_output,
        patient_metadata_ch,
        org_metadata_ch,
        practitioner_metadata_ch
    )

    validation_out = VALIDATE(merged_clinical_out.merged_fhir)

    UPLOAD_FHIR(validation_out.validated_fhir)

    VERSIONS()
}