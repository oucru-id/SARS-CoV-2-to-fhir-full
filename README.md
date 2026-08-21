# SARS-CoV-2 Genomics Mutation Analysis to FHIR Genomics Pipeline (CoV2toFHIR)

A platform-agnostic Nextflow pipeline for SARS-CoV-2 genomic analysis from raw sequencing data, producing HL7 FHIR R4 genomics bundles. Adapted from [nf-core/viralrecon 3.0.0](https://nf-co.re/viralrecon/3.0.0).

## Key Features

- **Multi-platform**: Illumina paired-end short reads, Oxford Nanopore (ONT) long reads, and pre-called VCF input.
- **Amplicon-aware**: Per-sample primer scheme handling.
- **Host Read Removal**: Human read depletion with `hostile`.
- **Variant Annotation**: snpEff against a `MN908947.3` database, giving gene, HGVS and Sequence Ontology effect per variant.
- **Lineage & Clade Assignment**: Nextclade Pango lineage, Nextstrain clade, WHO variant label and QC status.
- **Quality Control**: Aggregated to MultiQC.
- **FHIR Compliance**: HL7 FHIR R4 bundles with Variant, SARS-CoV-2 Panel, Lineage, Genome Quality Observations, and DiagnosticReport resources.
- **Clinical Integration**: Merges genomic results with patient, organization, and practitioner metadata.

## Installation

### Setup

```bash
git clone https://github.com/oucru-id/SARS-CoV-2-to-fhir-full.git
cd SARS-CoV-2-to-fhir-full

# Install Nextflow
curl -s https://get.nextflow.io | bash

# Verify
nextflow -v
```

### Dependencies

| Purpose | Tool |
|---|---|
| Read QC / trimming | fastqc, fastp, chopper |
| Host removal | hostile |
| Alignment | bowtie2, minimap2, samtools |
| Primer clipping | ivar, samtools ampliconclip |
| Variant calling | ivar, bcftools, medaka |
| Consensus | bcftools, bedtools |
| Annotation | snpEff, SnpSift |
| Lineage / clade | nextclade |
| Coverage / QC | mosdepth, picard, multiqc |
| Assembly | cutadapt, megahit, blast, abacas, quast |
| FHIR | python3 (pandas), java 17+ |

### Reference Preparation

One-off downloads the reference genome, primer schemes and Nextclade dataset, builds the snpEff database and Bowtie2 index:

```bash
bash scripts/prepare_references.sh
```

## Directory Structure

```
SARS-CoV-2-to-fhir-full
├── main.nf                             # Main workflow
├── nextflow.config                     # Configuration and parameters
├── workflows/
│   ├── input_check.nf                  # Samplesheet parsing and validation
│   ├── host_removal.nf                 # Human read depletion
│   ├── illumina.nf                     # Illumina sub-workflow
│   ├── nanopore.nf                     # Nanopore sub-workflow
│   ├── vcf.nf                          # VCF sub-workflow
│   ├── annotate.nf                     # snpEff annotation
│   ├── consensus.nf                    # Consensus genome calling
│   ├── nextclade.nf                    # Lineage and clade assignment
│   ├── assembly.nf                     # De novo assembly sub-workflow
│   ├── fhir.nf                         # FHIR variants generation
│   ├── validate_fhir.nf                # FHIR validation
│   ├── merge_clinical_data.nf          # Clinical metadata merge
│   ├── upload_fhir.nf                  # FHIR server upload
│   ├── report.nf                       # QC and sample report generation
│   └── utils.nf                        # Utility functions
├── scripts/
│   ├── prepare_references.sh           # Reference and database builder
│   ├── detect_primer_scheme.py         # Primer scheme detector
│   ├── ivar_variants_to_vcf.py         # iVar TSV to VCF converter
│   ├── consensus_stats.py              # Consensus quality metrics
│   ├── parse_nextclade.py              # Nextclade TSV to lineage JSON
│   ├── annotated_to_fhir.py            # VCF-to-FHIR converter
│   ├── clinical_metadata_parser.py     # Patient/org/practitioner parser
│   ├── generate_sample_report.py       # Per-sample text report
│   ├── merge_clinical_fhir.py          # FHIR genomics and clinical data merger
│   ├── upload_fhir.py                  # FHIR uploader
│   ├── get_access_token.py             # Standalone token fetcher
│   └── get_versions.py                 # Software version collector
├── data/
│   ├── NGS/                            # Input FASTQ files
│   ├── VCF/                            # Input VCF files
│   ├── samplesheet.csv                 # Sample definitions (required input)
│   ├── who_variant_classification.csv  # WHO variant tier lookup
│   ├── MN908947.3.fasta                # Reference genome (Wuhan-Hu-1)
│   ├── MN908947.3.gff                  # Reference annotation
│   ├── primer_schemes/                 # Amplicon primer BED and FASTA
│   ├── nextclade_db/                   # Nextclade SARS-CoV-2 dataset
│   ├── snpeff_db/                      # snpEff database
│   ├── patient_clinical_metadata.csv   # Patient metadata
│   ├── organization_metadata.csv       # Organization metadata
│   └── practitioner_metadata.csv       # Practitioner metadata
└── tools/
    └── fhir-validator.jar              # HL7 FHIR validator
```

### Samplesheet

`data/samplesheet.csv`:

```csv
sample,platform,protocol,primer_scheme,medaka_model,fastq_1,fastq_2,min_len,max_len
S1,illumina,amplicon,artic_v3,,data/NGS/S1_1.fastq.gz,data/NGS/S1_2.fastq.gz,,
S2,nanopore,amplicon,artic_v3,r941_prom_sup_variant_g507,data/NGS/S2.fastq.gz,,300,1200
```

Setting `protocol` to `metagenomic` forces `primer_scheme` to `none`.

### Illumina Reads
Place paired-end FASTQ files in `data/NGS/`, then reference them from the samplesheet.

### Nanopore Reads
Place single-end FASTQ files in `data/NGS/`, then reference them from the samplesheet.

### Pre-called VCFs
Place VCF files (`.vcf` or `.vcf.gz`) in `data/VCF/`. These are normalised to the reference contig, filtered and annotated alongside the read-derived samples.

### Determining the Primer Scheme

Clipping with the wrong scheme removes real sequence while leaving primer-derived bases in place, producing spurious variants at amplicon boundaries. Archives rarely record the scheme, so determine it empirically:

```bash
python3 scripts/detect_primer_scheme.py \
  --fastq data/NGS/SAMPLE.fastq.gz \
  --schemes data/primer_schemes \
  --reference data/MN908947.3.fasta \
  --platform nanopore
```

## Usage

### Get Access Token (FHIR Upload)

```bash
python scripts/get_access_token.py
```

### Basic Run

```bash
nextflow run main.nf
```

### Run with FHIR Upload

> Get the access token first before running with upload.

```bash
nextflow run main.nf \
  --fhir_server_url "https://<BASE_URL>/fhir"
```

### Common Options

```bash
nextflow run main.nf --skip_host_removal true

# Skip de novo assembly
nextflow run main.nf --skip_assembly true

# Use bcftools instead of iVar for Illumina variant calling
nextflow run main.nf --variant_caller bcftools

# Use an alternative samplesheet
nextflow run main.nf --input /path/to/samplesheet.csv
```

## WHO Variant Classification

| Classification | Criteria |
|---|---|
| VOC | Variant of concern |
| VOI | Variant of interest |
| VUM | Variant under monitoring |
| Previously designated | Formerly a WHO VOC or VOI |
| Not designated | Lineage carries no WHO variant designation |
| Designated variant, tier unknown | WHO label absent from the local classification table |

## Output Structure

```
results/
├── qc/
│   └── multiqc_report.html         # Aggregated QC report
├── host_removal/
│   └── *.hostile.json              # Host read depletion reports
├── alignment/
│   ├── *.bam                       # Primer-clipped alignments
│   └── *.mosdepth.*                # Coverage summaries
├── variants/
│   ├── *.annotated_variants.vcf.gz # snpEff-annotated variants
│   └── *.variants.tsv              # Flat variants table
├── consensus/
│   ├── *.consensus.fa              # Consensus genome
│   └── *.consensus_stats.tsv       # Length, %N, depth, coverage
├── nextclade/
│   └── *.lineage.json              # Per-sample lineage results
├── assembly/
│   ├── megahit/                    # De novo contigs
│   ├── blast/                      # Contig BLAST hits
│   └── quast/                      # Assembly assessment
├── fhir/
│   └── *.fhir.json                 # FHIR genomics bundles
├── fhir_merged/
│   └── *.merged.fhir.json          # FHIR bundles with clinical data
├── fhir_validated/
│   ├── *.validation.txt            # FHIR validation results
│   └── *.validation.summary.txt    # Error/warning counts by class
├── reports/
│   └── *.summary_report.txt        # Per-sample summary reports
├── runningstat/
│   ├── execution.html              # Nextflow execution report
│   ├── timeline.html               # Timeline report
│   └── dag.html                    # Workflow DAG
└── software_versions.yml           # Software version
```

## Support

[GitHub Issues](https://github.com/oucru-id/SARS-CoV-2-to-fhir-full/issues)
