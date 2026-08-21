#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys
import uuid
from collections import defaultdict
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clinical_metadata_parser import load_organization_metadata  # noqa: E402

DEFAULT_ORG_ID = '100007732'

LOINC = 'http://loinc.org'
SNOMED = 'http://snomed.info/sct'
SO_SYSTEM = 'http://www.sequenceontology.org'
HGVS_SYSTEM = 'http://varnomen.hgvs.org'
NCBI_GENE = 'https://www.ncbi.nlm.nih.gov/gene'
NCBI_REFSEQ = 'http://www.ncbi.nlm.nih.gov/refseq'
UCUM = 'http://unitsofmeasure.org'
KEMKES_SP = 'http://terminology.kemkes.go.id/sp'
OBS_CATEGORY = 'http://terminology.hl7.org/CodeSystem/observation-category'
V2_0074 = 'http://terminology.hl7.org/CodeSystem/v2-0074'
PANGO_SYSTEM = 'https://cov-lineages.org/lineage_list.html'
NEXTSTRAIN_CLADE = 'https://clades.nextstrain.org'
WHO_VARIANT_CLASS_SYSTEM = 'https://www.who.int/activities/tracking-SARS-CoV-2-variants'

WHO_CLASS_DISPLAY = {
    'VOC': 'Variant of concern',
    'VOI': 'Variant of interest',
    'VUM': 'Variant under monitoring',
    'previous_VOC': 'Previously designated variant of concern (de-escalated)',
    'previous_VOI': 'Previously designated variant of interest (de-escalated)',
    'not_designated': 'Not a WHO-designated variant',
}

NCBI_GENE_IDS = {
    'orf1ab': '43740578',
    'orf1a': '43740578',
    's': '43740568',
    'orf3a': '43740569',
    'e': '43740570',
    'm': '43740571',
    'orf6': '43740572',
    'orf7a': '43740573',
    'orf7b': '43740574',
    'n': '43740575',
    'orf10': '43740576',
    'orf8': '43740577',
}

GENE_DISPLAY = {
    'orf1ab': 'ORF1ab polyprotein',
    's': 'Surface glycoprotein (spike)',
    'orf3a': 'ORF3a protein',
    'e': 'Envelope protein',
    'm': 'Membrane glycoprotein',
    'orf6': 'ORF6 protein',
    'orf7a': 'ORF7a protein',
    'orf7b': 'ORF7b protein',
    'orf8': 'ORF8 protein',
    'n': 'Nucleocapsid phosphoprotein',
    'orf10': 'ORF10 protein',
}

SO_TERMS = {
    'missense_variant': 'SO:0001583',
    'synonymous_variant': 'SO:0001819',
    'stop_gained': 'SO:0001587',
    'stop_lost': 'SO:0001578',
    'start_lost': 'SO:0002012',
    'frameshift_variant': 'SO:0001589',
    'inframe_insertion': 'SO:0001821',
    'inframe_deletion': 'SO:0001822',
    'conservative_inframe_insertion': 'SO:0001823',
    'disruptive_inframe_insertion': 'SO:0001824',
    'conservative_inframe_deletion': 'SO:0001825',
    'disruptive_inframe_deletion': 'SO:0001826',
    'initiator_codon_variant': 'SO:0001582',
    'stop_retained_variant': 'SO:0001567',
    'splice_region_variant': 'SO:0001630',
    'splice_site_variant': 'SO:0001629',
    'upstream_gene_variant': 'SO:0001631',
    'downstream_gene_variant': 'SO:0001632',
    'intergenic_variant': 'SO:0001628',
    'intergenic_region': 'SO:0000605',
    'intron_variant': 'SO:0001627',
    '5_prime_utr_variant': 'SO:0001623',
    '3_prime_utr_variant': 'SO:0001624',
    'protein_altering_variant': 'SO:0001818',
    'coding_sequence_variant': 'SO:0001580',
    'non_coding_transcript_variant': 'SO:0001619',
    'regulatory_region_variant': 'SO:0001566',
}

IMPACT_DISPLAY = {
    'HIGH': 'High functional impact predicted',
    'MODERATE': 'Moderate functional impact predicted',
    'LOW': 'Low functional impact predicted',
    'MODIFIER': 'Non-coding or modifier variant',
}


def sanitize_id(id_string):
    sanitized = re.sub(r'[^A-Za-z0-9\-]', '-', str(id_string))
    sanitized = re.sub(r'-+', '-', sanitized)
    return sanitized.strip('-') or 'unknown'


def now_iso():
    return datetime.now(timezone.utc).isoformat()


def parse_info(info_str):
    info = {}
    if not info_str or info_str == '.':
        return info
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_ann(info):

    raw = info.get('ANN')
    if not raw or raw is True:
        return {}

    first = raw.split(',')[0]
    fields = first.split('|')

    def at(idx):
        return fields[idx].strip() if len(fields) > idx and fields[idx].strip() else None

    return {
        'effect': at(1),
        'impact': at(2),
        'gene': at(3),
        'gene_id': at(4),
        'hgvs_c': at(9),
        'hgvs_p': at(10),
    }


def parse_vcf(path):
    variants = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            variants.append({
                'chrom': fields[0],
                'pos': int(fields[1]),
                'ref': fields[3],
                'alt': fields[4].split(',')[0],
                'qual': fields[5],
                'filter': fields[6],
                'info': info,
                'ann': parse_ann(info),
            })
    return variants


def load_protein_accessions(gff_path):

    accessions = {}
    if not gff_path or not os.path.exists(gff_path):
        return accessions

    try:
        with open(gff_path) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 9 or fields[2] != 'CDS':
                    continue
                gene = protein = None
                for attr in fields[8].split(';'):
                    if attr.startswith('gene='):
                        gene = attr[5:]
                    elif attr.startswith('protein_id='):
                        protein = attr[11:]
                if gene and protein:
                    accessions.setdefault(gene.lower(), protein)
    except OSError as exc:
        sys.stderr.write('warning: could not read GFF {}: {}\n'.format(gff_path, exc))

    return accessions


def load_json(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError) as exc:
        sys.stderr.write('warning: could not read {}: {}\n'.format(path, exc))
        return {}


def load_lineage(path, sample_id):
    data = load_json(path)
    if not data:
        return {}
    if sample_id in data:
        return data[sample_id]
    if len(data) == 1:
        return list(data.values())[0]
    return {}


def load_consensus_stats(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            lines = [line.rstrip('\n') for line in handle if line.strip()]
        if len(lines) < 2:
            return {}
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        return dict(zip(header, values))
    except OSError:
        return {}


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value, default=None):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def format_ghgvs(accession, pos, ref, alt):

    if len(ref) == 1 and len(alt) == 1:
        return '{}:g.{}{}>{}'.format(accession, pos, ref, alt)

    if len(ref) > len(alt) and ref.startswith(alt):
        start = pos + len(alt)
        end = pos + len(ref) - 1
        if start == end:
            return '{}:g.{}del'.format(accession, start)
        return '{}:g.{}_{}del'.format(accession, start, end)

    if len(alt) > len(ref) and alt.startswith(ref):
        inserted = alt[len(ref):]
        return '{}:g.{}_{}ins{}'.format(accession, pos, pos + 1, inserted)

    end = pos + len(ref) - 1
    if pos == end:
        return '{}:g.{}delins{}'.format(accession, pos, alt)
    return '{}:g.{}_{}delins{}'.format(accession, pos, end, alt)


def coding(system, code, display=None):
    entry = {'system': system, 'code': code}
    if display:
        entry['display'] = display
    return entry


def component(loinc_code, loinc_display, value=None, value_key='valueCodeableConcept'):
    comp = {'code': {'coding': [coding(LOINC, loinc_code, loinc_display)]}}
    if value is not None:
        comp[value_key] = value
    return comp


def narrative(text):
    escaped = (text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))
    return {
        'status': 'generated',
        'div': '<div xmlns="http://www.w3.org/1999/xhtml">{}</div>'.format(escaped),
    }


def lab_categories():
    return [
        {'coding': [coding(OBS_CATEGORY, 'laboratory', 'Laboratory')]},
        {'coding': [coding(V2_0074, 'GE', 'Genetics')]},
    ]


def genomics_meta(profile=None):
    meta = {'tag': [coding(KEMKES_SP, 'genomics', 'Genomics')]}
    if profile:
        meta['profile'] = [profile]
    return meta


def build_variant_observation(idx, record, sample_id, org_id, accession,
                              protein_accessions):
    pos = record['pos']
    ref = record['ref']
    alt = record['alt']
    ann = record['ann']
    info = record['info']

    gene = ann.get('gene')
    gene_key = gene.lower() if gene else None
    effect = ann.get('effect')
    hgvs_c = ann.get('hgvs_c')
    hgvs_p = ann.get('hgvs_p')
    impact = ann.get('impact')

    g_hgvs = format_ghgvs(accession, pos, ref, alt)

    components = [
        component('81290-9', 'Genomic DNA change (gHGVS)',
                  {'coding': [coding(HGVS_SYSTEM, g_hgvs, g_hgvs)]}),
        {'code': {'coding': [coding(NCBI_REFSEQ, accession,
                                    'SARS-CoV-2 isolate Wuhan-Hu-1, complete genome')]}},
    ]

    depth = to_int(info.get('DP'))
    if depth is not None:
        components.append(component(
            '82121-5', 'Allelic read depth',
            {'value': depth, 'unit': 'reads per base pair',
             'system': UCUM, 'code': '1'},
            value_key='valueQuantity'))

    allele_freq = to_float(info.get('AF'))
    if allele_freq is not None:
        components.append(component(
            '81258-6', 'Sample variant allelic frequency',
            {'value': round(allele_freq, 4), 'system': UCUM, 'code': '1'},
            value_key='valueQuantity'))

    if gene:
        gene_code = NCBI_GENE_IDS.get(gene_key, gene)
        components.append(component(
            '48018-6', 'Gene studied [ID]',
            {'coding': [coding(NCBI_GENE, gene_code,
                               GENE_DISPLAY.get(gene_key, gene))],
             'text': gene}))

    if effect:
        so_code = SO_TERMS.get(effect.lower())
        value = {'text': effect}
        if so_code:
            value['coding'] = [coding(SO_SYSTEM, so_code, effect)]
        components.append(component('48019-4', 'DNA change type', value))

    if hgvs_p:
        protein_acc = protein_accessions.get(gene_key) if gene_key else None
        if protein_acc:
            qualified_p = '{}:{}'.format(protein_acc, hgvs_p)
            value = {'coding': [coding(HGVS_SYSTEM, qualified_p, qualified_p)],
                     'text': qualified_p}
        else:
            value = {'text': hgvs_p}
        components.append(component('48005-3', 'Amino acid change (pHGVS)', value))

    if impact and impact.upper() in IMPACT_DISPLAY:
        components.append(component(
            '53037-8', 'Genetic variation clinical significance [Imp]',
            {'text': IMPACT_DISPLAY[impact.upper()]}))

    components.append(component(
        '81254-5', 'Variant exact start-end',
        {'low': {'value': pos}, 'high': {'value': pos + len(ref) - 1}},
        value_key='valueRange'))

    text = 'Variant at position {}: {}>{}'.format(pos, ref, alt)
    if gene:
        text += ' in {}'.format(gene)
    if hgvs_p:
        text += ' ({})'.format(hgvs_p)
    if effect:
        text += ' - {}'.format(effect)

    return {
        'resourceType': 'Observation',
        'id': sanitize_id('{}-obs-{}'.format(sample_id, idx)),
        'meta': genomics_meta(
            'http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant'),
        'text': narrative(text),
        'status': 'final',
        'category': lab_categories(),
        'code': {'coding': [coding(LOINC, '69548-6', 'Genetic variant assessment')]},
        'valueCodeableConcept': {
            'coding': [coding(LOINC, 'LA9633-4', 'Present')],
            'text': 'Present',
        },
        'subject': {'reference': 'Patient/{}-patient'.format(sample_id)},
        'specimen': {'reference': 'Specimen/{}-specimen'.format(sample_id)},
        'effectiveDateTime': now_iso(),
        'performer': [{'reference': 'Organization/{}'.format(org_id)}],
        'component': components,
    }


def build_panel_observation(sample_id, org_id, lineage_info, spike_changes):

    lineage = lineage_info.get('lineage') or 'unassigned'
    clade = lineage_info.get('clade') or 'unassigned'
    who_label = lineage_info.get('who_label')

    components = [
        component('96741-4',
                  'SARS-CoV-2 (COVID-19) variant [Type] in Specimen by Sequencing',
                  {'coding': [coding(SNOMED, '840533007',
                                     'Severe acute respiratory syndrome coronavirus 2')],
                   'text': who_label or 'No WHO-designated variant'}),
        component('96895-8', 'SARS-CoV-2 (COVID-19) lineage [Identifier] in Specimen by Molecular genetics method',
                  {'coding': [coding(PANGO_SYSTEM, lineage, 'Pango lineage {}'.format(lineage))],
                   'text': lineage}),
        component('96896-6', 'SARS-CoV-2 (COVID-19) clade [Type] in Specimen by Molecular genetics method',
                  {'coding': [coding(NEXTSTRAIN_CLADE, clade, 'Nextstrain clade {}'.format(clade))],
                   'text': clade}),
    ]

    if spike_changes:
        components.append(component(
            '96751-3',
            'SARS-CoV-2 (COVID-19) S gene mutation detected [Identifier] in Specimen '
            'by Molecular genetics method',
            {'text': ', '.join(spike_changes)}))

    return {
        'resourceType': 'Observation',
        'id': sanitize_id('{}-sarscov2-panel'.format(sample_id)),
        'meta': genomics_meta(),
        'text': narrative(
            'SARS-CoV-2 sequencing and identification panel for {}: lineage {}, '
            'clade {}'.format(sample_id, lineage, clade)),
        'status': 'final',
        'category': lab_categories(),
        'code': {'coding': [coding(
            LOINC, '96894-1',
            'SARS-CoV-2 (COVID-19) sequencing and identification panel - Specimen '
            'by Molecular genetics method')]},
        'subject': {'reference': 'Patient/{}-patient'.format(sample_id)},
        'specimen': {'reference': 'Specimen/{}-specimen'.format(sample_id)},
        'effectiveDateTime': now_iso(),
        'performer': [{'reference': 'Organization/{}'.format(org_id)}],
        'component': components,
    }


def build_lineage_observation(sample_id, org_id, lineage_info):

    lineage = lineage_info.get('lineage') or 'unassigned'
    provenance = lineage_info.get('provenance') or {}
    source = provenance.get('lineage_source', 'nextclade')
    tool_version = provenance.get('lineage_tool_version', 'unknown')
    db_version = provenance.get('lineage_database_version', 'unknown')

    method_text = '{} {} (database {})'.format(source, tool_version, db_version)

    who_label = lineage_info.get('who_label')
    who_class = lineage_info.get('who_classification')
    confidence = lineage_info.get('confidence', 'unknown')

    text = 'SARS-CoV-2 Pango lineage {} assigned by {}'.format(lineage, method_text)
    if who_label:
        text += '; WHO label {} ({})'.format(who_label, who_class or 'unclassified')
    text += '; assignment confidence {}'.format(confidence)

    notes = [{'text': provenance.get(
        'note', 'Lineage assignment tool and database version recorded in '
                'Observation.method.')}]
    if lineage_info.get('who_status_as_of'):
        notes.append({'text': 'WHO variant classification current as of {}.'.format(
            lineage_info['who_status_as_of'])})

    observation = {
        'resourceType': 'Observation',
        'id': sanitize_id('{}-lineage'.format(sample_id)),
        'meta': genomics_meta(),
        'text': narrative(text),
        'status': 'final',
        'category': lab_categories(),
        'code': {'coding': [coding(
            LOINC, '96895-8',
            'SARS-CoV-2 (COVID-19) lineage [Identifier] in Specimen by Molecular genetics method')]},
        'valueCodeableConcept': {
            'coding': [coding(PANGO_SYSTEM, lineage, 'Pango lineage {}'.format(lineage))],
            'text': lineage,
        },
        'method': {'text': method_text},
        'note': notes,
        'subject': {'reference': 'Patient/{}-patient'.format(sample_id)},
        'specimen': {'reference': 'Specimen/{}-specimen'.format(sample_id)},
        'effectiveDateTime': now_iso(),
        'performer': [{'reference': 'Organization/{}'.format(org_id)}],
    }

    components = []
    if who_label:
        components.append(component(
            '96741-4',
            'SARS-CoV-2 (COVID-19) variant [Type] in Specimen by Sequencing',
            {'text': who_label}))
    if who_class:
        components.append({
            'code': {'text': 'WHO variant classification'},
            'valueCodeableConcept': {
                'coding': [coding(WHO_VARIANT_CLASS_SYSTEM, who_class,
                                  WHO_CLASS_DISPLAY.get(who_class, who_class))],
                'text': who_class,
            },
        })
    if components:
        observation['component'] = components

    return observation


def build_quality_observation(sample_id, org_id, stats, lineage_info):

    if not stats:
        return None

    def text_component(label, value=None, value_key='valueCodeableConcept'):
        comp = {'code': {'text': label}}
        if value is not None:
            comp[value_key] = value
        return comp

    components = []

    length = to_int(stats.get('consensus_length'))
    if length is not None:
        components.append(text_component(
            'Consensus genome length',
            {'value': length, 'unit': 'base pairs', 'system': UCUM, 'code': '1'},
            value_key='valueQuantity'))

    n_percent = to_float(stats.get('n_percent'))
    if n_percent is not None:
        components.append(text_component(
            'Percent of consensus genome masked as N',
            {'value': round(n_percent, 2), 'unit': '%', 'system': UCUM, 'code': '%'},
            value_key='valueQuantity'))

    mean_depth = to_float(stats.get('mean_depth'))
    if mean_depth is not None:
        components.append(text_component(
            'Mean sequencing depth',
            {'value': round(mean_depth, 1), 'unit': 'reads per base pair',
             'system': UCUM, 'code': '1'},
            value_key='valueQuantity'))

    pct_covered = to_float(stats.get('pct_genome_covered'))
    if pct_covered is not None:
        threshold = stats.get('min_coverage_threshold', '10')
        components.append(text_component(
            'Percent of genome covered at >= {}x'.format(threshold),
            {'value': round(pct_covered, 2), 'unit': '%', 'system': UCUM, 'code': '%'},
            value_key='valueQuantity'))

    qc_status = lineage_info.get('qc_status')
    if qc_status:
        components.append(text_component(
            'Sequence quality control status', {'text': qc_status}))

    if not components:
        return None

    return {
        'resourceType': 'Observation',
        'id': sanitize_id('{}-genome-quality'.format(sample_id)),
        'meta': genomics_meta(),
        'text': narrative(
            'Consensus genome quality for {}: {} bp, {}% N, {}x mean depth, '
            '{}% covered at >= {}x'.format(
                sample_id, stats.get('consensus_length', '?'),
                stats.get('n_percent', '?'), stats.get('mean_depth', '?'),
                stats.get('pct_genome_covered', '?'),
                stats.get('min_coverage_threshold', '?'))),
        'status': 'final',
        'category': lab_categories(),
        'code': {'text': 'Consensus genome quality metrics'},
        'subject': {'reference': 'Patient/{}-patient'.format(sample_id)},
        'specimen': {'reference': 'Specimen/{}-specimen'.format(sample_id)},
        'effectiveDateTime': now_iso(),
        'performer': [{'reference': 'Organization/{}'.format(org_id)}],
        'component': components,
    }


def bundle_entry(resource):
    return {
        'fullUrl': 'urn:uuid:{}'.format(str(uuid.uuid4()).lower()),
        'resource': resource,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='snpEff-annotated VCF')
    parser.add_argument('--output', required=True, help='Output FHIR JSON')
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--lineage_json', help='Per-sample lineage JSON from parse_nextclade.py')
    parser.add_argument('--consensus_stats', help='Consensus stats TSV')
    parser.add_argument('--organization_metadata', default='')
    parser.add_argument('--reference_accession', default='NC_045512.2',
                        help='RefSeq accession cited in the bundle')
    parser.add_argument('--gff', help='Reference GFF, for protein accessions')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.stderr.write('error: input VCF not found: {}\n'.format(args.input))
        return 1

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata) or {}
    org_id = org_data.get('org_id') or DEFAULT_ORG_ID

    sample_id = sanitize_id(args.sample_id)
    lineage_info = load_lineage(args.lineage_json, args.sample_id)
    stats = load_consensus_stats(args.consensus_stats)

    protein_accessions = load_protein_accessions(args.gff)
    variants = parse_vcf(args.input)
    sys.stderr.write('{}: {} variants read\n'.format(sample_id, len(variants)))

    entries = []
    observed_spike = []
    annotated = 0

    for idx, record in enumerate(variants, start=1):
        try:
            observation = build_variant_observation(
                idx, record, sample_id, org_id, args.reference_accession,
                protein_accessions)
            entries.append(bundle_entry(observation))

            ann = record['ann']
            if ann.get('gene'):
                annotated += 1
                if ann['gene'].lower() == 's' and ann.get('hgvs_p'):
                    observed_spike.append('S:{}'.format(ann['hgvs_p'].replace('p.', '')))
        except Exception as exc:  
            sys.stderr.write('warning: variant {} skipped: {}\n'.format(idx, exc))

    spike_changes = lineage_info.get('spike_changes') or observed_spike

    entries.append(bundle_entry(
        build_panel_observation(sample_id, org_id, lineage_info, spike_changes)))

    if lineage_info:
        entries.append(bundle_entry(
            build_lineage_observation(sample_id, org_id, lineage_info)))

    quality = build_quality_observation(sample_id, org_id, stats, lineage_info)
    if quality:
        entries.append(bundle_entry(quality))

    bundle = {
        'resourceType': 'Bundle',
        'type': 'collection',
        'timestamp': now_iso(),
        'entry': entries,
    }

    with open(args.output, 'w') as out:
        json.dump(bundle, out, indent=2)

    sys.stderr.write('{}: wrote {} entries ({} variants with gene annotation)\n'.format(
        sample_id, len(entries), annotated))
    return 0


if __name__ == '__main__':
    sys.exit(main())
