#!/usr/bin/env python3

import argparse
import json
import os
import sys
from collections import OrderedDict

WIDTH = 78
GENE_ORDER = ['S', 'ORF1ab', 'orf1ab', 'N', 'M', 'E', 'ORF3a', 'ORF6',
              'ORF7a', 'ORF7b', 'ORF8', 'ORF10']

NOTABLE_EFFECTS = {
    'missense_variant', 'stop_gained', 'stop_lost', 'start_lost',
    'frameshift_variant', 'disruptive_inframe_deletion',
    'conservative_inframe_deletion', 'disruptive_inframe_insertion',
    'conservative_inframe_insertion',
}


def rule(char='='):
    return char * WIDTH


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


def parse_annotated_vcf(path):
    variants = []
    if not os.path.exists(path):
        return variants

    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue

            info = parse_info(fields[7])
            ann_raw = info.get('ANN')
            gene = effect = hgvs_c = hgvs_p = impact = None

            if ann_raw and ann_raw is not True:
                parts = ann_raw.split(',')[0].split('|')

                def at(i):
                    return parts[i].strip() if len(parts) > i and parts[i].strip() else None

                effect, impact, gene = at(1), at(2), at(3)
                hgvs_c, hgvs_p = at(9), at(10)

            variants.append({
                'pos': int(fields[1]),
                'ref': fields[3],
                'alt': fields[4].split(',')[0],
                'gene': gene or 'intergenic',
                'effect': effect or 'unknown',
                'impact': impact or '',
                'hgvs_c': hgvs_c or '',
                'hgvs_p': hgvs_p or '',
                'depth': info.get('DP', ''),
                'af': info.get('AF', ''),
            })
    return variants


def load_lineage(path, sample_id):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return {}
    if sample_id in data:
        return data[sample_id]
    if len(data) == 1:
        return list(data.values())[0]
    return {}


def load_stats(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            lines = [l.rstrip('\n') for l in handle if l.strip()]
        if len(lines) < 2:
            return {}
        return dict(zip(lines[0].split('\t'), lines[1].split('\t')))
    except OSError:
        return {}


def quality_verdict(stats, lineage):

    try:
        n_percent = float(stats.get('n_percent', 100))
    except (TypeError, ValueError):
        n_percent = 100.0

    qc = (lineage.get('qc_status') or '').lower()

    if n_percent > 50:
        return ('UNRELIABLE',
                'More than half the genome is missing. Treat the lineage '
                'assignment below as unusable.')
    if n_percent > 10:
        return ('CAUTION',
                'Substantial genome loss ({:.1f}% N). The lineage assignment may '
                'be imprecise, particularly for closely related sublineages.'
                .format(n_percent))
    if qc == 'bad':
        return ('CAUTION',
                'Nextclade flagged this genome as failing QC despite adequate '
                'coverage; check the QC detail before reporting.')
    return ('ACCEPTABLE',
            'Genome recovery is sufficient to support the lineage assignment.')


def write_report(out, sample_id, args, variants, lineage, stats):
    out.write(rule() + '\n')
    out.write('SARS-CoV-2 GENOMIC SURVEILLANCE REPORT - {}\n'.format(sample_id))
    out.write(rule() + '\n\n')

    out.write('Platform      : {}\n'.format(args.platform))
    out.write('Protocol      : {}\n'.format(args.protocol))
    out.write('Primer scheme : {}\n'.format(
        args.primer_scheme if args.primer_scheme != 'none'
        else 'none (primer clipping not applied)'))
    out.write('Reference     : {}\n\n'.format(args.reference))

    out.write(rule('-') + '\n')
    out.write('1. CONSENSUS GENOME QUALITY\n')
    out.write(rule('-') + '\n\n')

    if stats:
        out.write('  Consensus length        : {} bp\n'.format(
            stats.get('consensus_length', 'n/a')))
        out.write('  Masked (N) bases        : {} ({}%)\n'.format(
            stats.get('n_count', 'n/a'), stats.get('n_percent', 'n/a')))
        out.write('  Ambiguous (IUPAC) bases : {}\n'.format(
            stats.get('ambiguous_count', 'n/a')))
        out.write('  Mean depth              : {}x\n'.format(
            stats.get('mean_depth', 'n/a')))
        out.write('  Median depth            : {}x\n'.format(
            stats.get('median_depth', 'n/a')))
        out.write('  Genome covered >= {}x   : {}%\n'.format(
            stats.get('min_coverage_threshold', '?'),
            stats.get('pct_genome_covered', 'n/a')))
        out.write('  Reads mapped            : {}\n'.format(
            stats.get('reads_mapped', 'n/a')))
    else:
        out.write('  No consensus statistics available.\n')

    if lineage.get('qc_status'):
        out.write('  Nextclade QC            : {}\n'.format(lineage['qc_status']))

    verdict, explanation = quality_verdict(stats, lineage)
    out.write('\n  VERDICT: {}\n'.format(verdict))
    for line in _wrap(explanation, WIDTH - 4):
        out.write('  {}\n'.format(line))
    out.write('\n')

    out.write(rule('-') + '\n')
    out.write('2. LINEAGE AND CLADE\n')
    out.write(rule('-') + '\n\n')

    if lineage:
        out.write('  Pango lineage     : {}\n'.format(lineage.get('lineage', 'unassigned')))
        out.write('  Nextstrain clade  : {}\n'.format(lineage.get('clade', 'unassigned')))
        out.write('  WHO label         : {}\n'.format(lineage.get('who_label') or 'none'))
        out.write('  WHO classification: {}\n'.format(
            lineage.get('who_classification') or 'not designated'))
        if lineage.get('who_status_as_of'):
            out.write('  WHO status as of  : {}\n'.format(lineage['who_status_as_of']))
        out.write('  Call confidence   : {}\n'.format(lineage.get('confidence', 'unknown')))

        prov = lineage.get('provenance') or {}
        out.write('\n  Assigned by       : {} {} (database {})\n'.format(
            prov.get('lineage_source', 'unknown'),
            prov.get('lineage_tool_version', '?'),
            prov.get('lineage_database_version', '?')))
        if prov.get('note'):
            for line in _wrap(prov['note'], WIDTH - 4):
                out.write('  {}\n'.format(line))
    else:
        out.write('  No lineage assignment available.\n')
    out.write('\n')

    out.write(rule('-') + '\n')
    out.write('3. SPIKE (S) PROTEIN CHANGES\n')
    out.write(rule('-') + '\n\n')

    spike = lineage.get('spike_changes') or []
    if spike:
        out.write('  {} amino-acid change(s) in spike:\n\n'.format(len(spike)))
        for i in range(0, len(spike), 4):
            out.write('    {}\n'.format('  '.join(spike[i:i + 4])))
    else:
        spike_variants = [v for v in variants
                          if v['gene'].lower() == 's' and v['hgvs_p']]
        if spike_variants:
            out.write('  {} spike change(s) from the annotated VCF:\n\n'.format(
                len(spike_variants)))
            for v in spike_variants:
                out.write('    {:<8} {:<24} {}\n'.format(
                    v['pos'], v['hgvs_p'], v['effect']))
        else:
            out.write('  No spike protein changes detected.\n')
    out.write('\n')

    out.write(rule('-') + '\n')
    out.write('4. VARIANTS BY GENE\n')
    out.write(rule('-') + '\n\n')

    if not variants:
        out.write('  No variants called.\n\n')
    else:
        by_gene = OrderedDict()
        for v in variants:
            by_gene.setdefault(v['gene'], []).append(v)

        def gene_sort_key(name):
            try:
                return (0, GENE_ORDER.index(name))
            except ValueError:
                return (1, name)

        for gene in sorted(by_gene, key=gene_sort_key):
            gene_variants = by_gene[gene]
            notable = sum(1 for v in gene_variants
                          if v['effect'] in NOTABLE_EFFECTS)
            out.write('  {} ({} variant{}, {} protein-altering)\n'.format(
                gene, len(gene_variants),
                '' if len(gene_variants) == 1 else 's', notable))
            out.write('  {}\n'.format('-' * (WIDTH - 4)))
            out.write('  {:<8} {:<12} {:<22} {:<26} {:>6}\n'.format(
                'POS', 'CHANGE', 'PROTEIN', 'EFFECT', 'DEPTH'))
            for v in sorted(gene_variants, key=lambda x: x['pos']):
                change = '{}>{}'.format(v['ref'], v['alt'])
                if len(change) > 11:
                    change = change[:8] + '...'
                out.write('  {:<8} {:<12} {:<22} {:<26} {:>6}\n'.format(
                    v['pos'], change, v['hgvs_p'] or '-',
                    v['effect'][:26], v['depth'] or '-'))
            out.write('\n')

    out.write(rule('-') + '\n')
    out.write('5. SUMMARY\n')
    out.write(rule('-') + '\n\n')

    protein_altering = sum(1 for v in variants if v['effect'] in NOTABLE_EFFECTS)
    out.write('  Total variants          : {}\n'.format(len(variants)))
    out.write('  Protein-altering        : {}\n'.format(protein_altering))
    out.write('  Spike changes           : {}\n'.format(
        len(spike) if spike else
        len([v for v in variants if v['gene'].lower() == 's' and v['hgvs_p']])))
    out.write('  Lineage                 : {}\n'.format(
        lineage.get('lineage', 'unassigned')))
    out.write('  Genome quality          : {}\n'.format(verdict))
    out.write('\n' + rule() + '\n')


def _wrap(text, width):
    words = text.split()
    lines, current = [], ''
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(current)
            current = word
        else:
            current = '{} {}'.format(current, word).strip()
    if current:
        lines.append(current)
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotated_vcf', required=True)
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--lineage_json')
    parser.add_argument('--consensus_stats')
    parser.add_argument('--platform', default='unknown')
    parser.add_argument('--protocol', default='unknown')
    parser.add_argument('--primer_scheme', default='none')
    parser.add_argument('--reference', default='MN908947.3 (Wuhan-Hu-1)')
    parser.add_argument('--output_dir', default='.')
    args = parser.parse_args()

    variants = parse_annotated_vcf(args.annotated_vcf)
    lineage = load_lineage(args.lineage_json, args.sample_id)
    stats = load_stats(args.consensus_stats)

    out_path = os.path.join(args.output_dir,
                            '{}.summary_report.txt'.format(args.sample_id))
    with open(out_path, 'w') as out:
        write_report(out, args.sample_id, args, variants, lineage, stats)

    sys.stderr.write('Wrote {} ({} variants)\n'.format(out_path, len(variants)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
