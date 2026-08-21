#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone


def load_who_table(path):
    """who_label -> {classification, status_as_of, notes}"""
    table = {}
    if not path or not os.path.exists(path):
        return table
    with open(path, newline='') as handle:
        rows = (line for line in handle if not line.startswith('#'))
        for row in csv.DictReader(rows):
            label = (row.get('who_label') or '').strip()
            if label:
                table[label.lower()] = {
                    'classification': (row.get('classification') or '').strip(),
                    'status_as_of': (row.get('status_as_of') or '').strip(),
                    'notes': (row.get('notes') or '').strip(),
                }
    return table


def read_dataset_version(dataset_dir):
    path = os.path.join(dataset_dir, 'pathogen.json')
    if not os.path.exists(path):
        return 'unknown'
    try:
        with open(path) as handle:
            data = json.load(handle)
        version = data.get('version') or {}
        return version.get('tag') or data.get('tag') or 'unknown'
    except (OSError, ValueError):
        return 'unknown'


def read_nextclade_row(tsv_path, sample):
    if not os.path.exists(tsv_path):
        return {}
    with open(tsv_path, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows:
        return {}
    for row in rows:
        if (row.get('seqName') or '').strip() == sample:
            return row
    return rows[0]


def read_pangolin_csv(path, sample):
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle):
            taxon = (row.get('taxon') or row.get('Sequence name') or '').strip()
            if taxon == sample or not taxon:
                return row
    return {}


def split_list(value):
    if not value:
        return []
    return [item for item in value.replace(' ', '').split(',') if item]


def spike_changes(aa_changes):
    return [c for c in aa_changes if c.upper().startswith('S:')]


def derive_confidence(row, n_percent=None):

    status = (row.get('qc.overallStatus') or '').strip().lower()
    base = {
        'good': 'high',
        'mediocre': 'medium',
        'bad': 'low',
    }.get(status, 'unknown')

    coverage = row.get('coverage')
    try:
        coverage = float(coverage) if coverage not in (None, '') else None
    except ValueError:
        coverage = None

    if coverage is not None:
        if coverage < 0.5 and base != 'low':
            return 'low'
        if coverage < 0.9 and base == 'high':
            return 'medium'
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextclade-tsv', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--dataset', help='Nextclade dataset directory')
    parser.add_argument('--who-table', help='WHO variant classification CSV')
    parser.add_argument('--nextclade-version', default='unknown')
    parser.add_argument('--pangolin_csv', help='Optional pangolin output; supersedes Nextclade')
    args = parser.parse_args()

    row = read_nextclade_row(args.nextclade_tsv, args.sample)
    if not row:
        sys.stderr.write(
            'warning: no Nextclade result for {}; emitting an empty call\n'.format(args.sample))

    who_table = load_who_table(args.who_table)
    dataset_tag = read_dataset_version(args.dataset) if args.dataset else 'unknown'

    clade = (row.get('clade') or '').strip()
    who_label = (row.get('clade_who') or '').strip()
    lineage = (row.get('Nextclade_pango') or '').strip()

    lineage_source = 'nextclade'
    lineage_tool_version = args.nextclade_version
    lineage_db_version = dataset_tag

    pangolin_row = read_pangolin_csv(args.pangolin_csv, args.sample)
    if pangolin_row:
        pango_lineage = (pangolin_row.get('lineage') or '').strip()
        if pango_lineage:
            lineage = pango_lineage
            lineage_source = 'pangolin'
            lineage_tool_version = (pangolin_row.get('pangolin_version') or 'unknown').strip()
            lineage_db_version = (pangolin_row.get('version')
                                  or pangolin_row.get('pangolin_data_version')
                                  or 'unknown').strip()

    who_info = who_table.get(who_label.lower(), {})

    aa_subs = split_list(row.get('aaSubstitutions'))
    aa_dels = split_list(row.get('aaDeletions'))
    aa_all = aa_subs + aa_dels

    result = {
        'sample_id': args.sample,
        'lineage': lineage or 'unassigned',
        'clade': clade or 'unassigned',
        'clade_display': (row.get('clade_display') or '').strip(),
        'clade_nextstrain': (row.get('clade_nextstrain') or '').strip(),
        'who_label': who_label or None,
        'who_classification': who_info.get('classification') or (
            'not_designated' if not who_label else 'unknown'),
        'who_status_as_of': who_info.get('status_as_of') or None,
        'who_notes': who_info.get('notes') or None,
        'confidence': derive_confidence(row),
        'qc_status': (row.get('qc.overallStatus') or '').strip() or 'unknown',
        'qc_score': row.get('qc.overallScore') or None,
        'coverage': row.get('coverage') or None,
        'total_missing': row.get('totalMissing') or None,
        'total_substitutions': row.get('totalSubstitutions') or None,
        'total_deletions': row.get('totalDeletions') or None,
        'frame_shifts': split_list(row.get('frameShifts')),
        'nucleotide_substitutions': split_list(row.get('substitutions')),
        'aa_substitutions': aa_subs,
        'aa_deletions': aa_dels,
        'spike_changes': spike_changes(aa_all),
        'provenance': {
            'lineage_source': lineage_source,
            'lineage_tool_version': lineage_tool_version,
            'lineage_database_version': lineage_db_version,
            'clade_source': 'nextclade',
            'clade_tool_version': args.nextclade_version,
            'clade_database_version': dataset_tag,
            'generated_at': datetime.now(timezone.utc).isoformat(),
            'note': (
                'Pango lineage assigned by Nextclade tree placement '
                '(Nextclade_pango), not by pangolin/UShER.'
                if lineage_source == 'nextclade' else
                'Pango lineage assigned by pangolin.'
            ),
        },
    }

    with open(args.output, 'w') as out:
        json.dump({args.sample: result}, out, indent=2)

    sys.stderr.write('{}: lineage={} clade={} who={} confidence={} (via {})\n'.format(
        args.sample, result['lineage'], result['clade'],
        result['who_label'] or 'none', result['confidence'], lineage_source))


if __name__ == '__main__':
    main()
