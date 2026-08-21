#!/usr/bin/env python3

import argparse
import bisect
import os
import re
import subprocess
import sys
import tempfile


AMPLICON_NAME = re.compile(r'^(?P<base>.+?)_(?P<side>LEFT|RIGHT)(?:_.*)?$', re.IGNORECASE)

GENOME_LENGTH = 29903


def load_primer_bounds(bed_path):

    fwd_outer, fwd_inner, rev_outer, rev_inner = [], [], [], []
    left, right = {}, {}

    with open(bed_path) as handle:
        for line in handle:
            if line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue

            name = parts[3] if len(parts) > 3 else ''
            strand = parts[5] if len(parts) > 5 else ''
            if strand not in ('+', '-'):
                strand = '-' if 'RIGHT' in name.upper() else '+'

            if strand == '+':
                fwd_outer.append(start)  
                fwd_inner.append(end)     
            else:
                rev_outer.append(end)
                rev_inner.append(start)

            match = AMPLICON_NAME.match(name)
            if match:
                base = match.group('base')
                if match.group('side').upper() == 'LEFT':
                    left.setdefault(base, start)
                else:
                    right.setdefault(base, end)

    amplicons = [right[b] - left[b] for b in left if b in right and right[b] > left[b]]
    median_amplicon = sorted(amplicons)[len(amplicons) // 2] if amplicons else 0

    return {
        'fwd_outer': sorted(fwd_outer),
        'fwd_inner': sorted(fwd_inner),
        'rev_outer': sorted(rev_outer),
        'rev_inner': sorted(rev_inner),
        'median_amplicon': median_amplicon,
        'n_primers': len(fwd_outer) + len(rev_outer),
        'n_amplicons': len(amplicons),
    }


def near(sorted_positions, value, tolerance):
    if not sorted_positions:
        return False
    idx = bisect.bisect_left(sorted_positions, value)
    for candidate in sorted_positions[max(0, idx - 1): idx + 2]:
        if abs(candidate - value) <= tolerance:
            return True
    return False


def subsample_fastq(fastq, n_reads, workdir):
    out = os.path.join(workdir, 'sub.fastq')
    n_lines = n_reads * 4
    with open(out, 'w') as handle:
        cmd = 'zcat -f {} | head -n {}'.format(fastq, n_lines)
        subprocess.run(['bash', '-c', cmd], stdout=handle, check=False)
    return out


def align(fastq, reference, workdir, preset, threads):
    bam = os.path.join(workdir, 'aln.bam')
    cmd = (
        'minimap2 -a -x {preset} -t {threads} {ref} {fq} 2>/dev/null '
        '| samtools sort -o {bam} - 2>/dev/null'
    ).format(preset=preset, threads=threads, ref=reference, fq=fastq, bam=bam)
    subprocess.run(['bash', '-c', cmd], check=True)
    subprocess.run(['samtools', 'index', bam], check=False)
    return bam


def read_alignment_bounds(bam, min_mapq=1):
    result = subprocess.run(
        ['samtools', 'view', '-F', '0x904', '-q', str(min_mapq), bam],
        capture_output=True, text=True, check=False)

    bounds = []
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        if len(fields) < 6:
            continue
        try:
            pos = int(fields[3]) - 1         
        except ValueError:
            continue
        cigar = fields[5]
        span, num = 0, ''
        for ch in cigar:
            if ch.isdigit():
                num += ch
            else:
                if ch in 'MDN=X' and num:
                    span += int(num)
                num = ''
        if span:
            bounds.append((pos, pos + span))
    return bounds


def score_scheme(bounds, bed_path, tolerance, median_read_len):

    info = load_primer_bounds(bed_path)
    if not info['fwd_outer'] or not bounds:
        return {'both': 0.0, 'either': 0.0, 'amplicon': 0, 'size_ratio': 0.0,
                'n_primers': 0, 'state': 'n/a'}

    total = len(bounds)
    best = None

    for state, fwd_key, rev_key in (('retained', 'fwd_outer', 'rev_outer'),
                                    ('trimmed', 'fwd_inner', 'rev_inner')):
        fwd, rev = info[fwd_key], info[rev_key]
        both = either = 0
        for read_start, read_end in bounds:
            start_hit = near(fwd, read_start, tolerance)
            end_hit = near(rev, read_end, tolerance)
            if start_hit and end_hit:
                both += 1
            if start_hit or end_hit:
                either += 1
        candidate = {'both': 100.0 * both / total,
                     'either': 100.0 * either / total,
                     'state': state}
        if best is None or candidate['both'] > best['both']:
            best = candidate

    median_amplicon = info['median_amplicon']
    if median_amplicon and median_read_len:
        size_ratio = min(median_amplicon, median_read_len) / float(
            max(median_amplicon, median_read_len))
    else:
        size_ratio = 0.0

    window = 2 * tolerance + 1
    p_fwd = min(1.0, len(info['fwd_outer']) * window / float(GENOME_LENGTH))
    p_rev = min(1.0, len(info['rev_outer']) * window / float(GENOME_LENGTH))
    chance_either = 100.0 * (1.0 - (1.0 - p_fwd) * (1.0 - p_rev))

    best.update({
        'amplicon': median_amplicon,
        'size_ratio': size_ratio,
        'n_primers': info['n_primers'],
        'chance': chance_either,
        'enrichment': (best['either'] / chance_either) if chance_either > 0 else 0.0,
    })
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq', required=True)
    parser.add_argument('--schemes', required=True, help='Directory of scheme subdirectories')
    parser.add_argument('--reference', required=True)
    parser.add_argument('--platform', choices=['nanopore', 'illumina'], default='nanopore')
    parser.add_argument('--reads', type=int, default=5000, help='Reads to subsample')
    parser.add_argument('--tolerance', type=int, default=10, help='Boundary tolerance in bp')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()

    preset = 'map-ont' if args.platform == 'nanopore' else 'sr'

    with tempfile.TemporaryDirectory() as workdir:
        sub = subsample_fastq(args.fastq, args.reads, workdir)
        bam = align(sub, args.reference, workdir, preset, args.threads)
        bounds = read_alignment_bounds(bam)

        if not bounds:
            sys.stderr.write('No reads mapped, is this the right reference?\n')
            return 1

        lengths = sorted(e - s for s, e in bounds)
        median_len = lengths[len(lengths) // 2]

        results = []
        for name in sorted(os.listdir(args.schemes)):
            bed = os.path.join(args.schemes, name, 'primer.bed')
            if os.path.exists(bed):
                score = score_scheme(bounds, bed, args.tolerance, median_len)
                score['name'] = name
                results.append(score)

        if not results:
            sys.stderr.write('No primer schemes found under {}\n'.format(args.schemes))
            return 1

        results.sort(key=lambda r: r['enrichment'], reverse=True)

        print('\nSample: {}'.format(os.path.basename(args.fastq)))
        print('Mapped reads analysed: {}   median aligned length: {} bp'.format(
            len(bounds), median_len))
        print('\n{:<16} {:>8} {:>8} {:>8} {:>10} {:>10} {:>6} {:>9}'.format(
            'scheme', 'either', 'chance', 'enrich', 'both', 'amplicon',
            'size', 'primers'))
        print('-' * 88)
        for r in results:
            print('{:<16} {:>7.1f}% {:>7.1f}% {:>7.1f}x {:>9.1f}% {:>9}bp '
                  '{:>6.2f} {:>9}'.format(
                      r['name'], r['either'], r['chance'], r['enrichment'],
                      r['both'], r['amplicon'], r['size_ratio'], r['state']))

        best = results[0]
        runner_up = results[1] if len(results) > 1 else None

        print('')
        if best['enrichment'] < 1.5:
            print('No scheme stands out (best enrichment {:.1f}x over chance).'
                  .format(best['enrichment']))
            print('Either this is metagenomic, or the real scheme is not bundled '
                  'here. Set primer_scheme=none: clipping with the wrong scheme '
                  'removes real sequence and leaves the actual primer bases in '
                  'place, which is worse than not clipping.')
        elif runner_up and best['enrichment'] < 1.3 * runner_up['enrichment']:
            print('Ambiguous: {} ({:.1f}x) does not clearly beat {} ({:.1f}x).'
                  .format(best['name'], best['enrichment'],
                          runner_up['name'], runner_up['enrichment']))
            print('Prefer the scheme whose amplicon length matches the {}bp '
                  'median read, and cross-check the sample\'s publication date '
                  'against the scheme\'s release date.'.format(median_len))
        else:
            print('Best match: {}  ({:.1f}x enrichment over chance, {}bp '
                  'amplicons vs {}bp median read, primers appear {})'.format(
                      best['name'], best['enrichment'], best['amplicon'],
                      median_len, best['state']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
