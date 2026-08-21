#!/usr/bin/env python3

import argparse
import subprocess
import sys

FIELDS = [
    'sample',
    'consensus_length',
    'n_count',
    'n_percent',
    'ambiguous_count',
    'mean_depth',
    'median_depth',
    'pct_genome_covered',
    'reads_mapped',
    'min_coverage_threshold',
]

UNAMBIGUOUS = set('ACGT')


def read_fasta_sequence(path):
    seq = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq).upper()


def run(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout
    except (OSError, subprocess.SubprocessError) as exc:
        sys.stderr.write('warning: {} failed: {}\n'.format(cmd[0], exc))
        return ''


def depth_stats(bam, min_coverage):
    output = run(['samtools', 'depth', '-a', bam])
    depths = []
    for line in output.splitlines():
        parts = line.split('\t')
        if len(parts) >= 3:
            try:
                depths.append(int(parts[2]))
            except ValueError:
                continue

    if not depths:
        return 0.0, 0, 0.0

    depths.sort()
    mean = sum(depths) / len(depths)
    median = depths[len(depths) // 2]
    covered = sum(1 for d in depths if d >= min_coverage)
    pct = 100.0 * covered / len(depths)
    return mean, median, pct


def mapped_reads(bam):
    output = run(['samtools', 'flagstat', bam])
    for line in output.splitlines():
        if ' mapped (' in line and 'primary' not in line:
            try:
                return int(line.split()[0])
            except (ValueError, IndexError):
                return 0
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True, help='Consensus FASTA')
    parser.add_argument('--bam', required=True, help='Alignment BAM')
    parser.add_argument('--sample', required=True, help='Sample id')
    parser.add_argument('--min-coverage', type=int, default=10)
    parser.add_argument('--output', required=True, help='Output TSV')
    args = parser.parse_args()

    seq = read_fasta_sequence(args.fasta)
    length = len(seq)
    n_count = seq.count('N')
    ambiguous = sum(1 for b in seq if b not in UNAMBIGUOUS and b != 'N')
    n_percent = (100.0 * n_count / length) if length else 0.0

    mean_depth, median_depth, pct_covered = depth_stats(args.bam, args.min_coverage)

    values = {
        'sample': args.sample,
        'consensus_length': length,
        'n_count': n_count,
        'n_percent': '{:.2f}'.format(n_percent),
        'ambiguous_count': ambiguous,
        'mean_depth': '{:.1f}'.format(mean_depth),
        'median_depth': median_depth,
        'pct_genome_covered': '{:.2f}'.format(pct_covered),
        'reads_mapped': mapped_reads(args.bam),
        'min_coverage_threshold': args.min_coverage,
    }

    with open(args.output, 'w') as out:
        out.write('\t'.join(FIELDS) + '\n')
        out.write('\t'.join(str(values[f]) for f in FIELDS) + '\n')

    sys.stderr.write(
        '{}: {} bp, {:.2f}% N, {:.1f}x mean depth, {:.2f}% >= {}x\n'.format(
            args.sample, length, n_percent, mean_depth, pct_covered, args.min_coverage))


if __name__ == '__main__':
    main()
