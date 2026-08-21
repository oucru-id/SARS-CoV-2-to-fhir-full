#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import OrderedDict

VCF_HEADER = """##fileformat=VCFv4.2
##source=ivar_variants_to_vcf.py
##contig=<ID={contig}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">
##INFO=<ID=REF_DP,Number=1,Type=Integer,Description="Depth supporting the reference allele">
##INFO=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth supporting the alternate allele">
##INFO=<ID=ALT_QUAL,Number=1,Type=Float,Description="Mean quality of bases supporting the alternate allele">
##FILTER=<ID=ft,Description="Fisher's exact test p-value above threshold, as reported by ivar">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def expand_indel(pos, ref, alt):
    if alt.startswith('+'):
        return pos, ref, ref + alt[1:]
    if alt.startswith('-'):
        return pos, ref + alt[1:], ref
    return pos, ref, alt


def parse_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='ivar variants TSV')
    parser.add_argument('--output', required=True, help='Output VCF path')
    parser.add_argument('--contig', required=True, help='Reference contig name')
    parser.add_argument('--sample', required=True, help='Sample id for the genotype column')
    parser.add_argument('--pass-only', action='store_true',
                        help='Emit only records ivar marked PASS')
    args = parser.parse_args()

    records = OrderedDict()

    with open(args.input, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or 'POS' not in reader.fieldnames:
            sys.stderr.write('ivar table is empty or headerless; writing empty VCF\n')
            reader = []

        for row in reader:
            pos = parse_int(row.get('POS'))
            ref = (row.get('REF') or '').strip()
            alt = (row.get('ALT') or '').strip()
            if not pos or not ref or not alt:
                continue

            ivar_pass = (row.get('PASS') or '').strip().upper() in ('TRUE', 'T', '1')
            if args.pass_only and not ivar_pass:
                continue

            pos, ref, alt = expand_indel(pos, ref, alt)

            total_dp = parse_int(row.get('TOTAL_DP'))
            ref_dp = parse_int(row.get('REF_DP'))
            alt_dp = parse_int(row.get('ALT_DP'))
            alt_freq = parse_float(row.get('ALT_FREQ'))
            alt_qual = parse_float(row.get('ALT_QUAL'))

            key = (pos, ref, alt)
            if key in records:
                continue

            info = (
                'DP={dp};AF={af:.4f};REF_DP={ref_dp};ALT_DP={alt_dp};ALT_QUAL={qual:.1f}'
            ).format(dp=total_dp, af=alt_freq, ref_dp=ref_dp, alt_dp=alt_dp, qual=alt_qual)

            records[key] = {
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': '{:.1f}'.format(alt_qual),
                'filter': 'PASS' if ivar_pass else 'ft',
                'info': info,
                'sample': '1:{dp}:{ref_dp},{alt_dp}'.format(
                    dp=total_dp, ref_dp=ref_dp, alt_dp=alt_dp),
            }

    with open(args.output, 'w') as out:
        out.write(VCF_HEADER.format(contig=args.contig, sample=args.sample))
        for rec in sorted(records.values(), key=lambda r: (r['pos'], r['ref'], r['alt'])):
            out.write('\t'.join([
                args.contig,
                str(rec['pos']),
                '.',
                rec['ref'],
                rec['alt'],
                rec['qual'],
                rec['filter'],
                rec['info'],
                'GT:DP:AD',
                rec['sample'],
            ]) + '\n')

    sys.stderr.write('Wrote {} variant records to {}\n'.format(len(records), args.output))


if __name__ == '__main__':
    main()
