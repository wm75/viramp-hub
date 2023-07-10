import argparse
import sys

from . import notice
from .scheme import Scheme


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Input file')
    parser.add_argument(
        '-o', '--output', required=True,
        help='Name of the output file'
    )
    parser.add_argument(
        '-a', '--amplicon-info', default=argparse.SUPPRESS,
        help='Amplicon info file for inputs from which primer groupings '
             'cannot be infered.'
    )
    parser.add_argument(
        '-t', '--to', required=True,
        choices=['bed', 'amplicon-info'],
        help='Type of output to be generated'
    )
    parser.add_argument(
        '-b', '--bed-type',
        choices=['ivar', 'cojac'],
        help='For "bed" output, the type of bed to be written; '
             'Currently, you can specify "ivar" to produce primer bed output '
             'compatible with the ivar suite of tools, or "cojac" to generate '
             'the amplicon insert bed expected by cojac.'
    )
    parser.add_argument(
        '-r', '--report-nested',
        choices=['full', 'inner', 'outer'], default='full',
        help='For amplicons formed by nested primers, report all primers, '
             'or just inner or outer ones. Applied only when writing '
             'amplicon info files (default: "full").'
    )
    parser.add_argument(
        '-f', '--from', choices=['bed'], default='bed',
        help='Format of the input file '
             '(only "bed" is supported in this version)'
    )

    if len(sys.argv)<2:
        print(notice, file=sys.stderr)
        print('\nPlease run with -h / --help for help.', file=sys.stderr)
        sys.exit(2)
    args = parser.parse_args()

    with open(args.input) as input_data:
        if 'amplicon_info' in args:
            with open(args.amplicon_info) as amplicon_info:
                scheme = Scheme.from_primers_and_amplicons(
                    input_data, amplicon_info
                )
        else:
            scheme = Scheme.infer_from_primer_scheme(
                input_data
            )

    with open(args.output, 'w') as out:
        if args.to == 'bed':
            if not args.bed_type:
                sys.exit(
                    'Please specify the type of bed output to be produced '
                    'through the -b / --bed-type option.'
                )
            if args.bed_type == 'ivar':
                scheme.write_sanitized_bed(out)
            elif args.bed_type == 'cojac':
                scheme.write_insert_bed(out)
        elif args.to == 'amplicon-info':
            scheme.write_amplicon_info(out, args.report_nested)

if __name__ == '__main__':
    main()
