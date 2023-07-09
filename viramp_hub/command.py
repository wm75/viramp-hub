import argparse

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
        choices=['ivar-bed', 'cojac-insert-bed', 'amplicon-info'],
        help='Type of output to be generated'
    )
    parser.add_argument(
        '-r', '--report-nested',
        choices=['full', 'inner', 'outer'], default='full',
        help='For amplicons formed by nested primers, report all primers, '
             'or just inner or outer ones. Only applied when writing '
             'amplicon info files (default: "full").'
    )
    parser.add_argument(
        '-f', '--from', choices=['bed'], default='bed',
        help='Format of the input file '
             '(only "bed" is supported in this version)'
    )

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
        if args.to == 'ivar-bed':
            scheme.write_sanitized_bed(out)
        elif args.to == 'cojac-insert-bed':
            scheme.write_insert_bed(out)
        elif args.to == 'amplicon-info':
            mode = args.report_nested or 'full'
            scheme.write_amplicon_info(out, mode)

if __name__ == '__main__':
    main()
