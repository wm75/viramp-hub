# VirAmp-Hub

A tool that lets you manipulate/convert viral amplicon/primer scheme information.

## What is it for?

The layout of (tiled-) amplicon schemes for amplifying and sequencing viral
genomes is typically specified through one of a number of non-standardized
scheme files. These can, generally, not be used interchangeably with downstream
analysis tools.

VirAmp-Hub can convert between commonly used scheme file formats, such as used,
e.g., by the [ARTIC network](https://github.com/artic-network/primer-schemes),
the [ivar](https://github.com/andersen-lab/ivar) and
[cojac](https://github.com/cbg-ethz/cojac) suites of tools, and
[Galaxy](https://galaxyproject.org/) tools and workflows for viral sequencing
data analysis.

## Usage

After installation (`pip install viramp-hub`), use the package via the
`scheme-convert` command.

For help run:

    scheme-convert --help

Prepare an existing primer scheme file for use with ivar:

    scheme-convert primer_scheme.bed --to bed --bed-type ivar -o ivar.bed

Generate an amplicon insert scheme file for cojac:

    scheme-convert primer_scheme.bed --to bed --bed-type cojac -o cojac_insert.bed

Generate an amplicon info file for use with ivar inside Galaxy:

    scheme-convert primer_scheme.bed --to amplicon-info -o info.tsv

Amplicon info files are simple tab-separated files listing primers contributing
to the same amplicon on a single line. For nested primer schemes, you can
control, which primers to report in the info file.

List only the inner primers:

    scheme-convert primer_scheme.bed --to amplicon-info -r inner -o info.tsv

List only the outer primers:

    scheme-convert primer_scheme.bed --to amplicon-info -r outer -o info.tsv

VirAmp-Hub tries to autodetect primer to amplicon assignments using the Python
re pattern: `r'(?P<prefix>(.*_)*)(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)'`.
If this fails for any given primer scheme, it is possible to override
autodetection by providing an amplicon info file directly, e.g.:

    scheme-convert primer_scheme.bed -a amplicon_info.tsv --to bed --bed-type cojac -o cojac_insert.bed

or:

    scheme-convert primer_scheme.bed -a amplicon_info.tsv --to amplicon-info -r inner -o info.tsv

