import re

from typing import NamedTuple, Literal


class Primer(NamedTuple):
    start: int
    end: int
    name: str
    strand: Literal['+', '-']


class Amplicon:
    def __init__(self, primers, ref_id, pool_id):
        self.primers = primers
        self.ref_id = ref_id
        self.pool_id = pool_id

    @property
    def outer_fw(self):
        return min(
            (primer for primer in self.primers if primer.strand == '+'),
            key=lambda x: x.start
        )

    @property
    def outer_rv(self):
        return max(
            (primer for primer in self.primers if primer.strand == '-'),
            key=lambda x: x.end
        )

    @property
    def inner_fw(self):
        return max(
            (primer for primer in self.primers if primer.strand == '+'),
            key=lambda x: x.start
        )

    @property
    def inner_rv(self):
        return min(
            (primer for primer in self.primers if primer.strand == '-'),
            key=lambda x: x.end
        )

    def append_primer(self, primer, ref_id, pool_id):
        """Safely add a primer to an amplicon.

        The required ref_id and pool_id are checked for a match with the
        amplicon's attributes and a ValueError gets raised if the primer cannot
        belong to the amplicon based on them.
        """

        if self.ref_id != ref_id:
            raise ValueError(
                'Reference mismatch between primer "{0}" and existing amplicon: '
                'Primer associated with reference "{1}", amplicon with "{2}".'
                .format(primer.name, ref_id, self.ref_id)
            )
        if self.pool_id != pool_id:
            raise ValueError(
                'Pool mismatch between primer "{0}" and existing amplicon: '
                'Primer associated with pool "{1}", amplicon with "{2}".'
                .format(primer.name, pool_id, self.pool_id)
            )
        self.primers.append(primer)


class Scheme:
    amplicon_pat = re.compile(
        r'(?P<prefix>(.*_)*)(?P<num>\d+).*_(?P<name>L(?:EFT)?|R(?:IGHT)?)'
    )

    @classmethod
    def infer_from_primer_scheme(cls, primer_scheme):
        amplicons = {}
        prefixes_seen = set()
        
        for primer_dat, ref_id, pool_id in cls.read_primer_bed(primer_scheme):
            re_match = cls.amplicon_pat.match(primer_dat.name)
            if re_match is None:
                raise ValueError(
                    '{} does not match expected amplicon name format'
                    .format(primer_dat.name)
                )
            prefix = re_match.group('prefix')[:-1]
            prefixes_seen.add(prefix)
            amplicon_id = int(re_match.group('num'))
            if amplicon_id in amplicons:
                amplicons[amplicon_id].append_primer(
                    primer_dat, ref_id, pool_id
                )
            else:
                amplicons[amplicon_id] = Amplicon(
                    [primer_dat], ref_id, pool_id
                )

        if len(prefixes_seen) == 1 and prefix:
            scheme_name = prefix
        else:
            scheme_name = None
            
        return cls(amplicons, scheme_name)

    @classmethod
    def from_primers_and_amplicons(
        cls, primer_scheme, amplicon_info, scheme_name=None
    ):
        primer_amplicon_mapping = {}
        amplicon_id = 1
        for line in amplicon_info:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            for primer_name in line.split('\t'):
                primer_amplicon_mapping[primer_name] = amplicon_id
            amplicon_id += 1    

        amplicons = {}
        for primer_dat, ref_id, pool_id in cls.read_primer_bed(primer_scheme):
            try:
                mapped_id = primer_amplicon_mapping[primer_dat.name]
            except KeyError:
                raise ValueError(
                    'BED file primer: "{}" not listed in amplicon info!'
                    .format(primer_dat.name)
                )
            if mapped_id in amplicons:
                amplicons[mapped_id].append_primer(
                    primer_dat, ref_id, pool_id
                )
            else:
                amplicons[mapped_id] = Amplicon(
                    [primer_dat], ref_id, pool_id
                )

        return cls(amplicons, scheme_name)

    @classmethod
    def read_primer_bed(cls, primer_bed):
        for record in primer_bed:
            if not record.strip() or record[0] == '#':
                continue
            fields = record.strip('\n').split('\t')
            
            primer_dat = Primer(
                start = int(fields[1]),
                end = int(fields[2]),
                name = fields[3],
                strand = fields[5]
            )
            ref_id = fields[0]
            pool_id = fields[4]

            yield primer_dat, ref_id, pool_id

    def __init__(self, amplicons, name=None):
        self.amplicons = amplicons
        self.name = name

    def write_sanitized_bed(self, o):
        records = sorted(
            ((primer, amplicon.ref_id, amplicon.pool_id)
            for amplicon in self.amplicons.values()
            for primer in amplicon.primers),
            # sort by ref, start, end
            key=lambda x: (x[1], x[0].start, x[0].end)
        )
        for primer, ref_id, pool_id in records:
            o.write(
                f'{ref_id}\t{primer.start}\t{primer.end}\t{primer.name}\t60\t{primer.strand}\n'
            )
            
    def write_amplicon_info(self, o, mode='full'):
        for amplicon in self.amplicons.values():
            if mode == 'full':
                names = [primer_dat.name for primer_dat in amplicon.primers]
            elif mode == 'outer':
                names = [amplicon.outer_fw.name, amplicon.outer_rv.name]
            elif mode == 'inner':
                names = [amplicon.inner_fw.name, amplicon.inner_rv.name]
            o.write('\t'.join(names) + '\n')

    def write_insert_bed(self, o):
        for amplicon_id, amplicon in sorted(self.amplicons.items()):
            insert_start = amplicon.inner_fw.end
            insert_end = amplicon.inner_rv.start
            if self.name:
                insert_name = f'{self.name}_INSERT_{amplicon_id}'
            else:
                insert_name = f'INSERT_{amplicon_id}'
            o.write(
                f'{amplicon.ref_id}\t{insert_start}\t{insert_end}\t{insert_name}\t{amplicon.pool_id}\t+\n'
                )

    def write_bedpe(self, o):
        """Write the primer scheme as 11-column BEDPE.

        The amplicon name will be used as the value of column 7 and
        the pool ID will be written to column 11. If there is more
        than one primer pair defining an amplicon all (coordinate-sorted)
        combinations of fw and rv primers will be written as separate lines
        with increasing index numbers in column 8 (as with
        https://github.com/rki-mf1/CoVpipe2).
        See https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
        for a definition of the BEDPE format."""

        # First, generate a totally sorted copy of the amplicons,
        # in which the primers for each amplicon are sorted by start, stop,
        # and where also the amplicons are sorted by their outermost
        # coordinates
        coord_sorted_amplicons = sorted((
            (
                amplicon_id,
                Amplicon(
                    sorted(
                        amplicon.primers,
                        key=lambda x: (x.start, x.end)
                    ),
                    amplicon.ref_id,
                    amplicon.pool_id
                )
            ) for amplicon_id, amplicon in self.amplicons.items()),
            key=lambda x: (x[1].ref_id, x[1].primers[0].start, x[1].primers[-1].end)
        )

        for amplicon_id, amplicon in coord_sorted_amplicons:
            fw_primers = [primer_dat for primer_dat in amplicon.primers if primer_dat[3] == '+']
            rv_primers = [primer_dat for primer_dat in amplicon.primers if primer_dat[3] == '-']
            pair_index = 0
            if self.name:
                amplicon_name = f'{self.name}_AMPLICON_{amplicon_id}'
            else:
                amplicon_name = f'AMPLICON_{amplicon_id}'
            for fw_p in fw_primers:
                for rv_p in rv_primers:
                    o.write(
                        f'{amplicon.ref_id}\t{fw_p.start}\t{fw_p.end}\t'
                        f'{amplicon.ref_id}\t{rv_p.start}\t{rv_p.end}\t'
                        f'{amplicon_name}\t{pair_index}\t{fw_p.strand}\t{rv_p.strand}\t'
                        f'{amplicon.pool_id}\n'
                    )
                    pair_index += 1
