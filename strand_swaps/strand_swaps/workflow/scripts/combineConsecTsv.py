#!/usr/bin/env python3


def parseTsv(file=None, category=None):
    sample = file.replace('.tsv','')
    with open(file, 'r') as f:
        for line in f:
            l = line.strip()
            yield f'{sample}\t{category}\t{l}\n'


with open(snakemake.output, 'w') as out:
    for file in snakemake.input.strand:
        for line in parseTsv(file=file, category='strand'):
            out.write(line)
    for file in snakemake.input.frame:
        for line in parseTsv(file=file, category='frame'):
            out.write(line)
    for file in snakemake.input.overlap:
        for line in parseTsv(file=file, category='overlap'):
            out.write(line)
