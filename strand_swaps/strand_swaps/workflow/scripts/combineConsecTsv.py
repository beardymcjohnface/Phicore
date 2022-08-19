#!/usr/bin/env python3


def parseTsv(sample=None, file=None, category=None):
    with open(file, 'r') as f:
        for line in f:
            l = line.strip()
            yield f'{sample}\t{category}\t{l}\n'


with open(snakemake.output[0], 'w') as out:
    for sample in snakemake.params.samples:
        for suffix in snakemake.params.suffixes.keys():
            file = os.path.join(snakemake.params.dir, sample + snakemake.params.suffixes[suffix])
            for line in parseTsv(file=file, category=suffix, sample=sample):
                out.write(line)
