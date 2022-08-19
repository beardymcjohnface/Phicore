#!/usr/bin/env python3

import re

with open(snakemake.input[0], 'r') as f:
    with open(snakemake.output[0], 'w') as o:
        for line in f:
            l = line.strip().split('\t')
            l[0] = re.sub('<|>', '', l[0])
            coords = l[0].split('..')
            ctg = re.sub('\[.*', '', l[5])
            o.write(f'{ctg}\t{coords[0]}\t{coords[1]}\t{l[1]}\n')
