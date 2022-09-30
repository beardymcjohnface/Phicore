#!/usr/bin/env python3

import sys

with open(snakemake.output[0],'w') as outfh:
    with open(snakemake.input[0], 'r') as infh:
        for line in infh:
            if line.startswith('>'):
                l = line.strip().replace('>','').split('_')
                seq = infh.readline().replace('*','').strip()
                outfh.write('\t'.join([
                    l[0],                       # contig
                    l[1],                       # start
                    str(int(l[1]) + len(seq)),  # stop
                    l[3],                       # orfmID
                    l[2]                        # frame
                ]))
            else:
                sys.stderr.write(f'ERROR: expecting fasta id, got {line}')
                sys.exit(1)