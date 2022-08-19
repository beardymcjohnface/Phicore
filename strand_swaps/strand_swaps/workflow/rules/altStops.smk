

rule prot_to_bed:
    input:
        os.path.join(inDir,'{file}.prot')
    output:
        os.path.join(inDir, '{file}.prot.bed')
    script:
        '../scripts/prot2bed.py'

rule genome_window:
    input:
        os.path.join(inDir, '{file}.gbk')
    output:
        temp(os.path.join(inDir, '{file}.genome.bed'))
    shell:
        """grep LOCUS {input} | awk '{{print $2"\t0\t"$3}}' > {output}"""

rule prot_gaps:
    input:
        prot = os.path.join(inDir, '{file}.prot.bed'),
        gen  = os.path.join(inDir, '{file}.genome.bed')
    output:
        os.path.join(inDir, '{file}.gaps.bed')
    conda:
        '../envs/bedtools.yaml'
    shell:
        """bedtools subtract -a {input.gen} -b {input.prot} > {output}"""

rule orfm:
    """Run ORFM on the input fasta genome file"""
    input:
        '{file}.fna'
    output:
        prot = '{file}.orfm.faa',
        nucl = '{file}.orfm.fna'
    conda:
        "../envs/orfm.yaml"
    shell:
        "orfm {input} -t {output.nucl} -p -s > {output.prot}"
