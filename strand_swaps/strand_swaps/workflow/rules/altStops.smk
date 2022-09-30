

rule prot_to_bed:
    input:
        os.path.join(outDir,'{file}.prot')
    output:
        os.path.join(outDir, '{file}.prot.bed')
    script:
        '../scripts/prot2bed.py'

rule genome_window:
    input:
        os.path.join(inDir, '{file}.gbk')
    output:
        temp(os.path.join(outDir, '{file}.genome.bed'))
    shell:
        """grep LOCUS {input} | awk '{{print $2"\t0\t"$3}}' > {output}"""

rule prot_gaps:
    input:
        prot = os.path.join(outDir, '{file}.prot.bed'),
        gen  = os.path.join(outDir, '{file}.genome.bed')
    output:
        os.path.join(outDir, '{file}.gaps.bed')
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

rule orfm_bed:
    """Get BED coords for ORFM fragments"""
    input:
        '{file}.orfm.faa'
    output:
        '{file}.orfm.bed'
    script:
        '../scripts/orfm2bed.py'

rule orfm_gap_filt:
    """Filter orfm fragments that overlap protein annotations"""
    input:
        orf = '{file}.orfm.bed',
        gap = '{file}.gaps.bed'
    output:
        '{file}.gaporfm.bed'
    conda:
        '../envs/bedtools.yaml'
    shell:
        "bedtools -a {input.orf} -b {input.gap} -wa > {output}"