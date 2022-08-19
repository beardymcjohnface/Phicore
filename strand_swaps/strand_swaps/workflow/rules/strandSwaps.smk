
rule convert_genbank:
    input:
        os.path.join(inDir,'{file}.gbk')
    output:
        prot = temp(os.path.join(inDir,'{file}.prot')),
        fna  = temp(os.path.join(inDir,'{file}.fna'))
    params:
        script = os.path.join(workflow.basedir, 'scripts', 'genbank2sequences.py'),
        prefix = os.path.join(inDir,'{file}')
    conda:
        '../envs/pystuff.yaml'
    shell:
        """python {params.script} -g {input} -p {output.prot} -n {params.prefix}"""

rule summarise_swaps:
    input:
        os.path.join(inDir,'{file}.prot')
    output:
        strnd = temp(os.path.join(inDir,'{file}.consec_strand.tsv')),
        frame = temp(os.path.join(inDir,'{file}.consec_frame.tsv')),
        ovlps = temp(os.path.join(inDir,'{file}.consec_overlap.tsv'))
    script:
        '../scripts/summarise_swaps.py'

rule summaries_to_table:
    input:
        strand = expand(os.path.join(inDir,'{sample}.consec_strand.tsv'), sample=samples),
        frame = expand(os.path.join(inDir,'{sample}.consec_frame.tsv'), sample=samples),
        overlp = expand(os.path.join(inDir,'{sample}.consec_overlap.tsv'), sample=samples)
    output:
        os.path.join(inDir, 'summary.consec.tsv')
    params:
        samples = samples,
        suffixes = {'strand':'.consec_strand.tsv','frame':'.consec_frame.tsv','overlap':'.consec_overlap.tsv'},
        dir = inDir
    script:
        '../scripts/combineConsecTsv.py'
