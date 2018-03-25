from os.path import join
import glob

configfile: "./config.yaml"

# globals
REF = config['genome']
INDEX = config['index']
GTF = config['gtf']
ADAPTORS = config['adaptors']
DIRS = ['bams/', 'raw_reads/', 'clean_reads/', 'logs/', 'counts/']

# key step to get sample names from R1 read, example: AG0069-01_R1_001.fastq.gz
SAMPLES, = glob_wildcards(join('raw_reads',
    '{samples,AG0069[^/]+}_R1_001.fastq.gz'))
print(SAMPLES)

PATTERN_R1 = config['pat_r1']
PATTERN_CLN_R1 = config['cln_r1']
print(PATTERN_R1)
rule all:
    input:
        expand('clean_reads/{sample}_R1_001.cln.fastq.gz', sample=SAMPLES),
        REF, GTF, ADAPTORS, DIRS,
        expand('{INDEX}.1.ht2', INDEX=INDEX),
        expand('bams/{sample}.sam', sample=SAMPLES),
        expand('bams/{sample}.sbn.bam', sample=SAMPLES),
        expand('counts/{sample}.sbn.counts', sample=SAMPLES)

rule project_setup:
    output: DIRS
    shell:
        "mkdir -p "+' '.join(DIRS)

rule make_index:
    input:
        expand('{REF}',REF=REF)
    output:
        expand('{INDEX}.1.ht2', INDEX=INDEX),
    threads: 12
    log:
        'logs/index_log.txt'
    shell:
        'hisat2-build -p {threads} {input} $(echo "{input}" | sed -e "s/.fa//") > {log} 2>&1'

rule qc_trim:
    input:
        r1 = join('raw_reads', PATTERN_R1),
    threads: 12
    params:
        minlen=config['minlen'],
        qtrim=config['qtrim'],
        trimq=config['trimq'],
        ktrim=config['ktrim'],
        kwin=config['kwin'],
        mink=config['mink'],
        hdist=config['hdist']
    output:
        r1_out = 'clean_reads/{sample}_R1_001.cln.fastq.gz',
    log:
        'logs/trim_log.txt'
    shell:
        "bbduk.sh in={input.r1} out={output.r1_out} "
        "minlen={params.minlen} qtrim={params.qtrim} trimq={params.trimq} "
        "ktrim={params.ktrim} k={params.kwin} mink={params.mink} "
        "ref={ADAPTORS} hdist={params.hdist} 2>&1 | tee -a {log}"

# not loggin output, need to check params
rule aln:
    input:
        ref=expand('{REF}',REF=REF),
        r1 = join('clean_reads', PATTERN_CLN_R1),
    log:
        'logs/aln_log.txt'
    output:
        sam='bams/{sample}.sam',
        splice='ref/{sample}.novel_splices.txt'
    params:
        index=expand('{INDEX}',INDEX=INDEX),
        strand=config['alnstrand']
    threads: 24
    shell:
        """
        echo {output.sam} >> {log}
        hisat2 -p {threads} -x {params.index} -q --dta -U {input.r1} --rna-strandness {params.strand} --novel-splicesite-outfile {output.splice} -S {output.sam} >> {log} 2>&1
        """

rule sam_to_bam:
    input:
        'bams/{sample}.sam'
    output:
        'bams/{sample}.sbn.bam'
    shell:
        "samtools view -bS {input} | samtools sort -n - bams/{wildcards.sample}.sbn"

# does not work on name sorted files, bugger
#rule aln_index:
#    input:
#        'bams/{sample}.sbn.bam'
#    output:
#        'bams/{sample}.sbn.bam.bai'
#    shell:
#        "samtools index {input}"

rule do_counts:
    input:
        'bams/{sample}.sbn.bam'
    output:
        'counts/{sample}.sbn.counts'
    log:
        'logs/{sample}_count.log'
    params:
        strand=config['countstrand'],
        countmode=config['mode']
    shell:
        'htseq-count -r name -s {params.strand} -f bam -m {params.countmode} {input} '
        '{GTF} > {output} 2> {log}'
