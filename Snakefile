from os.path import join
import glob

configfile: "./config.yaml"

# globals
REF = config['genome']
THREADS = config['threads']
INDEX = config['index']
GTF = config['gtf']
ADAPTORS = config['adaptors']
DIRS = ['bams/', 'raw_reads/', 'clean_reads/', 'logs/', 'counts/', 'ref/', 'tables']

# key step to get sample names from R1 read
SAMPLES, = glob_wildcards('raw_reads/{sample}_R1_001.fastq.gz')
#CC9ARANXX-3492-01_R1_001.fastq.gz
PATTERN_R1 = config['pat_r1']
PATTERN_R2 = config['pat_r2']
PATTERN_CLN_R1 = config['cln_r1']
PATTERN_CLN_R2 = config['cln_r2']

print('using the following sample wildcards')
print(SAMPLES)

rule all:
    input:
        expand('clean_reads/{sample}_R1_001.cln.fastq.gz', sample=SAMPLES),
        REF, GTF, ADAPTORS, DIRS,
        expand('{INDEX}.1.ht2', INDEX=INDEX),
        expand('bams/{sample}.bam', sample=SAMPLES),
        expand('bams/{sample}.sbn.bam', sample=SAMPLES),
        expand('counts/{sample}.sbn.counts', sample=SAMPLES),
        'logs/count_results.log'

rule project_setup:
    output: DIRS
    shell:
        "mkdir -p "+' '.join(DIRS)

rule make_index:
    input:
        expand('{REF}',REF=REF)
    output:
        expand('{INDEX}.1.ht2', INDEX=INDEX),
    threads: THREADS
    log:
        'logs/index_log.txt'
    shell:
        'hisat2-build -p {threads} {input} $(echo "{input}" | sed -e "s/.fa//") > {log} 2>&1'

rule qc_trim:
    input:
        r1 = join('raw_reads', PATTERN_R1),
        r2 = join('raw_reads', PATTERN_R2)
    threads: THREADS
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
        r2_out = 'clean_reads/{sample}_R2_001.cln.fastq.gz'
    log:
        'logs/trim_log.txt'
    shell:
        "/opt/bbmap/bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1_out} out2={output.r2_out} "
        "minlen={params.minlen} qtrim={params.qtrim} trimq={params.trimq} "
        "ktrim={params.ktrim} k={params.kwin} mink={params.mink} "
        "ref={ADAPTORS} hdist={params.hdist} 2>&1 | tee -a {log}"

# not loggin output, need to check params
rule aln:
    input:
        ref=expand('{REF}',REF=REF),
        r1 = join('clean_reads', PATTERN_CLN_R1),
        r2 = join('clean_reads', PATTERN_CLN_R2)
    log:
        'logs/aln_log.txt'
    output:
        bam='bams/{sample}.bam',
        splice='ref/{sample}.novel_splices.txt'
    params:
        index=expand('{INDEX}',INDEX=INDEX),
        strand=config['alnstrand']
    threads: THREADS
    shell:
        """
        echo {output.bam} >> {log}
        hisat2 -p {threads} -x {params.index} -q --dta -1 {input.r1} -2 {input.r2} --rna-strandness {params.strand} --novel-splicesite-outfile {output.splice} 2>> {log} | samtools view -Sb - > {output.bam}
        """
# don't delete the pos sorted bam or rerunning will try to regenerate it (ie
# aln again...) 
#### WARNING BROKEN ####
rule sam_to_bam:
    input:
        'bams/{sample}.bam'
    output:
        'bams/{sample}.sbn.bam'
    shell:
        #"samtools sort -@ 4 -m 20G -n {input} -o {output}"
        # think about filtering non-uniuq mapped reads here
        # but see discussion here on not removing them: https://www.biostars.org/p/55648/
        # and a dicusssion on removing them: http://www.bioconductor.org/help/course-materials/2015/CSAMA2015/lect/L05-deseq2-anders.pdf
        # MY typical workflow is not too remove them unless lots in the data QC

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

rule log_count_result:
    input:
        expand('counts/{sample}.sbn.counts', sample=SAMPLES)
    output:
        'logs/count_results.log',
    shell:
        'tail -n5 counts/*counts > {output}'

rule make_latex_tables:
    input:
        aln='logs/aln_log.txt',
        qc='logs/trim_log.txt'
    output:
        aln='tables/aln_table.tex',
        qc='tables/qc_table.tex'
    shell:
        """
        python scripts/make_aln_tab.py {input.aln} > {output.aln}
        python scripts/make_qc_tab.py {input.qc} > {output.qc}
        """

rule clean:
    shell:
        """
        rm -f bams/*
        rm -f tables/*
        rm -f counts/*
        rm -f logs/*
        rm -f clean_reads/*
        rm -f ref/*novel_splices.txt
        """
