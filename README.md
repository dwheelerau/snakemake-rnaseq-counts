# snakemake-rnaseq-counts on single end RNA-seq data   
This snakemake workflow will run a hisat2 - htseq count based gene count
workflow. Trimming and filtering is also carried out using bbduk.  

## Details  
Open the `config.yaml` file in a text editor and change the params as required.  

Need to create the directory structure.  
Run ```snakemake project_setup```  

Need to create an index BEFORE running all.  
Run ```snakemake make_index --cores 12 > logs/snakemake_run.log```  

Run the script ```snakemake all --cores 12 >> logs/snakemake_run.log```  

## Requirements   
1.  The raw sequenced reads need to be found in a directory called `raw_reads`, the
script expects the reads to be fastq.gz files, though it may work with standard
fastq files, however these will be converted to gz files after the cleanup
script. The file pattern needs to be set in the config file to recognise the
sample names that should be incorporated into the files.   

2.  The file `contams_forward_rev.fa` contains adaptors etc that will be
removed during the cleaning process.  

3.  Run the script using ```nohup snakemake all --cores 24 > snakemake_run.log
    &```, using as many cores are you have available.  

4.  The final output will be ht-seq based counts for each sample. These can be
    used as raw reads for DESeq2.  

5.  Latex tables can be generated for the qc and alnment data see rules below.   

## Rules  
rule all - run the entire pipeline, run ```project_setup``` and ```make_index``` first.  
rule project_setup - setup the diretory structure  
rule make_index - make hisat2 index of reference genome, genome must be .fa  
rule qc_trim - run bbduk.sh script to quality filter reads  
rule aln - aln qc reads to reference using hisat2
rule sam_to_bam - convert sam to bam using samtools1.7  
rule do_counts - geneate count files using htseq-count  
rule log_count_result - output final 5 lines of count files  
rule make_latex_tables - make latex tables from key log files for aln and qc  
rule clean - delete everything except ```raw_reads```  
