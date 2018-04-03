# snakemake-rnaseq-counts  
This snakemake workflow will run a hisat2 - htseq count based gene count
workflow. Trimming and filtering is also carried out using bbduk.  

## envs  
conda env available to run this pipeline ```conda create --name myenv --file envs/myenv.yaml```

## Details  
Open the `config.yaml` file in a text editor and change the params as required.  
You need to run ```project_setup``` and ```make_index``` before make all.   

Run the make_latex_tables rule to generate latex tables from key log files  

## Requirements   
1.  The raw sequenced reads need to be found in a directory called `raw_reads`, the
script expects the reads to be fastq.gz files, though it may work with standard
fastq files, however these will be converted to gz files after the cleanup
script.  

2.  The file `contams_forward_rev.fa` contains adaptors etc that will be
removed during the cleaning process. 

3.  All directories need to be in the sample place as the Snakemake file.  

4.  Run the script using `nohup snakemake all --cores 24 > snakemake_run.log
    &`, using as many cores are you have available.  

5.  The final output will be ht-seq based counts for each sample. These can be
    used as raw reads for DESeq2.

## Rules  
rule all - Run entire pipeline, need to have run ```project_setup``` and
```make_index```  
rule project_setup - setup directory structure   
rule make_index - create hisat2 index from reference.fa file  
rule qc_trim - bbduk.sh qc script  
rule aln - hisat2 alnment  
rule sam_to_bam  - convert sam to bam  
rule do_counts - generate read counts  
rule log_count_result - log last 5 lines of each count file  
rule make_latex_tables - generate latex tables from key log files for aln and
qc  

## Tests  
## Tests  
To test the pipeline:  
1.  Copy the tests/config.yaml to the base directory  
2.  Copy tests/*.gz to raw_reads/  
3.  Download (from ensembl) ```Danio_rerio.GRCz10.dna.toplevel.fa``` and copy
    to
    ref/  
4.  Download (from ensembl) ```ref/Danio_rerio.GRCz10.91.gtf``` and copy to
    ref/  
5.  Run ```snakemake project_setup```  
6.  Run ```snakemake make_index```  
7.  Run ```snakemake all```  
8.  Check the result:  
 ```
diff counts/CC9ARANXX-3492-01.sbn.counts tests/CC9ARANXX-3492-01.sbn.counts  
diff counts/CC9ARANXX-3492-02.sbn.counts tests/CC9ARANXX-3492-02.sbn.counts  
```
