# snakemake-rnaseq-counts  
This snakemake workflow will run a hisat2 - htseq count based gene count
workflow. Trimming and filtering is also carried out using bbduk.  

# Details  
Open the `config.yaml` file in a text editor and change the params as required.  

Need to create an index first.  
Run ```snakemake make_index --cores 12```  

Run the script ```snakemake all --cores 12```  

# Requirements   
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
6.  Latex tables are generated automatically for the qc and alnment data  
