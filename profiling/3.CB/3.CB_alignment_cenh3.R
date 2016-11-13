# Jinliang Yang
# updated: 5.9.2014

### the only thing you need to do is to specify the fasta file
myfasta <- "data/cenh3.fasta"
mydir <- "largedata"


source("lib/setup_pipe.R")
########################
fq1_file <- "~/NGS/BD/7-31-14/rep2/Sample_Bri-57-1/Bri-57-1_GTCCGC.trimmed.fq";
fq2_file <- "~/dbcenter/CB/CB_CCGTCC_L006_R1_001.fastq";

setup_pipe(sh="slurm-script/my_program.sh", database="CB_gene_gsnap", 
           fasta=myfasta,
           fq1=fq1_file, fq2=fq2_file,
           dir=mydir,
           align_output1="bri", align_output2="cb")

setup_pipe(sh="my_program.sh",
           gmap_build=TRUE, dbdir="~/db", dbnm="CB_fro_gsnap",
           fasta="data/cenh3.fasta",
           pwd="",
           fq=c("~/dbcenter/CB/CB_CCGTCC_L006_R1_001.fastq"))