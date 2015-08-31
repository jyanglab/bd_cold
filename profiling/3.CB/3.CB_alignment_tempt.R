# Jinliang Yang
# updated: 5.9.2014

### the only thing you need to do is to specify the fasta file
myfasta <- "~/NGS/BD/JJ-RNA-seq/7.8.2014_more_genes/BR_signaling_gene.txt"
mydir <- "~/NGS/BD/JJ-RNA-seq/7.8.2014_more_genes"


source("~/Documents/BDproj/lib/setup_pipe.R")
########################
fq1_file <- "~/NGS/BD/7-31-14/rep2/Sample_Bri-57-1/Bri-57-1_GTCCGC.trimmed.fq";
fq2_file <- "~/NGS/BD/7-31-14/rep2/Sample_CB/CB_CCGTCC.trimmed.fq";

setup_pipe(sh="my_program.sh", database="CB_gene_gsnap", 
           fasta=myfasta,
           fq1=fq1_file, fq2=fq2_file,
           dir=mydir,
           align_output1="bri", align_output2="cb")
