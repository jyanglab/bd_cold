# Jinliang Yang
# updated: 5.9.2014

### the only thing you need to do is to specify the fasta file
myfasta <- "~/NGS/BD/Bri_drought/Bri1_drought_related_genes.txt"
mydir <- "~/NGS/BD/Bri_drought/"


source("~/Documents/BDproj/lib/setup_pipe.R")
########################
Bri_rep1 <- "~/NGS/BD/7-31-14/rep1/Sample_Bri-57-1/Bri-57-1_GTCCGC.trimmed.fq";
Bri_rep2 <- "~/NGS/BD/7-31-14/rep2/Sample_Bri-57-1/Bri-57-1_GTCCGC.trimmed.fq";
WT_rep1 <- "~/NGS/BD/7-31-14/rep1/Sample_CB/CB_CCGTCC.trimmed.fq";
WT_rep2 <- "~/NGS/BD/7-31-14/rep2/Sample_CB/CB_CCGTCC.trimmed.fq";

fqfiles <- c(Bri_rep1, Bri_rep2, WT_rep1, WT_rep2)
setup_pipe(sh="drought_run.sh",
           gmap_build=TRUE, dbdir="~/db", dbnm="CB_drought_gsnap",
           fasta="~/NGS/BD/Bri_drought/Bri1_drought_related_genes.txt",
           pwd="~/NGS/BD/Bri_drought/",
           fq=fqfiles)




