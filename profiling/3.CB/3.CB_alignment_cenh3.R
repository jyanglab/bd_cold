# Jinliang Yang
# updated: 5.9.2014

### the only thing you need to do is to specify the fasta file


source("lib/setup_pipe.R")
########################

setup_pipe(sh="slurm-script/my_program.sh",
           gmap_build=TRUE, dbdir="largedata/", dbnm="CB_fro_gsnap",
           fasta="data/cenh3.fasta",
           pwd=".",
           fq=c("~/dbcenter/CB/CB_CCGTCC_L006_R1_001.fastq"))

# module load gmap/2015-06-23