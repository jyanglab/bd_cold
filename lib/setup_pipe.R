# Jinliang Yang
# a function to setup db and align NGS to it.
# updated: version2 8/28/2014

setup_pipe <- function(sh="my_program.sh",
                       gmap_build=TRUE, dbdir="~/db", dbnm="CB_fro_gsnap",
                       fasta="~/DBcenter/CB_FROs/fro_arabidopsis.fasta",
                       pwd="",
                       fq=c("Bri-57-1_GTCCGC.trimmed.fq", "CB_CCGTCC.trimmed.fq")){
  
  setwd(pwd)
  #### setup database
  #gmap_build -D ~/db -d CB_Bri1_gsnap ~/DBcenter/CB_Bri1/Bri1_gene.fasta
  if(gmap_build == TRUE){
    cat(paste("gmap_build -D", dbdir, "-d", dbnm, fasta),
        file=sh, sep="\n");
  }
  
  #### alignment
  ### gsnap allow 2 mismatches every 36 bp and 3bp tails, parameters learned form eddy!
  for(i in 1:length(fq)){
    idbase <- gsub("\\..*", "", fq[i])
    cat(
      paste("gsnap -D", dbdir, "-d", dbnm, 
            "-m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output",
            idbase, fq[i], sep=" "),
      "",
      file=sh, sep="\n", append=TRUE)
  }
    
  ### samtools
  for(i in 1:length(fq)){
    idbase <- gsub("\\..*", "", fq[i])
    uniq <- paste(idbase, ".unpaired_uniq", sep="")
    bam <- paste(idbase, ".bam", sep="")
    
    cat(
      #CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
      # samtools view -bS Bri_pro.unpaired_uniq > Bri_pro_uniq.bam
      paste("samtools view -bS", uniq, ">", bam, sep=" "),
      #CREATING A BAM INDEX FILE
      #samtools sort Bri_pro_uniq.bam Bri_pro_uniq.sorted
      paste("samtools sort ", bam, " ", idbase, ".sorted", sep=""),
      #samtools index Bri_pro_uniq.sorted.bam Bri_pro_uniq.sorted.bam.bai
      paste("samtools index ", idbase, ".sorted.bam ", idbase, ".sorted.bam.bai", sep=""),
      "",
      file=sh, sep="\n", append=TRUE);
  }
  message(sprintf("go to [cd %s ] and run $sh %s", pwd, sh))
}

