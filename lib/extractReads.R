# Jinliang Yang
# extract reads from bam file
# updated: 8/29/2014

#############
extractReads <- function(pwd="", sh="my_program.sh", 
                         bedfile="~/DBcenter/CB_FROs/fro_arabidopsis.fasta",
                         bamfiles=c("Bri-57-1_GTCCGC.trimmed.fq", "CB_CCGTCC.trimmed.fq")){
  
  setwd(pwd)   
  ### use bedtools and samtools
  for(i in 1:length(bamfiles)){
    
    idbase <- gsub("\\..*", "", bamfiles[i])
    idbase <- gsub(".*\\/", "", idbase)
    idbase <- paste(i, idbase, sep="_")
    bam <- paste(idbase, ".bam", sep="")
    
    cat(
      # bedtools intersect -abam 3-2-12-1_CAGATC.bam -b cbf1.txt -wa > reads_in_cbf1.bam
      paste("bedtools intersect -abam", bamfiles[i], "-b", bedfile,  "-wa >", bam, sep=" "),
      #samtools sort reads_in_cbf1.bam reads_in_cbf1.sorted
      paste("samtools sort ", bam, " ", idbase, ".sort", sep=""),
      #samtools index Bri_pro_uniq.sorted.bam Bri_pro_uniq.sorted.bam.bai
      paste("samtools index ", idbase, ".sort.bam ", idbase, ".sort.bam.bai", sep=""),
      "",
      file=sh, sep="\n", append=TRUE);
  }
  message(sprintf("go to [cd %s ] and run $sh %s", pwd, sh))
}

