# Jinliang Yang
# extract reads from bam file
# updated: 8/29/2014

getBamFiles <- function(mypath="~/NGS/BD/7-31-14/rep2"){
  
  files <- list.files(path = mypath , pattern="^Sample")
  #cat(paste("#cbind fq and trim bd RNA-seq", Sys.time(), sep=" "),
  #    file=shfile, sep="\n")
  out <- c()
  for(i in 1:length(files)){
    tmppath <- paste(mypath, files[i], sep="/")
    bam <- list.files(path=tmppath, pattern=".*-.*bam$")
    if(length(bam) > 1){
      print(tmppath)
      stop("more than one bam files found!")
    }
    out <- c(out, paste(tmppath, bam, sep="/"))
  }
  return(out)
}
 
###### => rawreads
rep1 <- getBamFiles(mypath="~/NGS/BD/7-31-14/rep1")
rep1 <- rep1[-13:-14]
rep2 <- getBamFiles(mypath="~/NGS/BD/7-31-14/rep2")
rep2 <- rep2[-13:-14]
bamfiles <- c(rep1, rep2)
bedfile <- "~/NGS/BD/7-31-14/CBF3/cbf1.txt"

source("~/Documents/BDproj/lib/extractReads.R")
extractReads(pwd="~/NGS/BD/7-31-14/CBF3", sh="extractCBF3.sh", 
             bedfile=bedfile, bamfiles=bamfiles)

###### => reads-300
source("~/Documents/BDproj/lib/extractReads.R")
bam1 <- c(rep1[1:6], rep2[1:6])
bedfile <- "~/NGS/BD/7-31-14/CBF3/cbf1-300bp.txt"
extractReads(pwd="~/NGS/BD/7-31-14/CBF3", sh="extractCBF3.sh", 
             bedfile=bedfile, bamfiles=bamfiles)

