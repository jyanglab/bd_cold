### Jinliang Yang
### 8/29/2014

### sort and index bam files before call reads count
### using the unpaired_mult as well
source("~/Documents/BDproj/lib/SamSortIndex.R")
SamSortIndex(pwd="~/NGS/BD/7-31-14/CBF3", shfile="sortbam_rep1.sh", mypattern="unpaired_mult$",
             mypath="~/NGS/BD/7-31-14/rep1")
SamSortIndex(pwd="~/NGS/BD/7-31-14/CBF3", shfile="sortbam_rep2.sh", mypattern="unpaired_mult$",
             mypath="~/NGS/BD/7-31-14/rep2")

cbf8 <- read.csv("~/NGS/BD/7-31-14/CBF3/CBF_8genes_bed.csv", header=FALSE)
write.table(cbf8, "~/NGS/BD/7-31-14/CBF3/CBF_8genes.bed", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)


#############
col_files <- function(mypattern="unpaired_mult.sort.bam$",
                      searchpath="~/NGS/BD/7-31-14/rep1"){
  
  files <- list.files(path = searchpath, pattern="^Sample")
  cat(paste("#collect read count", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  out <- c()
  for(i in 1:length(files)){
    tmppath <- paste(mypath, files[i], sep="/")
    bam <- list.files(path=tmppath, pattern=mypattern)
    bam <- paste(tmppath, bam, sep="/")
    if(length(bam) > 1){
      print(tmppath)
      print(bam)
      stop("more than one sam files found!")
    }else{
      out <- c(out, bam)
    }
  }
  return(out)
}





samtools sort ~/NGS/BD/7-31-14/rep1/Sample_3-2-12-1/3-2-12-1_CAGATC.unpaired_mult ~/NGS/BD/7-31-14/rep1/Sample_3-2-12-1/3-2-12-1_CAGATC.unpaired_mult.sort
samtools index ~/NGS/BD/7-31-14/rep1/Sample_3-2-12-1/3-2-12-1_CAGATC.unpaired_mult.sort.bam ~/NGS/BD/7-31-14/rep1/Sample_3-2-12-1/3-2-12-1_CAGATC.unpaired_mult.sort.bam.bai
#####
cd ~/NGS/BD/7-31-14/CBF3
sh count_cbf.txt

cbf1 <- read.table("~/NGS/BD/7-31-14/CBF3/cbf_rc-300_rep1.txt", header=FALSE)
cbf2 <- read.table("~/NGS/BD/7-31-14/CBF3/cbf_rc-300_rep2.txt", header=FALSE)

wt1 <- read.table("~/NGS/BD/7-31-14/CBF3/cbf_rc_rep1.txt", header=FALSE)
wt2 <- read.table("~/NGS/BD/7-31-14/CBF3/cbf_rc_rep2.txt", header=FALSE)

cbf <- cbf1
cbf[, 4:9] <- cbf1[, 4:9] +cbf2[, 4:9]
#V1       V2       V3 V4 V5 V6 V7 V8 V9
#1 Bd4 81954700 81955838 21 19 18  4  0  5

wt <- wt1
wt[, 4:9] <- wt1[, 4:9] + wt2[, 4:9]
#V1       V2       V3 V4 V5 V6  V7  V8 V9
#1 Bd4 40976869 40977919 90 45 66 279 110 13
