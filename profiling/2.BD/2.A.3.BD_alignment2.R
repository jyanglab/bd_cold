# Jinliang Yang
# updated: 5.9.2014


########################################################################################
set_alignment2 <- function(dir="/mnt/02/yangjl/NGS/BD/5-14-14/", cpu=8,
                           shfile="alignment.sh"){
  
  files <- list.files(path = dir, pattern="^Sample")
  cat(paste("#alignment of 1st and 2nd run separatedly", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    fastq.gz <- list.files(path=tmppath, pattern="fastq.gz$")
    if(length(fastq.gz) != 2){
      stop(paste(tempath, "do not have two fastq.gz files!"))
    }
    
    prefix <- gsub("\\.fastq.gz$", "", fastq.gz)
    uniq <- paste(prefix, "unpaired_uniq", sep=".")
    bam <- paste(prefix, "unpaired_uniq.bam", sep=".")
    
    cat(paste("cd", files[i]),
        
        #### GSNAP alignment
        paste("#gsnap -D ~/db/brachy1.0 -d brachy1.0_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t", cpu,
              "-n 3 --quality-protocol=sanger --nofails --gunzip",
              fastq.gz[1], "--split-output", prefix[1], sep=" "),
        ### extract the unique (or reliable) aligned reads
        #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
        paste("samtools view -bS", uniq[1], ">", bam[1], sep=" "),
        paste("#gsnap -D ~/db/brachy1.0 -d brachy1.0_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t", cpu,
              "-n 3 --quality-protocol=sanger --nofails --gunzip",
              fastq.gz[2], "--split-output", prefix[2], sep=" "),
        ### extract the unique (or reliable) aligned reads
        #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
        paste("samtools view -bS", uniq[2], ">", bam[2], sep=" "),
        paste("cd .."),
        "",
        file=shfile, sep="\n", append=TRUE)
  }
  cat(paste("python ~/bin/send_email.py -s 'Alignment job finished' "),
      file=shfile, sep="\n", append=TRUE)s
  print(paste("Go to: ", dir, " | RUN: ", shfile, sep=""))  
}

### http://research-pub.gene.com/gmap/src/README
set_alignment2(dir="/mnt/02/yangjl/NGS/BD/5-14-14/", cpu=8,
               shfile="align2.sh")


align_stat2 <- function(dir="/mnt/02/yangjl/NGS/BD/5-14-14/"){
  library(ShortRead); library(Rsamtools)
  files <- list.files(path = dir, pattern="^Sample")
  
  out <- data.frame()
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    bam <- list.files(path=tmppath, pattern="uniq.bam$")
    if(length(bam)!=2){
      stop(paste("I do not seed two bam files! in", tmppath))
    }
    
    bfl1 <- BamFile(paste(tmppath, bam[1], sep="/"), yieldSize=50000,
                    index=character())
    Nalign1 <- countBam(bfl1)
    bfl2 <- BamFile(paste(tmppath, bam[2], sep="/"), yieldSize=50000,
                    index=character())
    Nalign2 <- countBam(bfl2)
    tem <- data.frame(File=c(Nalign1$file, Nalign2$file), 
                      Nalign=c(Nalign1$records,Nalign2$records), 
                      NT=c(Nalign1$nucleotides, Nalign2$nucleotides))
    out <- rbind(out, tem)
  }
  return(out)
}


###########################################################
fq_wc <- function(dir="/mnt/02/yangjl/NGS/BD/5-14-14/", shfile="fastq_wc.sh"){
  
  files <- list.files(path = dir, pattern="^Sample")
  cat(paste("#collect raw reads info of 1st and 2nd separatedly", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  outputfile <- paste(dir, "fastq_wc.txt", sep="")
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    fq.gz <- list.files(path=tmppath, pattern="fastq.gz$")
    if(length(fq.gz)!=2){
      stop(paste("I do not seed two bam files! in", tmppath))
    }
    fq <- gsub("\\.gz", "", fq.gz)
    cat(paste("cd", files[i]),
        
        #### GSNAP alignment
        paste("gunzip -cd", fq.gz[1], ">", fq[1]),
        paste("wc -l", fq[1], ">>", outputfile),
        paste("gunzip -cd", fq.gz[2], ">", fq[2]),
        paste("wc -l", fq[2], ">>", outputfile),
        paste("rm", fq[1], fq[2]),
        paste("cd .."),
        "",
        file=shfile, sep="\n", append=TRUE)  
  }
  cat(paste("python ~/bin/send_email.py -s 'fastq wc finished' "),
      file=shfile, sep="\n", append=TRUE)
  print(paste("Go to: ", dir, " | RUN: ", shfile, sep=""))   
}

########## collect info:
main <- function(){
  fq_wc(dir="/mnt/02/yangjl/NGS/BD/5-14-14/", shfile="fastq_wc.sh")
  reads <- read.table("~/NGS/BD/5-14-14/fastq_wc.txt", header=FALSE)
  names(reads) <- c("reads", "file")
  reads$file <- gsub("\\.fastq", "", reads$file)
  
  align_table <- align_stat2(dir="/mnt/02/yangjl/NGS/BD/5-14-14/")
  align_table$File <- gsub("\\.unpaired_uniq.bam", "", align_table$File)
  
  align_table2 <- merge(reads, align_table, by.x="file", by.y="File")
  align_table2$PercAligned <- round(align_table2$Nalign*4 / align_table2$reads,3)
  write.table(align_table2, "~/NGS/BD/BDproj/reports/align_table.csv",
              sep=",", row.names=FALSE, quote=FALSE)
}
#################
main()