# Jinliang Yang
# 5/15/2014
# trimming and collect trimming results
# to do: change to trimmatic

################
uncompress <- function(dir="/mnt/02/yangjl/NGS/BD/7-31-14", fileptn="tar$", 
                       command="cat -xvf", shfile="uncompress.sh"){
  files <- list.files(path = dir, pattern=fileptn)
  ## file line of the shell file:
  setwd(dir)
  cat(paste("#uncompress", Sys.time(), sep=" "),
      file=shfile, sep="\n"
  )
  for(i in 1:length(files)){
    cat(paste(command, files[i], sep=" "),
        file=shfile, sep="\n", append=TRUE
    )
  }
  print(paste("Go to: cd ", dir, "| RUN: sh ", shfile, sep=""))
}
##########
uncompress(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep1", fileptn="tar$", 
           command="tar -xvf", shfile="uncompress.sh")
uncompress(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep2", fileptn="tar$", 
           command="tar -xvf", shfile="uncompress.sh")

###############################
trim_fq <- function(dir="/mnt/02/yangjl/NGS/BD/4-29-14", cpu=4,
                    shfile="trimming.sh"){
  
  setwd(dir)
  files <- list.files(path = dir, pattern="^Sample")
  cat(paste("#cbind fq and trim bd RNA-seq", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    fq.gz <- list.files(path=tmppath, pattern="fastq.gz$")
    
    
    if(length(fq.gz)==1){
      #decompress
      prefix <- gsub("_L0.*", "", fq.gz)
      #fq2 <- gsub("\\.gz", "", fq.gz)    
      
      fq <- paste(prefix, "fastq", sep=".")
      trimmed.fq <- paste(prefix, "trimmed.fq", sep=".")
      log <- paste(prefix, "log", sep=".")
      out <- paste(prefix, "stat", sep=".")
      
      cat(paste("cd", files[i]),
          
          #### merge two fastq files
          paste("gzip -cd", fq.gz, ">",  fq, sep=" "),
          #paste("gzip -cd", fq.gz[2], ">", fq2[2], sep=" "),
          #paste("cat", fq2[1], fq2[2], ">", fq, sep=" "),
          
          # Trimming
          paste("trim_fastq.pl -i", fq, "-it fastq-sanger -ot fastq-sanger",
                "-o", trimmed.fq, "-l", log, "-p", cpu, sep=" "),
          paste("rm -r", fq, sep=" "),
          paste("trim_stats.pl", "-l", log, "-o", out),
          paste("cd .."),
          "",
          file=shfile, sep="\n", append=TRUE
      )  
    }else{
      stop(paste(tmppath, "donot have fq.gz files!"))
    }
  }
  cat(paste("python ~/bin/send_email.py -s 'trimming job finished' "),
      file=shfile, sep="\n", append=TRUE)
  print(paste("Go to: cd ", dir, " | RUN: sh ", shfile, sep=""))
}

#########
trim_fq(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep1", cpu=2, shfile="trimming.sh")
trim_fq(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep2", cpu=2, shfile="trimming.sh")


###############################
collect_trim_stat <- function(dir="/mnt/02/yangjl/NGS/BD/4-29-14"){
  files <- list.files(path = dir, pattern="^Sample")
  res <- data.frame()
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    statfile <- list.files(path=tmppath, pattern="stat$")
    if(length(statfile) != 1){
      #return(res)
      stop(paste("Error!!!", tmppath))
    }else{
      statfile2 <- paste(tmppath, statfile, sep="/")
      stat <- read.delim(statfile2, header=TRUE)
      names(stat) <- c("File", "NoRawReads", "RawBasePairs", "ReadLength", "NoTrimmedReads", 
                       "TrimmedReadsRemaining", "TrimmedBasePairs", "TrimmedBasesRemaining",
                       "ReadLength", "NoTrimmedReadsP", "TrimmedBasePairsP","ReadLength2"
                       )
      res <- rbind(res, stat)
    }
  }
  return(res)  
}

#####
res1 <- collect_trim_stat(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep1")
write.table(res1, "~/Documents/BDproj/reports/bd_rep1_report.txt", sep="\t",
            row.names=FALSE, quote=FALSE)
#####
res2 <- collect_trim_stat(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep2")
write.table(res2, "~/Documents/BDproj/reports/bd_rep2_report.txt", sep="\t",
            row.names=FALSE, quote=FALSE)