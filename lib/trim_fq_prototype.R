###############################
trim_fq <- function(dir="/mnt/02/yangjl/NGS/BD/4-29-14", cpu=4,
                    shfile="trimming.sh"){
  
  files <- list.files(path = dir, pattern="^Sample")
  cat(paste("#cbind fq and trim bd RNA-seq", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    fq.gz <- list.files(path=tmppath, pattern="fastq.gz$")
    
    
    if(length(fq.gz)==2){
      #decompress
      prefix <- gsub("_L0.*", "", fq.gz[1])
      fq2 <- gsub("\\.gz", "", fq.gz)    
      
      fq <- paste(prefix, "fastq", sep=".")
      trimmed.fq <- paste(prefix, "trimmed.fq", sep=".")
      log <- paste(prefix, "log", sep=".")
      out <- paste(prefix, "stat", sep=".")
      
      cat(paste("cd", files[i]),
          
          #### merge two fastq files
          paste("gzip -cd", fq.gz[1], ">",  fq2[1], sep=" "),
          paste("gzip -cd", fq.gz[2], ">", fq2[2], sep=" "),
          paste("cat", fq2[1], fq2[2], ">", fq, sep=" "),
          
          # Trimming
          paste("trim_fastq.pl -i", fq, "-it fastq-sanger -ot fastq-sanger",
                "-o", trimmed.fq, "-l", log, "-p", cpu, sep=" "),
          paste("rm -r", fq2[1], fq2[2], sep=" "),
          paste("trim_stats.pl", "-l", log, "-o", out),
          paste("cd .."),
          "",
          file=shfile, sep="\n", append=TRUE
      )  
    }else{
      stop(paste(tmppath, "donot have two fq.gz files!"))
    }
  }
  cat(paste("python ~/bin/send_email.py -s 'trimming job finished' "),
      file=shfile, sep="\n", append=TRUE)
  print(paste("Go to: ", dir, " | RUN: ", shfile, sep=""))
}
