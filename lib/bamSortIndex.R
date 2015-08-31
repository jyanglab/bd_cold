# Jinliang Yang
# bam sort and index
# updated: 8/29/2014

#############
bamSortIndex <- function(pwd="~", shfile="my_program.sh", mypattern=".*-.*bam$",
                         mypath="~"){
  
  setwd(pwd)
  
  files <- list.files(path = mypath , pattern="^Sample")
  cat(paste("#bam sort and index", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  for(i in 1:length(files)){
    tmppath <- paste(mypath, files[i], sep="/")
    bam <- list.files(path=tmppath, pattern=mypattern)
    bam <- paste(tmppath, bam, sep="/")
    if(length(bam) > 1){
      print(tmppath)
      stop("more than one bam files found!")
    }else{
      idbase <- bam
      cat(
        #samtools sort reads_in_cbf1.bam reads_in_cbf1.sorted
        paste("samtools sort ", bam, " ", idbase, ".sort", sep=""),
        #samtools index Bri_pro_uniq.sorted.bam Bri_pro_uniq.sorted.bam.bai
        paste("samtools index ", idbase, ".sort.bam ", idbase, ".sort.bam.bai", sep=""),
        "",
        file=shfile, sep="\n", append=TRUE);
    } 
  }
  message(sprintf("go to [cd %s ] and run $sh %s", pwd, shfile))
}

