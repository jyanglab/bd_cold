# Jinliang Yang
# bam sort and index
# updated: 8/29/2014

#############
SamSortIndex <- function(pwd="~", shfile="my_program.sh", mypattern=".*-.*bam$",
                         mypath="~"){
  
  setwd(pwd)
  
  files <- list.files(path = mypath , pattern="^Sample")
  cat(paste("#bam sort and index", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  for(i in 1:length(files)){
    tmppath <- paste(mypath, files[i], sep="/")
    sam <- list.files(path=tmppath, pattern=mypattern)
    sam <- paste(tmppath, sam, sep="/")
    if(length(sam) > 1){
      print(tmppath)
      print(sam)
      stop("more than one sam files found!")
    }else{
      bam <- paste(sam, ".bam", sep="")
      cat(
        paste("samtools view -bS", sam, ">", bam),
        #samtools sort reads_in_cbf1.bam reads_in_cbf1.sorted
        paste("samtools sort ", bam, " ", sam, ".sort", sep=""),
        #samtools index Bri_pro_uniq.sorted.bam Bri_pro_uniq.sorted.bam.bai
        paste("samtools index ", sam, ".sort.bam ", sam, ".sort.bam.bai", sep=""),
        "",
        file=shfile, sep="\n", append=TRUE);
    } 
  }
  message(sprintf("go to [cd %s ] and run $sh %s", pwd, shfile))
}

