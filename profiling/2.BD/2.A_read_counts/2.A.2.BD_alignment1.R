# Jinliang Yang
# updated: 5.9.2014


set_alignment <- function(dir="/mnt/02/yangjl/NGS/BD/5-14-14/", cpu=8,
                          shfile="alignment.sh"){
  
  setwd(dir)
  files <- list.files(path = dir, pattern="^Sample")
  cat(paste("#cbind fq and trim bd RNA-seq", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    trimmed.fq <- list.files(path=tmppath, pattern="trimmed.fq$")
    
    prefix <- gsub(".trimmed.fq$", "", trimmed.fq)
    #sam <- paste(prefix, "unpaired_uniq", sep="")
    bam <- paste(prefix, "bam", sep=".")
  
    cat(paste("cd", files[i]),
        
        #### GSNAP alignment
        paste("gsnap -D ~/db/brachy1.0 -d brachy1.0_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t",
              cpu,
              "-n 3 --quality-protocol=sanger --nofails --split-output",
              prefix, trimmed.fq, sep=" "),
        paste("samtools view -bS ", prefix, ".unpaired_uniq > ", bam, sep=""),
        ### extract the unique (or reliable) aligned reads
        #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
        #paste("samtools view -bSq 1", sam, ">", bam, sep=" "),
        paste("cd .."),
        "",
        file=shfile, sep="\n", append=TRUE)
  }
  cat(paste("python ~/bin/send_email.py -s 'Alignment job finished' "),
      file=shfile, sep="\n", append=TRUE)
  print(paste("Go to: ", dir, " | RUN: ", shfile, sep=""))  
}

### http://research-pub.gene.com/gmap/src/README
### Setting up to build a GMAP/GSNAP database (one chromosome per FASTA entry)
#gmap_build -D ~/db/brachy1.0 -d brachy1.0_gsnap -g brachy1.0_wholegenome_unmasked.mfa  

### gsnap allow 2 mismatches every 36 bp and 3bp tails, parameters learned form eddy!
#gsnap -D ~/db/brachy1.0 -d brachy1.0_gsnap -m 10 -i 2 -N 1 -w 10000 -t 8 -n 3 -A sam --quality-protocol=sanger --nofails --gunzip 3-2-12-1_CAGATC_L006_R1_001.fastq.gz > 3-2-12-1_CAGATC_L006_R1_001.sam

set_alignment(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep1/", cpu=8,
              shfile="alignment.sh")

set_alignment(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep2/", cpu=8,
              shfile="alignment.sh")



align_stat <- function(dir="/mnt/02/yangjl/NGS/BD/5-14-14/"){
  library(ShortRead); library(Rsamtools)
  files <- list.files(path = dir, pattern="^Sample")
  
  out <- data.frame()
  for(i in 1:length(files)){
    tmppath <- paste(dir, files[i], sep="/")
    bam <- list.files(path=tmppath, pattern="bam$")
    
    bfl <- BamFile(paste(tmppath, bam, sep="/"), yieldSize=50000,
                   index=character())
    Nalign <- countBam(bfl)
    tem <- data.frame(File=Nalign$file, Nalign=Nalign$records, NT=Nalign$nucleotides)
    out <- rbind(out, tem)
    }
  return(out)
}

##################
align_table1 <- align_stat(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep1")
write.table(align_table1, "~/Documents/BDproj/reports/align_table_rep1.csv",
            sep=",", row.names=FALSE, quote=FALSE)


align_table <- align_stat(dir="/mnt/02/yangjl/NGS/BD/7-31-14/rep2")
write.table(align_table, "~/Documents/BDproj/reports/align_table_rep2.csv",
            sep=",", row.names=FALSE, quote=FALSE)



