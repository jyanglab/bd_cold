### Jinliang Yang
### Dec 13th, 2016
### Off-target analysis


library("Biostrings")

fa <- readDNAStringSet(filepath = "data/RNAi_CBF3_Target.fa", format="fasta")
class(fa)
length(fa[[1]])
width(fa)


out <- c()
for(i in 1:(width(fa)-18)){
   tmp <-  fa[[1]][i:(i+18)]
   out <- c(out, tmp)
}






################### summarize the uniquely aligned reads!
library(ShortRead); library(Rsamtools)
Nreads <- 17613904/4

write.table(read_statsDF, "results/read_statsDF.xls", row.names=FALSE, quote=FALSE, sep="\t")


############# fastq QC report, not working!!!
source("~/Documents/BDproj/lib/fastqQuality.R")
## myfiles <- paste0("data/", targets$FileName); names(myfiles) <- targets$SampleName

myfiles <- "/mnt/02/yangjl/NGS/BD/7-31-14/rep2/Sample_3-2-12-1/3-2-12-1_CAGATC.trimmed.fq"
#names(myfiles) <- "wt"
setwd("/mnt/02/yangjl/NGS/BD/BDproj/data/")
system("wget http://biocluster.ucr.edu/~tgirke/HTML_Presentations/Manuals/Rngsapps/chipseqBioc2012/data.zip");
system("unzip data.zip")
fqlist <- seeFastq(fastq=myfiles, batchsize=50000, klength=8)
#pdf("~/NGS/BD/BDproj/reports//fastqReport.pdf", height=18, width=4*length(myfiles)); 
seeFastqPlot(fqlist); 
#dev.off()



# source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/fastqQuality.R")
library(ggplot2)
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
# names(fastq) <- paste("flowcell_6_lane", 1:4, sep="_")
# fqlist <- seeFastq(fastq=fastq, batchsize=100000, klength=8)
# pdf("fastqQuality.pdf", height=16, width=4*length(fastq))
# seeFastqPlot(fqlist, arrange=seq(along=fqlist))
# dev.off()


library(qrqc)
s.fastq <- readSeqFile(filename="/mnt/02/yangjl/NGS/BD/7-31-14/rep2/Sample_3-2-12-1/3-2-12-1_CAGATC.trimmed.fq",
                       max.length=110)
trimed <- readSeqFile(filename="/mnt/02/yangjl/NGS/BD/4-29-14/Sample_3-2-12-1/3-2-12-1_CAGATC_L006_R1_001.trimmed.fq",
                      max.length=150)
qualPlot(s.fastq)
qualPlot(list("trimmed"=trimed, "untrimmed"=s.fastq))



readCountsPerGene.pl -g ~/DBcenter/BD_v1.0/annot_v1.2/brachypodium1.2_20100223_MIPSGFF.gff -o 3-2-12-1_CAGATC_L006_R1_001.count  --gff3 3-2-12-1_CAGATC_L006_R1_001.gff3

setup_gsnap <- function(){
  
}




### modify the gff file, nine columns
bd1 <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/brachypodium1.2_20100223_MIPSGFF.gff", header=FALSE)
bd1$V3 <- gsub("_predicted", "", bd1$V3)
bd1$V3 <- gsub("predicted_", "", bd1$V3)
#CDS   exon   gene   mRNA 
#167295 173325  26552  31029
write.table(bd1, "~/DBcenter/BD_v1.0/annot_v1.2/brachypodium1.2_modified.gff", sep="\t",col.names=FALSE,
            row.names=FALSE, quote=FALSE)

samplefile <- system.file("extdata", "UCSC_knownGene_sample.sqlite", package = "GenomicFeatures")
txdb <- loadDb(samplefile)
txdb
###################################################
### prepare genomic features
###################################################
library(GenomicFeatures)
txdb <- makeTranscriptDbFromGFF(file="~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_gene_exons.gff3",
                                format="gff3",
                                dataSource="ftp://ftp.jgi-psf.org/pub/compgen/phytozome/v9.0/Bdistachyon/",
                                species="brachypodium")
saveDb(txdb, file="~/NGS/BD/BDproj/cache/Bd192.sqlite")
txdb <- loadDb("~/NGS/BD/BDproj/cache/Bd192.sqlite") 
columns(txdb)
keytypes(txdb)

#select(txdb, keys = keys, columns = cols, keytype = "GENEID")
TXCHROM"
select(txdb, keys = keys, columns = cols, keytype = "GENEID")

# then for each gene, reduce all the exons to a set of non overlapping exons,
# calculate their lengths (widths) and sum then

eByg <- exonsBy(txdb, by="gene")
length(eByg)
names(eByg)[10:13]
eByg[['Bradi0009s00250']]

test <- read.table("~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_gene_exons.gff3", header=FALSE)


subset(test, V1=="Bd4" & V4 >= 40965956 & V5 <= 40967004)

###################################################
### collect read count
###################################################

countDF <- data.frame(row.names=names(eByg))
for(i in samplespath) {
aligns <- readGAlignmentsFromBam("~/NGS/BD/4-29-14/Sample_3-2-12-1/3-2-12-1_CAGATC_L006_R1_001.bam") # Substitute next two lines with this one.
counts <- countOverlaps(eByg, aligns, ignore.strand=TRUE)
countDF <- cbind(countDF, counts)
}
colnames(countDF) <- samples
countDF[1:4,]
write.table(countDF, "./results/countDF", quote=FALSE, sep="\t", col.names = NA)
countDF <- read.table("./results/countDF")

###################################################
### Compute RPKM
###################################################
returnRPKM <- function(counts, gffsub) {
geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
return(rpkm)
}
countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
countDFrpkm[1:4,]
