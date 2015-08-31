### Jinliang Yang
### count reads and conduct a chi-square test




### manually input the refgene file
bed3 <- read.csv("~/Documents/BDproj/data/refgenes.csv")
write.table(bed3, "~/NGS/BD/JJ-RNA-seq/7.8.2014_more_genes/refgene.bed3",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


bedtools intersect -abam bri.sorted.bam -b refgene.bed3 -bed > bri_rc.txt
bedtools intersect -abam cb.sorted.bam -b refgene.bed3 -bed > cb_rc.txt



read1 <- read.table("~/NGS/BD/JJ-RNA-seq/7.8.2014_more_genes/bri_rc.txt")
r1 <- subset(read1, V1=="BRI1_Triticum")
nrow(r1) #72, 527
#118
read2 <- read.table("~/NGS/BD/JJ-RNA-seq/7.8.2014_more_genes/cb_rc.txt")
r2 <- subset(read2, V1=="BRI1_Triticum")
nrow(r2) #130, 284
#197

totbri = 11926415;
totcb = 7357495;

M <- as.table(rbind(c(726, 480), c(16294052, 10295858)))
dimnames(M) <- list(val = c("gene","total"),
                    type = c("Bri","CB"))

M
chisq.test(M)
