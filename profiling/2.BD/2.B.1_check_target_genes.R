### Jinliang Yang
### check target genes


ob <- load("~/Documents/BDproj/cache/count.RData")


targets <- read.csv("~/Documents/BDproj/data/target.csv")
conds <- targets$Factor

cbf <- read.csv("~/Documents/BDproj/data/Bd_CBF_related.csv")

names(countDF1) <- conds
names(countDF2) <- conds
colnames(rpkm) <- conds
names(rpkmrep1) <- conds
names(rpkmrep2) <- conds

idx <- grep("Bradi2g59500", row.names(rpkm))
idx


cbf_exp <- merge(cbf, rpkm, by.x="geneid", by.y="row.names", all.x=TRUE)

gff <- read.table("~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_gene_exons.gff3", header=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")
#2  NC_016132.1 (57180604..57183197, complement)
subset(gff, seqname=="Bd2" & start >= 57180604 & end <= 57183197 )


countDF1[idx:(idx+5),]
countDF2[idx,]
rpkmrep1[idx,]
rpkmrep2[idx,]
rpkm[idx,]




library(lattice)
library(gplots)
targets <- read.csv("~/Documents/BDproj/data/target.csv")
colnames(y0) <- targets$Factor
##########
tcbf <- t(as.matrix(cbf_exp[, 6:17]))
colnames(tcbf) <- cbf_exp$name
tcbf <- tcbf[,-4]

source("~/Documents/Rcodes/rescale.R")
for(i in 1:ncol(tcbf)){
  tcbf[, i] <- rescale(tcbf[, i], c(0,1))
}

levelplot(tcbf, height=0.5, xlab="", ylab="", colorkey=list(space="top"),
          col.regions=colorpanel(40, "yellow", "red")) 


