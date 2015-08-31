## Jinliang Yang
## 10.15.2014
## get venn gene sets

ob <- load("~/Documents/BDproj/cache/deg_res.RData")
source("~/Documents/Rcodes/overLapper.R")
source("~/Documents/BDproj/lib/getMapMan.R")

######### venn plot 1
myres2 <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 1)
myres3 <- subset(res3, padj < 0.05 & abs(log2FoldChange) > 1)
setlist23 <- list(CBF_CBF4C=row.names(myres2), WT_WT4C=row.names(myres3))
OLlist23 <- overLapper(setlist=setlist23, sep="_", type="vennsets")
counts <- sapply(OLlist23$Venn_List, length)


### CBF-independent gene set: res2xres3
mylist <- OLlist23[[2]]
mylist <- as.data.frame(mylist)
cbf2x3 <- row.names(subset(mylist, CBF_CBF4C == 1 & WT_WT4C ==1))

subres2 <- subset(res2, row.names(res2) %in% cbf2x3)
subres3 <- subset(res3, row.names(res3) %in% cbf2x3)

mmcbf23 <- getMapMan(DEG=rbind(subres2, subres3), idcol = TRUE)
#>>> DEG total [ 2714 ], mapped [ 1357 ]
write.table(mmcbf23, "~/Documents/BDproj/reports/res5_cbf_independent_1357.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)

### CBF-compensation gene set: res2-res3
cbf2no3 <- row.names(subset(mylist, CBF_CBF4C == 1 & WT_WT4C !=1))
length(cbf2no3) #[1] 1482

subres2 <- subset(res2, row.names(res2) %in% cbf2no3)
subres3 <- subset(res3, row.names(res3) %in% cbf2no3)

mmcbf2no3 <- getMapMan(DEG=rbind(subres2, subres3), idcol = FALSE)
#>>> DEG total [ 1482 ], mapped [ 1481 ]
write.table(mmcbf2no3, "~/Documents/BDproj/reports/res5_cbf_compensation_1482.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)


### CBF-dependent gene set: res3-res2
cbf3no2 <- row.names(subset(mylist, CBF_CBF4C != 1 & WT_WT4C ==1))
length(cbf3no2) #[1] 1856

subres2 <- subset(res2, row.names(res2) %in% cbf3no2)
subres3 <- subset(res3, row.names(res3) %in% cbf3no2)

mmcbf3no2 <- getMapMan(DEG=rbind(subres3, subres2))
#>>> DEG total [ 1856 ], mapped [ 1856 ]
write.table(mmcbf3no2, "~/Documents/BDproj/reports/res5_cbf_dep_1856.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)


### res23 x res4
myres4 <- subset(res4, padj < 0.05 & abs(log2FoldChange) > 1)
cbfdep <- list(res23=cbf3no2, res4=row.names(myres4))
OLlist4 <- overLapper(setlist=cbfdep, sep="_", type="vennsets")
counts4 <- sapply(OLlist4$Venn_List, length)

mylist <- as.data.frame(OLlist4[[2]])

###### confident set
confid <- subset(mylist, res23==1 & res4==1)

subres3 <- subset(res3, row.names(res3) %in% row.names(confid))
subres4 <- subset(res4, row.names(res4) %in% row.names(confid))
deg34 <- rbind(subres3, subres4)
deg34 <- deg34[order(deg34$pvalue),]

mmcbf34 <- getMapMan(DEG=deg34)
#>>> DEG total [ 1464 ], mapped [ 733 ]
write.table(mmcbf_ol, "~/Documents/BDproj/reports/res6_cbf_confid_732.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)

###### less confident set
subres3 <- as.data.frame(subset(res3, row.names(res3) %in% row.names(mylist)))
subres3$id <- row.names(subres3)
subres4 <- as.data.frame(subset(res4, row.names(res4) %in% row.names(mylist)))
subres4$id <- row.names(subres4)
deg43 <- rbind(subres3, subres4)
#deg34 <- deg34[order(deg34$padj),]

mmcbf34 <- getMapMan(DEG=deg43, idcol = TRUE)
#>>> DEG total [ 4055 ], mapped [ 2995 ]
write.table(mmcbf34, "~/Documents/BDproj/reports/res6_cbf_noconfid_2995.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)









################################################


source("~/Documents/BDproj/lib/getGO.r")

gos <- getGO(de.genes= row.names(set2))
#By default, goseq tests all three major Gene Ontology branches; Cellular Components, Biological
#Processes and Molecular Functions. However, it is possible to limit testing to any combination
#of the major branches by using the test.cats argument to the goseq function. This is done by
#specifying a vector consisting of some combination of the strings“GO:CC”, “GO:BP”and“GO:MF”.
#For example, to test for only Molecular Function GO categories:

enriched <- gos[gos$over_represented_pvalue < 0.05, ]
nrow(enriched)
enriched2 <- gos[gos$under_represented_pvalue < 0.05, ]
nrow(enriched2)

library(GO.db)
for(gi in enriched$category){
  print(GOTERM[[gi]])
  cat("--------------------------------------\n")
}


#############
ob2 <- load("~/Documents/BDproj/cache/count.RData")
library(lattice)
library(gplots)






hmpplot <- function(res=res1, cols=c(1:3,7:9), 
                    outfile="~/Documents/BDproj/reports/Fig5.hmp_CBF_WT.pdf"){
  
  
}

res=res1; ht=20
cols=c(1:3,7:9)
outfile="~/Documents/BDproj/reports/Fig5.hmp_CBF_WT.pdf"
####

res=res2; ht=150
cols=c(4:6,10:12)
outfile="~/Documents/BDproj/reports/Fig6.hmp_CBFC4_WTC4.pdf"

res=res3; ht=200
cols=1:6
outfile="~/Documents/BDproj/reports/Fig7.hmp_CBF_CBFC4.pdf"

res=res4; ht=250
cols=7:12
outfile="~/Documents/BDproj/reports/Fig8.hmp_WT_WTC4.pdf"


###############################
y0 <- rpkm[rownames(rpkm) %in% res$id,]
targets <- read.csv("~/Documents/BDproj/data/target.csv")
colnames(y0) <- targets$Factor
##########
y <- y0[, cols]
y <- t(scale(t(as.matrix(y))))
y <- y[order(y[,1]),]
pdf(outfile, width=5, heigh=ht)
levelplot(t(y), height=0.5, xlab="", ylab="", colorkey=list(space="top"),
          col.regions=colorpanel(40, "darkblue", "yellow", "white")) 
dev.off()


