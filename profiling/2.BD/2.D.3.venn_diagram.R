## Jinliang Yang
## 10.5.2014
## venn diagram

ob <- load("~/Documents/BDproj/cache/deg_res.RData")

source("~/Documents/Rcodes/overLapper.R")

#########
myres2 <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 1)
myres3 <- subset(res3, padj < 0.05 & abs(log2FoldChange) > 1)
setlist <- list(CBF_CBF4C=row.names(myres2), WT_WT4C=row.names(myres3))
OLlist <- overLapper(setlist=setlist, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)

par(mfrow=c(1,2))
vennPlot(counts=counts, mymain="")

#########
myres4 <- subset(res4, padj < 0.05 & abs(log2FoldChange) > 1)
ls2 <- as.data.frame(OLlist[[2]])
cbf <- subset(ls2, CBF_CBF4C !=1 & WT_WT4C ==1)
setlist <- list(Overlapped=row.names(cbf), CBF4C_WT4C=row.names(myres4))
OLlist <- overLapper(setlist=setlist, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts, mymain="")





###############################################
source("~/Documents/BDproj/lib/getMapMan.R")

### res3 + res4
mylist <- OLlist[[2]]
nmall <- row.names(mylist)

subres3 <- subset(res3, row.names(res3) %in% nmall)
subres4 <- subset(res4, row.names(res4) %in% nmall)

mmcbf_all <- getMapMan(DEG=rbind(subres3, subres4))
#>>> DEG total [ 4055 ], mapped [ 2569 ]
write.table(mmcbf_all, "~/Documents/BDproj/reports/res5_cbf_all_2569.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)

### res3 x res4
mylist <- OLlist[[2]]
mylist <- as.data.frame(mylist)
nmol <- row.names(subset(mylist, res23==1 & res3 ==1))

subres3 <- subset(res3, row.names(res3) %in% nmol)
subres4 <- subset(res4, row.names(res4) %in% nmol)

mmcbf_ol <- getMapMan(DEG=rbind(subres3, subres4))
#>>> DEG total [ 1464 ], mapped [ 732 ]
write.table(mmcbf_ol, "~/Documents/BDproj/reports/res6_cbf_ol_732.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)
