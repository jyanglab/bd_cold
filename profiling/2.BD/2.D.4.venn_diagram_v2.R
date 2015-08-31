## Jinliang Yang
## 8.11.2014
##

ob <- load("~/Documents/BDproj/cache/DEgenes.RData")

source("~/Documents/Rcodes/overLapper.R")
setlist <- list(res2=res2$id, res3=res3$id, res4=res4$id)
OLlist <- overLapper(setlist=setlist, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts, mymain="DEG Comparison")
                                                                                                              


#### confident set of CBF dependent DE genes
sets <- as.data.frame(OLlist[[2]])
set1 <- subset(sets, res2==1 & res4==1 & res3==0)

#### extended set of CBF dependent DE genes
set2 <- subset(sets, (res2==1 | res4==1) & res3==0)


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


