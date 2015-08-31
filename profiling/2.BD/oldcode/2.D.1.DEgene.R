## Jinliang Yang
## 8.11.2014
## Using DESeq to analyze differentially expressed genes

ob <- load("~/Documents/BDproj/cache/count.RData")
#Raw count data are expected here!
library(DESeq)
#countDF <- read.table("./results/countDF")
targets <- read.csv("~/Documents/BDproj/data/target.csv")
conds <- targets$Factor
## Creates object of class CountDataSet derived from eSet class > counts(cds)[1:4, ]
# CountDataSet has similar accessor methods as eSet class.
cds <- newCountDataSet(countDF, conds) 

cds <- estimateSizeFactors(cds) # Estimates library size
cds <- estimateDispersions(cds) # Estimates the variance

#################
getDEgene <- function(treat1="CBF", treat2="WT", 
                      goplot=TRUE, plotfile="~/Documents/BDproj/reports/plotMA_CBF_WT.png"){
  # Calls DEGs with nbinomTest
  res <- nbinomTest(cds, treat1, treat2)
  res <- na.omit(res)
  if(goplot){
    par(mfrow=c(1,1))
    png(filename=plotfile, 
        width=480, height=480)
    plotMA(res, main=paste(treat1, "vs.", treat2))
    dev.off()
  }
  res2fold <- res[res$log2FoldChange >= 1 | res$log2FoldChange <= -1,]
  res2foldpadj <- res2fold[res2fold$padj <= 0.05, ]
  message(sprintf("# of DE genes [ %s ]", nrow(res2foldpadj)))
  #[1] 69  8
  return(res2foldpadj)
}

res1 <- getDEgene(treat1="CBF", treat2="WT", goplot=TRUE,
                  plotfile="~/Documents/BDproj/reports/plotMA_CBF_WT.png")

res2 <- getDEgene(treat1="CBF4C", treat2="WT4C", goplot=TRUE,
                  plotfile="~/Documents/BDproj/reports/plotMA_CBF4C_WT4C.png")

res3 <- getDEgene(treat1="CBF", treat2="CBF4C", goplot=TRUE,
                  plotfile="~/Documents/BDproj/reports/plotMA_CBF_CBF4C.png")

res4 <- getDEgene(treat1="WT", treat2="WT4C", goplot=TRUE,
                  plotfile="~/Documents/BDproj/reports/plotMA_WT_WT4C.png")

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("res1", "res2", "res3", "res4"), 
            file="~/Documents/BDproj/cache/DEgenes.RData",
            description=c("CBF vs. WT", "CBF4C vs. WT4C", "CBF vs. CBF4C", "WT vs. WT4C"))

write.table(res1, "~/Documents/BDproj/reports/Table1.DEgenes_CBF_WT.csv", sep=",", 
            row.names=FALSE, quote=FALSE)
write.table(res2, "~/Documents/BDproj/reports/Table2.DEgenes_CBF4C_WT4C.csv", sep=",", 
            row.names=FALSE, quote=FALSE)
write.table(res3, "~/Documents/BDproj/reports/Table3.DEgenes_CBF_CBF4C.csv", sep=",", 
            row.names=FALSE, quote=FALSE)
write.table(res4, "~/Documents/BDproj/reports/Table4.DEgenes_WT_WT4C.csv", sep=",", 
            row.names=FALSE, quote=FALSE)



