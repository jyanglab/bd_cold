## Jinliang Yang
## 8.11.2014
## Using DESeq to analyze differentially expressed genes

ob <- load("~/Documents/Github/BDproj/cache/count.RData")
#Raw count data are expected here!
library(DESeq2)


#countDF <- read.table("./results/countDF")
targets <- read.csv("~/Documents/Github/BDproj/data/target.csv")
tgs1 <- tgs2 <- targets
tgs1$Rep <- "Rep1"
tgs2$Rep <- "Rep2"
tgs <- rbind(tgs1, tgs2)

rc2tech <- cbind(countDF1,countDF2)

##########################
DESeq2_pipe <- function(rcdata=rc2tech[, c(7:12, 21:26)], design_table=tgs[c(7:12, 21:26),], 
                        design_model= ~ Rep + Treatment){
  ########## WT and WT4C
  dds <- DESeqDataSetFromMatrix(as.matrix(rcdata), colData=design_table,
                                design = design_model)
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  res <- as.data.frame(res)
  
  sig <- subset(res, padj < 0.05 & abs(log2FoldChange)>1)
  #myres <- subset(res, padj < 0.05 & abs(log2FoldChange)>2)
  message(sprintf("#>>> DE genes: [ %s ], up: [%s], down: [%s]", nrow(sig), 
                  nrow(subset(sig, log2FoldChange >0)), nrow(subset(sig, log2FoldChange <0)) ))
  return(resOrdered)
}


######## CBF and WT under nonstress
cols <- rows <- c(1:3,7:9, 15:17, 21:23)
res1 <- DESeq2_pipe(rcdata=rc2tech[, cols], design_table=tgs[rows,], 
                    design_model= ~ Rep + Type)
#>>> DE genes: [ 460 ], up: [299], down: [161]
sig <- subset(res1, padj < 0.05 & abs(log2FoldChange)>1)


########## CBF and CBF4C
cols <- rows <- c(1:6, 15:20)
res2 <- DESeq2_pipe(rcdata=rc2tech[, cols], design_table=tgs[rows,], 
                    design_model= ~ Rep + Treatment)
#>>> DE genes: [ 2839 ], up: [1335], down: [1504]
sig2 <- subset(res2, padj < 0.05 & abs(log2FoldChange)>1)


######### WT and WT4C
cols <- rows <- c(7:12, 21:26)
res3 <- DESeq2_pipe(rcdata=rc2tech[, cols], design_table=tgs[rows,], 
                    design_model= ~ Rep + Treatment)
#>>> DE genes: [ 3213 ], up: [1409], down: [1804]
sig3 <- subset(res3, padj < 0.05 & abs(log2FoldChange)>1)

######### CBF4C and WT4C
cols <- rows <- c(3:6, 10:12, 18:20, 24:26)
res4 <- DESeq2_pipe(rcdata=rc2tech[, cols], design_table=tgs[rows,], 
                    design_model= ~ Rep + Type)
#>>> DE genes: [ 1871 ], up: [1088], down: [783]
sig3 <- subset(res3, padj < 0.05 & abs(log2FoldChange)>1)

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("res1", "res2", "res3", "res4"), file="~/Documents/BDproj/cache/deg_res.RData",
            description=c("CBF vs. WT, DEG=460", "CBF vs. CBF4C, DEG=2,839",
                   "WT vs. WT4C, DEG=3,213", "CBF4C vs. WT4C, DEG=1,871"))




