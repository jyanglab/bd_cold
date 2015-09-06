## Jinliang Yang
## 8.11.2014
##

ob <- load("~/Documents/Github/BDproj/cache/count.RData")


countDF1["Bradi4g35650", 1:6] <- cbf[,4:9]
countDF1["Bradi4g35650", 7:12] <- wt[,4:9]


###################################################
txdb <- loadDb("~/Documents/BDproj/cache/Bd192.sqlite") 
columns(txdb)
keytypes(txdb)

eByg <- exonsBy(txdb, by="gene")
length(eByg)
#[1] 26552

###################################################
### Compute RPKM
###################################################
returnRPKM <- function(counts, gffsub=eByg) {
  geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
  millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
  rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
  rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
  return(rpkm)
}

rpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
rpkm1 <- apply(countDF1, 2, function(x) returnRPKM(counts=x, gffsub=eByg))

cbf_1 <- rpkm1["Bradi4g35650", ]

plot(cbf)
tbcbf <- data.frame(cbf, cov, )

### manually changed the format
cbf <- read.csv("~/Documents/Github/BDproj/reports/cbf_rc.csv")
cbf <- cbf[1:4,]

t.test(cbf[1,2:4], cbf[2, 2:4])
t.test(cbf[1,2:4], cbf[2, 2:4])

t.test(cbf[3,2:4], cbf[4, 2:4])
t.test(cbf[1,2:4], cbf[2, 2:4])



#QC check of the sample reproducibility by computing a correlating matrix and plotting it as a tree.
#Note: the plotMDS function from edgeR is a more robust method for this task.
library(ape)
d <- cor(rpkm, method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("countDF1", "countDF2", "countDF", "rpkmrep1", "rpkmrep2", "rpkm"), 
            file="~/Documents/BDproj/cache/count.RData",
            description="")
