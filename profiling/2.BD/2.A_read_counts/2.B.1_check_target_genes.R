### Jinliang Yang
### check target genes
### purpose: checking CBF3 and other 7 CBF homologs


ob <- load("~/Documents/BDproj/cache/count.RData")


targets <- read.csv("~/Documents/BDproj/data/target.csv")
conds <- targets$Factor

### CBF genes
cbf <- read.csv("~/Documents/BDproj/data/Bd_CBF_related.csv")
cbf <- cbf[1:8,]

colnames(rpkmrep1) <- conds
colnames(rpkmrep2) <- conds
colnames(rpkm) <- conds

cbf_exp1 <- merge(cbf, rpkmrep1, by.x="geneid", by.y="row.names", all.x=TRUE)
cbf_exp2 <- merge(cbf, rpkmrep2, by.x="geneid", by.y="row.names", all.x=TRUE)
cbf_exp <- merge(cbf, rpkm, by.x="geneid", by.y="row.names", all.x=TRUE)

### avergae three replications
exp <- cbf_exp[, 1:2]
exp$CBF <- (cbf_exp$CBF + cbf_exp$CBF.1 + cbf_exp$CBF.2)/3
exp$CBF4C <- (cbf_exp$CBF4C + cbf_exp$CBF4C.1 + cbf_exp$CBF4C.2)/3
exp$WT <- (cbf_exp$WT + cbf_exp$WT.1 + cbf_exp$WT.2)/3
exp$WT4C <- (cbf_exp$WT4C + cbf_exp$WT4C.1 + cbf_exp$WT4C.2)/3


library(lattice)
library(gplots)
targets <- read.csv("~/Documents/BDproj/data/target.csv")
colnames(y0) <- targets$Factor
##########
tcbf <- t(as.matrix(exp[-8, 3:6]))
colnames(tcbf) <- cbf_exp$geneid[-8]

### original
levelplot(tcbf, height=0.5, xlab="", ylab="", colorkey=list(space="top"),
          col.regions=colorpanel(40, "blue", "yellow", "red")) 


### rescaled
source("~/Documents/Rcodes/rescale.R")
for(i in 1:ncol(tcbf)){
  tcbf[, i] <- rescale(tcbf[, i], c(0,1))
}

levelplot(tcbf, height=0.5, xlab="", ylab="", colorkey=list(space="top"),
          col.regions=colorpanel(40, "blue", "yellow", "red")) 



