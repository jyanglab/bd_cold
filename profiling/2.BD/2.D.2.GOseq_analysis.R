### Jinliang Yang
### Go term enrichment test!


### Category mapping format
getGO <- function(de.genes=res2foldpadj$id){
  library(goseq)
  ### length data
  txdb <- loadDb("~/Documents/BDproj/cache/Bd192.sqlite") 
  txsByGene <- transcriptsBy(txdb,"gene")
  lengthData <- median(width(txsByGene))
  assayed.genes <- names(eByg)
  
  gene.vector <- as.integer(assayed.genes %in% de.genes)
  names(gene.vector) <- assayed.genes
  
  ## GO file
  gos <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/brachy_V1.2_GOs_fromSIMAP.txt", header=FALSE)
  names(gos) <- c("txid", "GOterm", "function", "N")
  gosin <- nrow(gos)
  ## txid to geneid
  info <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_annotation_info.txt", header=F)
  names(info)[2:3] <- c("geneid", "txid")
  gos <- merge(info[,2:3], gos, by="txid" )
  message(sprintf("Input go terms/found geneid [ %s | %s ]", gosin, nrow(gos)))
  
  #Fitting the Probability Weighting Function (PWF)
  pwf <- nullp(gene.vector, bias.data=lengthData, plot.fit=TRUE)
  go <- goseq(pwf, gene2cat=gos[, 2:3], method="Wallenius")
  
  return(go)
}
#######

#supportedGenomes()

ob <- load("~/Documents/BDproj/cache/DEgenes.RData")
gos1 <- getGO(de.genes=res1$id) #44
gos2 <- getGO(de.genes=res2$id) #61 4
gos3 <- getGO(de.genes=res3$id) #74 5
gos4 <- getGO(de.genes=res4$id) #90 6

gos <- getGO(de.genes= row.names(set1))
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


merge(enriched, mygo[, 2:3], by.x="category", by.y="GOterm" )


## GO file
mygo <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/brachy_V1.2_GOs_fromSIMAP.txt", header=FALSE)
names(mygo) <- c("txid", "GOterm", "function", "N")
mygo <- mygo[!duplicated(mygo$GOterm),]



gosin <- nrow(gos)
## txid to geneid
info <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_annotation_info.txt", header=F)
names(info)[2:3] <- c("geneid", "txid")
gos <- merge(info[,2:3], gos, by="txid" )
write.table(gos, "~/Documents/BDproj/reports/brachy_V1.2_GOs_fromSIMAP+geneid.csv",
            sep=",", row.names=FALSE, quote=FALSE)

