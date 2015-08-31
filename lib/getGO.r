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