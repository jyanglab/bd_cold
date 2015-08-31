#####################################
getMapMan <- function(DEG=res1, idcol = TRUE){
  mmid <- read.csv("~/Documents/BDproj/data/Bdistachyon_192_MapMan.csv")
  id <- mmid[, 2:3]
  id$gene <- gsub("\\..*", "", id$IDENTIFIER)
  message(">>> Annotation version: [ Bdistachyon_192 ]")
  
  #########
  myres <- subset(DEG, padj < 0.05 & abs(log2FoldChange) > 1)
  message(">>> Filtering padj < 0.05 and log2FoldChange > 1")
  
  if(idcol){
    wt4 <- as.data.frame(myres)
    wt4$id <- tolower(wt4$id)
  }else{
    wt4 <- as.data.frame(myres)
    wt4$id <- row.names(wt4)
    wt4$id <- tolower(wt4$id)
    
  }
  
  wt42 <- merge(id[, c("gene", "IDENTIFIER")], wt4, by.x="gene", by.y="id", all.y=TRUE)
  wt42 <- wt42[, c(2, 4, 1,3,5:8)]
  
  ####### find annotation
  ## txid to geneid
  info <- read.delim("~/DBcenter/BD_v1.0/annot_v1.2/Bdistachyon_192_annotation_info.txt", header=F)
  
  #1: Phytozome internal transcript ID (potentially useful to connect to biomart datasets)
  #2: Phytozome gene locus name
  #3: Phytozome transcript name
  #4: Phytozome protein name (often same as transcript name, but this can vary)
  #5: PFAM
  #6: Panther
  #7: KOG
  #8: KEGG ec
  #9: KEGG Orthology
  #10: Gene Ontology terms (NOTE: these are automated results from interpro2go in most genomes, *not* empirically derived)
  #11: best arabidopsis TAIR10 hit name
  #12: best arabidopsis TAIR10 hit symbol
  #13: best arabidopsis TAIR10 hit defline
  #14: best rice hit name
  #15: best rice hit symbol
  #16: best rice hit defline
  info <- info[, c(2, 10, 13, 16)]
  names(info) <- c("gene", "GO", "arabidopsis", "rice")
  info$gene <- tolower(info$gene)
  wt42in <- merge(wt42, info, by="gene", all.x=TRUE)
  
  #### formatting!!!
  wt42in[is.infinite(wt42in$log2FoldChange), ] <- "X"
  wt42in[is.na(wt42in$log2FoldChange), ] <- "X"
  wt42in <- wt42in[!duplicated(wt42in$gene),]
  wt42in[is.na(wt42in$IDENTIFIER),]$IDENTIFIER <- wt42in[is.na(wt42in$IDENTIFIER),]$gene
  wt42in <- wt42in[!duplicated(wt42in$IDENTIFIER),]
  wt42in <- wt42in[, -1]
  
  message(sprintf("#>>> DEG total [ %s ], mapped [ %s ]", nrow(myres), nrow(wt42in)))
  return(wt42in)
}
