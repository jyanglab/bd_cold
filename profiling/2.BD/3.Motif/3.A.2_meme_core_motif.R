### Jinliang Yang
### 9/4/2015

#http://floresta.eead.csic.es/footprintdb/index.php?tf=25f213b64bbd2f3012ccb61060dd270c
AykrCCGACmT #cbf1,2
TACTrCCGACAtGA #CBF3
### Motif Conversion Utilities
# largedata
system("iupac2meme -named rCCGAC cbf > cbf_core.meme")
system("iupac2meme -named TACTrCCGACAtGA cbf3 > cbf3.meme")


###
system("fimo --oc fimo_core_files --max-stored-scores 10000 --no-qvalue --parse-genomic-coord cbf_core.meme Bdistachyon_192.fa ")
#?mast

system("mast cbf_core.meme Bdistachyon_192.fa -oc mast_core_files -nohtml -comp")


#########################################
motif12 <- read.table("largedata/fimo_output_files/fimo.txt", header=FALSE)
motif3 <- read.table("largedata/cbf3_files/fimo.txt", header=FALSE)
motif <- rbind(motif12, motif3)
motif <- subset(motif, V2 %in% c("Bd1", "Bd2", "Bd3", "Bd4", "Bd5"))
#38780
motif$bin1 <- paste(motif$V2, round(motif$V3/10, 0), sep="_")
motif$bin2 <- paste(motif$V2, round(motif$V4/10, 0), sep="_")
motif <- motif[!duplicated(motif$bin1), ]
motif <- motif[!duplicated(motif$bin2), ]
write.table(motif[, 1:9], "data/core_motif.csv", sep=",", row.names=FALSE, quote=FALSE)


gff <- read.table("largedata/annot_v1.2/Bdistachyon_192_gene_exons.gff3", header=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")
gene <- subset(gff, feature == "gene")
gene$attribute <- gsub(".*Name=", "", gene$attribute)

cbf <- c("Bradi4g35630", "Bradi4g35650","Bradi4g35570", "Bradi4g35580",
         "Bradi4g35590", "Bradi4g35600", "Bradi4g35610", "Bradi4g35620")

cbfgene <- subset(gene, attribute %in% cbf)

r <- range(c(cbfgene$start, cbfgene$end))



find_gene_motif <- function(motif, cbfgene, upstream=7000){
    
    res <- data.frame()
    for(i in 1:6){
        tem <- subset(motif, V2 == "Bd4" & V3 > cbfgene$start[i]- upstream & V4 < cbfgene$start[i])
        print(sprintf("###>>> gene [ %s ] contains [ %s ] motif", cbfgene$attribute[i], nrow(tem)))
        if(nrow(tem)>0){
            tem$gene <- cbfgene$attribute[i]
            res <- rbind(res, tem)
        }
    }
    for(i in 7:8){
        tem <- subset(motif, V2 == "Bd4" & V3 > cbfgene$end[i] & V4 < cbfgene$end[i] + upstream)
        print(sprintf("###>>> gene [ %s ] contains [ %s ] motif", cbfgene$attribute[i], nrow(tem)))
        if(nrow(tem)>0){
            tem$gene <- cbfgene$attribute[i]
            res <- rbind(res, tem)
        }
    }
    return(res)
}

#### CBF3 motif
find_gene_motif(motif3, cbfgene, upstream=7000)

find_gene_motif(motif12, cbfgene, upstream=7000)

out <- find_gene_motif(motif, cbfgene, upstream=7000)
names(out) <- c("type", "chr", "start", "end", "strand", "value", "p", "cp", "seq", "gene")
write.table(out[, -8], "reports/motifs_in_cbfgene.csv", sep=",", row.names=FALSE, quote=FALSE)
###########
find_reg <- function(BP=5000, gene, motif){
    gene$seqname <- as.character(gene$seqname)
    gene$CRT <- 0
    for(i in 1:nrow(gene)){
        print(i)
        if(gene$strand[i] == "+"){
            tem <- subset(motif, V2 == gene$seqname[i] & V3 > gene$start[i]- BP & V3 <gene$start[i])
            gene$CRT[i] <- nrow(tem)
        }
        if(gene$strand[i] == "-"){
            tem <- subset(motif, V2 == gene$seqname[i] & V3 > gene$end[i] & V3 < gene$end[i] + BP)
            gene$CRT[i] <- nrow(tem)
        }
    }
    return(gene)
}
#########
motif3 <- read.table("largedata/cbf3_files/fimo.txt", header=FALSE)
gene_reg <- find_reg(BP=5000, gene, motif3)
write.table(gene_reg, "data/gene_with_motif.csv", sep=",", row.names=FALSE, quote=FALSE)

sum(gene_reg$CRT >0)
#13385

plot(cbf8$V2, cbf8$V3)

gm <- read.csv("data/gene_motif.csv")
names(gm)[10] <- "motif"
gm$motif <- as.numeric(as.character(gm$motif))
write.table(subset(gm, motif > 0), "data/Table_gene_with_motif.csv", sep=",", row.names=FALSE, quote=FALSE)

cbf <- read.delim("data/Table_S4_cbf_cbf4_2839.txt", header=TRUE)
wt <- read.delim("data/Table_S3_wt_wt4_3213.txt", header=TRUE)
cbf$ID <- gsub("\\..*", "", cbf$IDENTIFIER)
wt$ID <- gsub("\\..*", "", wt$IDENTIFIER)

gm$attribute <- tolower(gm$attribute)

sum(gm$attribute %in% cbf$ID)
sum(gm$attribute %in% wt$ID)
