### Jinliang Yang
### 9/4/2015

#http://floresta.eead.csic.es/footprintdb/index.php?tf=25f213b64bbd2f3012ccb61060dd270c
AykrCCGACmT #cbf1,2
TACTrCCGACAtGA #CBF3
### Motif Conversion Utilities
# largedata
system("iupac2meme -named GCCGAC cbf > cbf_core.meme")
system("iupac2meme -named TACTrCCGACAtGA cbf3 > cbf3.meme")


###
system("fimo --oc fimo_core_files --parse-genomic-coord cbf_core.meme Bdistachyon_192.fa")


#########################################
motif12 <- read.table("largedata/fimo_output_files/fimo.txt", header=FALSE)
motif3 <- read.table("largedata/cbf3_files/fimo.txt", header=FALSE)
motif <- rbind(motif12, motif3)
motif <- subset(motif, V2 %in% c("Bd1", "Bd2", "Bd3", "Bd4", "Bd5"))
#38780

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
    for(i in 1:6){
        tem <- subset(motif, V2 == "Bd4" & V3 > cbfgene$start[i]- upstream & V4 < cbfgene$start[i])
        print(sprintf("###>>> gene [ %s ] contains [ %s ] motif", cbfgene$attribute[i], nrow(tem)))
    }
    for(i in 7:8){
        tem <- subset(motif, V2 == "Bd4" & V3 > cbfgene$end[i] & V4 < cbfgene$end[i] + upstream)
        print(sprintf("###>>> gene [ %s ] contains [ %s ] motif", cbfgene$attribute[i], nrow(tem)))
    }
}

#### CBF3 motif
find_gene_motif(motif3, cbfgene, upstream=7000)

find_gene_motif(motif12, cbfgene, upstream=7000)


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
gene_reg <- find_reg(BP=5000, gene, motif)
write.table(gene_reg, "data/gene_motif.csv", sep=",", row.names=FALSE, quote=FALSE)

sum(gene_reg$CRT >0)
#13385


plot(cbf8$V2, cbf8$V3)

