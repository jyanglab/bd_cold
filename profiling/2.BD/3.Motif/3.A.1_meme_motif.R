### Jinliang Yang
### 9/4/2015

#http://floresta.eead.csic.es/footprintdb/index.php?tf=25f213b64bbd2f3012ccb61060dd270c
AykrCCGACmT

### Motif Conversion Utilities
# largedata
system("iupac2meme -named AykrCCGACmT cbf > cbf.meme")

###
system("fimo --oc fimo_output_files --parse-genomic-coord cbf.meme Bdistachyon_192.fa")


#########################################

motif <- read.table("largedata/fimo_output_files/fimo.txt", header=FALSE)
motif <- subset(motif, V2 %in% c("Bd1", "Bd2", "Bd3", "Bd4", "Bd5"))
#39231

gff <- read.table("largedata/annot_v1.2/Bdistachyon_192_gene_exons.gff3", header=FALSE)
names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")
gene <- subset(gff, feature == "gene")
gene$attribute <- gsub(".*Name=", "", gene$attribute)

cbf <- c("Bradi4g35630", "Bradi4g35650","Bradi4g35570", "Bradi4g35580",
         "Bradi4g35590", "Bradi4g35600", "Bradi4g35610", "Bradi4g35620")

cbfgene <- subset(gene, attribute %in% cbf)

cbf8 <- read.csv("data/CBF_8genes_bed.csv", header=FALSE)

r <- range(c(cbf8$V2, cbf8$V3))
subset(motif, V2 == "Bd4" & V3 > r[1]-5000 & V4 < r[2]+10000)



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

