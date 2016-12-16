### Jinliang Yang
### Dec 13th, 2016
### Off-target analysis


library("Biostrings")

fa <- readDNAStringSet(filepath = "data/RNAi_CBF3_Target.fa", format="fasta")
class(fa)
length(fa[[1]])
width(fa)


out <- data.frame()
for(i in 1:(width(fa)-18)){
   tmp <-  fa[[1]][i:(i+18)]
   
   d <- data.frame(id=i, seq=as.character(tmp))
   
   out <- rbind(out, d)
}

write.table(out, "data/RNAi_CBF3_19bp.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#########################

library("Biostrings")
ds_off_target <- function(target=fa, genome=cds){
    out <- data.frame()
    for(i in 1:(width(target)-18)){
        tmp <-  target[[1]][i:(i+18)]
        res <- vmatchPattern(pattern=tmp, subject=genome, max.mismatch=1)
        out1 <- data.frame(gene=names(res), mymatch=unlist(lapply(res, length)))
        
        out2 <- subset(out1, mymatch > 0)
        message(sprintf("[ds_off_target] running for read [ %s ] and found [ %s ] matches ...", i, nrow(out2)))
        if(nrow(out2) > 0){
            out <- rbind(out, out2)
        }
    }
    return(out)
}
#####
cds <- readDNAStringSet(filepath = "largedata/Bradi_1.0.cds.fa.gz", format="fasta")
length(cds)
# [1] 32255
fa <- readDNAStringSet(filepath = "data/RNAi_CBF3_Target.fa", format="fasta")

output <- ds_off_target(target=fa, genome=cds)
myres <- subset(output, gene != "Bradi4g35650.1")

out <- data.frame(table(myres$gene))
subset(out, Freq > 0)
