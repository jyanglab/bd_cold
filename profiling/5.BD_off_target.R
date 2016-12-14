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


