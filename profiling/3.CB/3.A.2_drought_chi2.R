### Jinliang Yang
### count reads and conduct a chi-square test

### manually input the refgene file
bed3 <- read.csv("~/Documents/BDproj/data/drought_refgen_bed3.csv")
write.table(bed3, "/mnt/02/yangjl/NGS/BD/Bri_drought/drought_refgene.bed3",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


#####
bedtools multicov -bams \
~/NGS/BD/7-31-14/rep1/Sample_Bri-57-1/Bri-57-1_GTCCGC.sorted.bam \
~/NGS/BD/7-31-14/rep2/Sample_Bri-57-1/Bri-57-1_GTCCGC.sorted.bam \
~/NGS/BD/7-31-14/rep1/Sample_CB/CB_CCGTCC.sorted.bam \
~/NGS/BD/7-31-14/rep2/Sample_CB/CB_CCGTCC.sorted.bam \
-bed drought_refgene.bed3 > drought_rc.txt




###########
rc <- read.table("~/NGS/BD/Bri_drought/drought_rc.txt", header=FALSE)
names(rc) <- c("gene", "start", "end", "rep1_bri", "rep2_bri", "rep1_wt", "rep2_wt")


tot <- data.frame(rep1_bri=11957506, rep2_bri=16294052, 
                  rep1_wt=7373802, rep2_wt=10295858)

getChisq <- function(rc=rc, tot=tot){
  rc$pval <- NA
  for(i in 1:nrow(rc)){
    M <- as.table(rbind(c(rc$rep1_bri[i]+rc$rep2_bri[i], rc$rep1_wt[i]+rc$rep2_wt[i]), 
                        c(tot$rep1_bri+tot$rep2_bri, tot$rep1_wt+tot$rep2_wt)))
    dimnames(M) <- list(val = c("gene","total"),
                        type = c("Bri","CB"))
    rc$pval[i] <- chisq.test(M)$p.value
    
  }
  return(rc)
}

res <- getChisq(rc=rc, tot=tot)
write.table(res, "~/Documents/BDproj/reports/Bri_droght_DEgene.csv", sep=",",
            row.names=FALSE, quote=FALSE)  
  
  