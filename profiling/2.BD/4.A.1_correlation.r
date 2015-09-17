### Jinliang Yang
### Sep 13th, 2015


QPCR <- function(doplot=FALSE){
    qpcr <- read.csv("data/qPCR_results_validation.csv")
    qpcr <- subset(qpcr, sample != "cbf")
    #qpcr$SampleName <- toupper(qpcr$SampleName)
    
    control <- subset(qpcr, Assay == "GAPDH")
    assay <- subset(qpcr, Assay != "GAPDH")
    
    ex <- merge(assay, control[, 2:3], by="SampleName")
    ex$exp <- 2^( ex$CqMean.y - ex$CqMean.x)
    
    ex$Assay <- as.character(ex$Assay)
    genes <-unique(ex$Assay)
    
    ####
    ob <- load("~/Documents/Github/BDproj/cache/count.RData")
    
    #countDF <- read.table("./results/countDF")
    targets <- read.csv("~/Documents/Github/BDproj/data/target.csv")
    tem <- as.data.frame(rpkm[genes, ])
    names(tem) <- targets$SampleName
    
    out <- data.frame()
    for(i in 1:10){
        sub1 <- subset(ex, Assay == genes[i])
        sub1$method <- "qpcr"
        
        sub2 <- tem[genes[i], ]
        sub2 <- as.data.frame(t(sub2))
        sub2$sampleid <- row.names(sub2)
        sub2$sampleid <- gsub("CBF3_sample1", "CBF_A", sub2$sampleid)
        sub2$sampleid <- gsub("CBF3_sample2", "CBF_B", sub2$sampleid)
        sub2$sampleid <- gsub("CBF3_sample3", "CBF_C", sub2$sampleid)
        sub2$sampleid <- gsub("CBF3_4C_sample1", "CBF_A_4", sub2$sampleid)
        sub2$sampleid <- gsub("CBF3_4C_sample2", "CBF_B_4", sub2$sampleid)
        sub2$sampleid <- gsub("CBF3_4C_sample3", "CBF_C_4", sub2$sampleid)
        sub2$sampleid <- gsub("WT_sample1", "BD21_A", sub2$sampleid)
        sub2$sampleid <- gsub("WT_sample2", "BD21_B", sub2$sampleid)
        sub2$sampleid <- gsub("WT_sample3", "BD21_C", sub2$sampleid)
        sub2$sampleid <- gsub("WT_4C_sample1", "BD21_A_4", sub2$sampleid)
        sub2$sampleid <- gsub("WT_4C_sample2", "BD21_B_4", sub2$sampleid)
        sub2$sampleid <- gsub("WT_4C_sample3", "BD21_C_4", sub2$sampleid)
        sub2 <- subset(sub2, sampleid %in% sub1$SampleName)
        
        c1 <- cbind(subset(sub1, temp==4)$exp, subset(sub1, temp==23)$exp)
        
        sub2$temp <- 23
        sub2$temp <- gsub(".*A_|.*B_|.*C_", "",  sub2$sampleid)
        sub2$temp <- gsub(".*A|.*B|.*C", "23",  sub2$temp)
        names(sub2)[1] <- "exp"
        sub2$method <- "rnaseq"
        
        temp <- merge(sub2, sub1[, c("SampleName", "Assay", "Category", "exp")], by.x="sampleid", by.y="SampleName")
        c2 <- cbind(subset(sub2, temp==4)$exp, subset(sub2, temp==23)$exp)
        
        out <- rbind(out, temp)
        if(doplot){
            par(mfrow=c(2, 5))
            barplot(c1, beside=TRUE, main=genes[i])
            barplot(c2, beside=TRUE, main=genes[i])
        }
        
    }
    return(out)
    
}


######
res <- QPCR(doplot=FALSE)
write.table(res, "data/qpcr_rnaseq_cv.csv", sep=",", row.names=FALSE, quote=FALSE)

cor.test(subset(res, temp==4)$exp.x, subset(res, temp==4)$exp.y)
cor.test(subset(res, temp==23)$exp.x, subset(res, temp==23)$exp.y)

cor.test(res$exp.x, res$exp.y)

