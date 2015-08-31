# Jinliang Yang
# Oct 5th, 2014
# MapMan annotation

ob <- load("~/Documents/BDproj/cache/deg_res.RData")

### new version with annotation
source("~/Documents/BDproj/lib/getMapMan.R")
######
setwd("~/Documents/BDproj/reports/")
mm1 <- getMapMan(DEG=res1)
#>>> DEG total [ 460 ], mapped [ 460 ]
write.table(mm1, "res1_cbf_wt_mapman_460.txt", sep="\t", row.names=FALSE, quote=FALSE)

mm2 <- getMapMan(DEG=res2)
#>>> DEG total [ 2839 ], mapped [ 2837 ]
write.table(mm2, "res2_cbf_cbf4c_mapman_2839.txt", sep="\t", row.names=FALSE, quote=FALSE)

mm3 <- getMapMan(DEG=res3)
#>>> DEG total [ 3213 ], mapped [ 3213 ]
write.table(mm3, "res3_wt_wt4c_mapman_3213.txt", sep="\t", row.names=FALSE, quote=FALSE)

mm4 <- getMapMan(DEG=res4)
#>>> DEG total [ 1871 ], mapped [ 1870 ]
write.table(mm4, "res4_cbf4c_wt4c_mapman_1871.txt", sep="\t", row.names=FALSE, quote=FALSE)

#source("~/Documents/Rcodes/save.append.R")
#save.append(file="~/Documents/BDproj/cache/deg_res.RData", mm1)



