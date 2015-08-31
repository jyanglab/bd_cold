# Jinliang Yang
# 4/7/2014
# tested the EBARDenovo and found it could not handle large data


install.packages("ProjectTemplate")
library('ProjectTemplate')
setwd("~/NGS/BD")
create.project('BDproj')

### config

### http://research-pub.gene.com/gmap/src/README
# cd /mnt/02/yangjl/DBcenter/BD_v1.0
### Setting up to build a GMAP/GSNAP database (one chromosome per FASTA entry)
gmap_build -D ~/db/brachy1.0 -d brachy1.0_gsnap brachy1.0_wholegenome_unmasked.mfa
  
### gsnap allow 2 mismatches every 36 bp and 3bp tails
gsnap -d ~/db/AGPv2/RefGenv2_gsnap/ -B 2 -m 10 -i 2 -N 1 -w 10000 -t 8 -n 3 --quality-protocol=illumina ???nofails <input> > <output>
  
  


