# Jinliang Yang
# updated: 5.9.2014

### http://research-pub.gene.com/gmap/src/README
### Setting up to build a GMAP/GSNAP database (one chromosome per FASTA entry)
gmap_build -D ~/db -d CB_fro_gsnap ~/DBcenter/CB_FROs/fro_arabidopsis.fasta

### gsnap allow 2 mismatches every 36 bp and 3bp tails, parameters learned form eddy!
gsnap -D ~/db -d CB_fro_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output Bri_pro  Bri-57-1_GTCCGC.trimmed.fq
gsnap -D ~/db -d CB_fro_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output CB_pro  CB_CCGTCC.trimmed.fq

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS Bri_pro.unpaired_uniq > Bri_pro_uniq.bam
#CREATING A BAM INDEX FILE
samtools sort Bri_pro_uniq.bam Bri_pro_uniq.sorted
samtools index Bri_pro_uniq.sorted.bam Bri_pro_uniq.sorted.bam.bai

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS CB_pro.unpaired_uniq > CB_pro_uniq.bam
samtools sort CB_pro_uniq.bam CB_pro_uniq_sorted 
#CREATING A BAM INDEX FILE
samtools index CB_pro_uniq_sorted.bam CB_pro_uniq_sorted.bam.bai



############################################################################################################
## RUN de novo assembly
Trinity.pl --seqType fq --JM 10G --single  --CPU 6

"~/bin/trinity/trinityrnaseq_r20131110/util/RSEM_util/"
gmap -D ~/db -d CB_Bri1_only_gsnap  --split-output test  Trinity.fasta
