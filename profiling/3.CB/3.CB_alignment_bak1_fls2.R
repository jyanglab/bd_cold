# Jinliang Yang
# updated: 5.9.2014

### http://research-pub.gene.com/gmap/src/README
### Setting up to build a GMAP/GSNAP database (one chromosome per FASTA entry)
gmap_build -D ~/db -d CB_bak1_fls2_gsnap ~/DBcenter/CB_bak1_fls2/BAK1_FLS2.fasta

### gsnap allow 2 mismatches every 36 bp and 3bp tails, parameters learned form eddy!
gsnap -D ~/db -d CB_bak1_fls2_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output Bri_bak1_fls2  Bri-57-1_GTCCGC.trimmed.fq
gsnap -D ~/db -d CB_bak1_fls2_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output CB_bak1_fls2  CB_CCGTCC.trimmed.fq

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS Bri_bak1_fls2.unpaired_uniq > Bri_bak1_fls2_uniq.bam
#CREATING A BAM INDEX FILE
samtools sort Bri_bak1_fls2_uniq.bam Bri_bak1_fls2_uniq.sorted
samtools index Bri_bak1_fls2_uniq.sorted.bam Bri_bak1_fls2_uniq.sorted.bam.bai

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS CB_bak1_fls2.unpaired_uniq > CB_bak1_fls2_uniq.bam
samtools sort CB_bak1_fls2_uniq.bam CB_bak1_fls2_uniq_sorted 
#CREATING A BAM INDEX FILE
samtools index CB_bak1_fls2_uniq_sorted.bam CB_bak1_fls2_uniq_sorted.bam.bai


############################################################################################################
## RUN de novo assembly
Trinity.pl --seqType fq --JM 10G --single  --CPU 6

"~/bin/trinity/trinityrnaseq_r20131110/util/RSEM_util/"
gmap -D ~/db -d CB_Bri1_only_gsnap  --split-output test  Trinity.fasta
