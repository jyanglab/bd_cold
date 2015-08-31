# Jinliang Yang
# updated: 5.9.2014

### http://research-pub.gene.com/gmap/src/README
### Setting up to build a GMAP/GSNAP database (one chromosome per FASTA entry)
# gmap_build -D ~/db -d CB_Bri1_gsnap ~/DBcenter/CB_Bri1/Bri1_gene.fasta
# gmap_build -D ~/db -d CB_Bri1_only_gsnap ~/DBcenter/CB_Bri1/Bri1_only.fasta

### gsnap allow 2 mismatches every 36 bp and 3bp tails, parameters learned form eddy!
gsnap -D ~/db -d CB_Bri1_only_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output Bri-57-1  Bri-57-1_GTCCGC.trimmed.fq
gsnap -D ~/db -d CB_Bri1_only_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output CB  CB_CCGTCC.trimmed.fq

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS Bri-57-1.unpaired_uniq > Bri-57-1_uniq.bam
#CREATING A BAM INDEX FILE
samtools sort Bri-57-1_uniq.bam Bri-57-1_uniq.sorted
samtools index Bri-57-1_uniq.sorted.bam Bri-57-1_uniq.sorted.bam.bai

#CONVERTING SAM DIRECTLY TO A SORTED BAM FILE
samtools view -bS CB.unpaired_uniq > CB_uniq.bam
samtools sort CB_uniq.bam CB_uniq_sorted 
#CREATING A BAM INDEX FILE
samtools index CB_uniq_sorted.bam CB_uniq_sorted.bam.bai



############################################################################################################
## RUN de novo assembly
Trinity.pl --seqType fq --JM 10G --single  --CPU 6

"~/bin/trinity/trinityrnaseq_r20131110/util/RSEM_util/"
gmap -D ~/db -d CB_Bri1_only_gsnap  --split-output test  Trinity.fasta
