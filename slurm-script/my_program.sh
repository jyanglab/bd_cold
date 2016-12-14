gmap_build -D largedata/ -d CB_fro_gsnap data/cenh3.fasta
gsnap -D largedata/ -d CB_fro_gsnap -m 10 -i 2 -N 1 -w 10000 -t 4 -n 3 -A sam --quality-protocol=sanger --nofails --split-output ~/dbcenter/CB/CB_CCGTCC_L006_R1_001 ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.fastq

samtools view -bS ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.unpaired_uniq > ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.bam
samtools sort ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.bam ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.sorted
samtools index ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.sorted.bam ~/dbcenter/CB/CB_CCGTCC_L006_R1_001.sorted.bam.bai

