### Jinliang Yang
### 8/29/2014
### extract reads four CBF3 gene only



#####
bedtools intersect [OPTIONS] -a <BED/BAM/GFF/VCF> -b <BED/BAM/GFF/VCF> -wa > output


bedtools intersect -abam 3-2-12-1_CAGATC.bam -b cbf1.txt -wa > reads_in_cbf1.bam
samtools sort reads_in_cbf1.bam reads_in_cbf1.sorted
samtools index reads_in_cbf1.sorted.bam reads_in_cbf1.sorted.bam.bai




getIntersectFromBam <- function(){
  
}





