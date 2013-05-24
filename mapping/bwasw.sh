 bwa bwasw -b 9 -q 16 -r 1 -t 6 /compbio/data/bwaIndex/hg19.fa haha.usearch.nr >haha.sam
 samtools view -bS haha.sam > haha.bam
 samtools sort haha.bam haha.sorted
 samtools index haha.sorted.bam
 bamToBed -i haha.sorted.bam >haha.bed

