 bwa aln -t 6 /compbio/data/bowtieIndex/hg19.fa MR64_R1.fastq.gz > sample2R1.sai 
 bwa aln -t 6 /compbio/data/bowtieIndex/hg19.fa MR64_R2.fastq.gz > sample2R2.sai
 bwa sampe /home/mzhao6/data/bowtieIndex/hg19.fa sample2R1.sai sample2R2.sai MR64_R1.fastq.gz MR64_R2.fastq.gz > MR64.bwa.sam

