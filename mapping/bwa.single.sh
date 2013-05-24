bwa aln -t 6 /compbio/data/bowtieIndex/hg19.fa rawdata/Foxa1sc101058.fastq.gz >sample1.sai
bwa samse /compbio/data/bowtieIndex/hg19.fa sample1.sai rawdata/Fox1sc101058.fastq.gz >foxa1sc101058.sam
python get_flagstat_summary.py foxa1sc101058.sam >summary1.txt
