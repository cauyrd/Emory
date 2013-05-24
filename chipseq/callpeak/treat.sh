bwa aln -t 6 /compbio/data/bwaIndex/ryang/hg19 /scratch/bioinfo2/chipseq-denovo/single-end/rawdata/NRSF/wgEncodeHaibTfbsH1hescNrsfV0416102RawDataRep1.fastq.gz >sample.sai
bwa samse /compbio/data/bwaIndex/ryang/hg19 sample.sai /scratch/bioinfo2/chipseq-denovo/single-end/rawdata/NRSF/wgEncodeHaibTfbsH1hescNrsfV0416102RawDataRep1.fastq.gz >nrsf_h1.sam
python get_flagstat_summary.py nrsf_h1.sam >summary.txt
java -d64 -Xmx81920k -jar /home/ryang24/Software/picard-tools-1.61/MarkDuplicates.jar I=nrsf_h1.sam.sorted.bam O=nrsf_h1.noDup.sam M=nrsf_h1.noDup.txt AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
epd_python get_flagstat_summary.py nrsf_h1.noDup.sam
samtools view -bq 15 nrsf_h1.noDup.sam.sorted.bam > reliable.bam
bamToBed -i reliable.bam > nrsf_h1.noDup.uniq.bed
#/usr/bin/perl
perl /compbio/software/HPeak3/HPeak.pl -sp HUMAN -format bed -t t.inp -n nrsf.h1.rep1 -fmin 150 -fmax 250 
