# Variable definition
REFERENCE=/compbio/data/bwaIndex/ryang/hg19
READ_FILE1=/bioinfo/data/wuyimi/RawData/mctp_SI_4609_C07ECACXX_8_trim_1.fq  
READ_FILE2=/bioinfo/data/wuyimi/RawData/mctp_SI_4609_C07ECACXX_8_trim_2.fq 
SAM_PREFIX=foxa1_lncap_dht_pe
READ_LEN=100

# BWA mapping pair-end reads
bwa aln -t 6 $REFERENCE $READ_FILE1>sample1.sai
bwa aln -t 6 $REFERENCE $READ_FILE2>sample2.sai
bwa sampe $REFERENCE sample1.sai sample2.sai $READ_FILE1 $READ_FILE2 > $SAM_PREFIX.sam

# SAM to BAM format
samtools view -bS $SAM_PREFIX.sam > $SAM_PREFIX.bam
samtools sort $SAM_PREFIX.bam $SAM_PREFIX.sorted
samtools index $SAM_PREFIX.sorted.bam
samtools flagstat $SAM_PREFIX.sorted.bam > $SAM_PREFIX.summary.txt
rm $SAM_PREFIX.bam $SAM_PREFIX.sam

# Remove PCR duplicate reads
java -d64 -Xmx81920k -jar /home/ryang24/Software/picard-tools-1.61/MarkDuplicates.jar I=$SAM_PREFIX.sorted.bam O=$SAM_PREFIX.noDup.sam M=$SAM_PREFIX.noDup.txt AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
samtools view -bS $SAM_PREFIX.noDup.sam > $SAM_PREFIX.noDup.bam
samtools sort $SAM_PREFIX.noDup.bam $SAM_PREFIX.noDup.sorted
samtools index $SAM_PREFIX.noDup.sorted.bam
rm $SAM_PREFIX.noDup.bam $SAM_PREFIX.noDup.sam

# Get both end uniquely mapped reads
samtools view -bq 10 $SAM_PREFIX.noDup.sorted.bam > $SAM_PREFIX.unique.bam
samtools index $SAM_PREFIX.unique.bam
samtools flagstat $SAM_PREFIX.unique.bam > $SAM_PREFIX.unique.txt
bamToBed -i $SAM_PREFIX.unique.bam > $SAM_PREFIX.unique.bed
sort -k4,4 $SAM_PREFIX.unique.bed > $SAM_PREFIX.unique.sorted
epd_python /scratch/bioinfo2/ryang24/NGStoolbox/chipseq/pipeline/get_pair_end.py $SAM_PREFIX.unique.sorted
epd_python /scratch/bioinfo2/ryang24/NGStoolbox/chipseq/pipeline/trim.py $SAM_PREFIX.unique.sorted.pe.bed 0.5
mv $SAM_PREFIX.unique.sorted.pe.bed.trimmed.bed $SAM_PREFIX.unique.bed
rm $SAM_PREFIX.unique.sorted.pe.bed

# Generate coverage TDF file for IGV
igvtools count $SAM_PREFIX.unique.bam $SAM_PREFIX.unique.cov.tdf hg19

# Peak calling using HPeak
echo $SAM_PREFIX.unique.bed > $SAM_PREFIX.treat.inp
perl /compbio/software/HPeak3/HPeak.pl -sp HUMAN -format BED -t $SAM_PREFIX.treat.inp -n $SAM_PREFIX.peak -pe -r $READ_LEN -ann -wig
epd_python /scratch/bioinfo2/ryang24/NGStoolbox/chipseq/callpeak/peakToBed.py $SAM_PREFIX.peak.hpeak.out > $SAM_PREFIX.peak.hpeak.bed
