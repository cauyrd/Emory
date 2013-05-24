java -d64 -Xmx81920k -jar /home/ryang24/Software/picard-tools-1.61/MarkDuplicates.jar I=MR64.bwa.reverse.sam.sorted.bam O=MR64.bwa.reverse.sam.sorted.noDup.sam M=MR64.bwa.reverse.sam.sorted.noDup.txt AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
python get_flagstat_summary.py MR64.bwa.reverse.sam.sorted.noDup.sam >summary_noDup.txt
