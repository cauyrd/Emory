import sys
import os

genome = sys.argv[1]
pe = sys.argv[2]
if pe == '1':
	tr1 = sys.argv[3]
	tr2 = sys.argv[4]
	prefix = sys.argv[5]
	os.system('echo "bwa aln -t 8 /compbio/data/bwaIndex/' + genome + '.fa ' + tr1 + ' > ' + prefix + '.1.sai" > ' + prefix + '.MapDedup.sh')
	os.system('echo "bwa aln -t 8 /compbio/data/bwaIndex/' + genome + '.fa ' + tr2 + ' > ' + prefix + '.2.sai" >> ' + prefix + '.MapDedup.sh')
	os.system('echo "bwa sampe /compbio/data/bwaIndex/' + genome + '.fa ' + prefix + '.1.sai ' + prefix + '.2.sai ' + tr1 + ' ' + tr2 + ' > ' + prefix + '.sam" >> ' + prefix + '.MapDedup.sh')
	os.system('echo "python /compbio/software/pipeline/sam2bed.pe.py ' + prefix + '.sam" >> ' + prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/trim.py ' + prefix + '.sam.bed 0.5" >> ' + prefix + '.MapDedup.sh')
	os.system('echo "mv ' + prefix + '.sam.bed.trimmed.bed ' + prefix + '.dup.bed" >> ' + prefix + '.MapDedup.sh')
	os.system('echo "samtools view -bS ' + prefix + '.sam > ' + prefix + '.bam" >> ' + prefix + '.MapDedup.sh') 
	os.system('echo "samtools sort ' + prefix + '.bam ' + prefix + '.sorted" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "java -d64 -Xmx1024m -jar /home/ryang24/Software/picard-tools-1.61/MarkDuplicates.jar I=' + prefix + '.sorted.bam ' + 'O=' + prefix + '.nodup.sam ' + 'M=' + prefix + '.nodup.sum AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true" >> '+ prefix + '.MapDedup.sh') 
	os.system('echo "samtools view -bS ' + prefix + '.nodup.sam ' + '|bamToBed -i stdin > ' + prefix + '.nodup.bed" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "sort -k 4,4 ' + prefix + '.nodup.bed > ' + prefix + '.nodup.sorted" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "python /compbio/software/pipeline/get_pair_end.py ' + prefix + '.nodup.sorted" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/trim.py ' + prefix + '.nodup.sorted.pe.bed 0.5" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "mv ' + prefix + '.nodup.sorted.pe.bed.trimmed.bed  ' + prefix + '.nodup.bed" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/hm_matrix.pe.py ' + prefix + '.nodup.bed /compbio/software/pipeline/' + genome + '.tss 10000 250" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/hm_matrix.pe.py ' + prefix + '.dup.bed /compbio/software/pipeline/' + genome + '.tss 10000 250" >> '+ prefix + '.MapDedup.sh') 
	os.system('echo "gzip ' + prefix + '*.sam ' + prefix + '*.sam.bed" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "rm ' + prefix + '.*.sai ' + prefix + '*.sorted ' + prefix + '*.bam ' + prefix + '*.pe.bed" >> '+ prefix + '.MapDedup.sh')

if pe == '0':
	tr = sys.argv[3]
	prefix = sys.argv[4]
	os.system('echo "bwa aln -t 8 /compbio/data/bwaIndex/' + genome + '.fa ' + tr + ' > ' + prefix + '.sai" > '+ prefix + '.MapDedup.sh')
	os.system('echo "bwa samse /compbio/data/bwaIndex/' + genome + '.fa ' + prefix + '.sai ' +  tr + ' > ' + prefix + '.sam" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "python /compbio/software/pipeline/sam2bed.se.py ' + prefix + '.sam" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "mv ' + prefix + '.sam.bed ' + prefix + '.dup.bed" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "samtools view -bS ' + prefix + '.sam  > ' + prefix + '.bam" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "samtools sort ' + prefix + '.bam ' + prefix + '.sorted" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "java -d64 -Xmx1024m -jar /home/ryang24/Software/picard-tools-1.61/MarkDuplicates.jar I=' + prefix + '.sorted.bam ' + 'O=' + prefix + '.nodup.sam ' + 'M=' + prefix + '.nodup.sum AS=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true" >> '+ prefix + '.MapDedup.sh') 
	os.system('echo "python /compbio/software/pipeline/sam2bed.se.py ' + prefix + '.nodup.sam" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "mv ' + prefix + '.nodup.sam.bed ' + prefix + '.nodup.bed" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/hm_matrix.se.py ' + prefix + '.nodup.bed /compbio/software/pipeline/' + genome + '.tss 5000 50" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "/home/ryang24/Software/epd-7.3-2-rh5-x86_64/bin/python /compbio/software/pipeline/hm_matrix.se.py ' + prefix + '.dup.bed /compbio/software/pipeline/' + genome + '.tss 5000 50" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "gzip ' + prefix + '*.sam" >> '+ prefix + '.MapDedup.sh')
	os.system('echo "rm ' + prefix + '.sai ' + prefix + '*.bam" >> '+ prefix + '.MapDedup.sh')

os.system('qsub -V -S /bin/bash -cwd -N ' + prefix + ' ' + prefix + '.MapDedup.sh')
os.system('sleep 1')


