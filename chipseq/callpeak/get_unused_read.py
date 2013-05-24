# usage: program peakfile reliable_bamfile raw_data.fa output.fa
import sys
import os
from Bio import SeqIO
used_set = set()
ifp = open(sys.argv[1])
total = 0
length = 0
read_len = 36
for i,line in enumerate(ifp):
	item = line.rstrip().split()
	region = item[0]+':'+item[1]+'-'+item[2]
	os.system('samtools view -b '+sys.argv[2]+' '+region+' >test.bam')
	os.system('bamToBed -i test.bam >test.bed')
	ifp2 = open('test.bed')
	for i,line in enumerate(ifp2):
		name = line.split()[3]
		used_set.add(name)
	total += i*read_len
	length += int(item[2]) - int(item[1])
	ifp2.close()
ifp.close()
avg = total/float(length)
ofp = open(sys.argv[4],'w')
for record in SeqIO.parse(sys.argv[3],'fasta'):
	if record.id not in used_set:
		print >> ofp, '>'+record.id
		print >> ofp, record.seq
print 'peak average coverage:',avg
