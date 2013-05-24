from Bio import SeqIO
import sys
ofp = open(sys.argv[1]+'.unique.fa','w')
nameset = set()
for record in SeqIO.parse(sys.argv[1],'fasta'):
	if record.id in nameset:
		continue
	nameset.add(record.id)
	print >> ofp, '>'+record.id
	print >> ofp, record.seq

