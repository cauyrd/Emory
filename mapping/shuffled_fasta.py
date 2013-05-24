from Bio import SeqIO
from Bio.Seq import Seq
import sys,random
ofp = open(sys.argv[1]+'.shuffled.fa','w')
length = 0
count = 0
for record in SeqIO.parse(sys.argv[1],'fasta'):
	print >> ofp, '>'+record.id
	length += len(record.seq)
	count += 1
	sequence = list(record.seq)
	random.shuffle(sequence)
	new_seq = ''.join(sequence)
	print >> ofp, new_seq
print length/float(count)
	
