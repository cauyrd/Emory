#usage: sam2fastx.py -x[a--fasta|q--fastq] your.sam
import sys
type = sys.argv[1]
ifp = open(sys.argv[2])
out = sys.argv[2]
if type == '-a':
	ofp = open(sys.argv[2]+'.fa','w')
	for line in ifp:
		if '@' == line[0]:
			continue
		items = line.split()
		if int(items[1]) & 0x1:
			if int(items[1])&0x40:
				names = items[0]+'/1'
			else:
				names = items[0]+'/2'
		else:
			name = items[0]
		seq = items[9]
		print >> ofp, '>'+names
		print >> ofp, seq
elif type == '-q':
	ofp = open(sys.argv[2]+'.fq','w')
	for line in ifp:
		if '@' == line[0]:
			continue
		items = line.split()
		if int(items[1]) & 0x1:
			if int(items[1])&0x40:
				names = items[0]+'/1'
			else:
				names = items[0]+'/2'
		else:
			name = items[0]
		seq = items[9]
		qual = items[10]
		print >> ofp, '@'+names
		print >> ofp, seq
		print >> ofp, '+'
		print >> ofp, qual
else:
	print 'input wrong file type. Exit!'
	exit(0)

