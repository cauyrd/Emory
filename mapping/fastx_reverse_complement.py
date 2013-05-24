import gzip,sys
from Bio.Seq import Seq
ifp = gzip.open(sys.argv[1])
ofp = open(sys.argv[1]+'.reverse.fastq','w')
while True:
	line = ifp.readline().rstrip()
	if not line:
		break
	print >> ofp, line
	line = ifp.readline().rstrip()
	seq = Seq(line).reverse_complement()
	print >> ofp, seq
	line = ifp.readline().rstrip()
	print >> ofp, line
	line = ifp.readline().rstrip()
	print >> ofp, line[::-1]
print 'finish reverse complement'
