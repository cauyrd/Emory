#program peakfile orderfile
import sys
peakdict = {}
ifp = open(sys.argv[1])
for line in ifp:
	item = line.rstrip().split()
	peakdict[item[3]] = line.rstrip()
ifp.close()
ifp = open(sys.argv[2])
ofp = open(sys.argv[1]+'.reorder.bed','w')
for line in ifp:
	print >> ofp,peakdict[line.rstrip()]
ifp.close()
ofp.close()
	

