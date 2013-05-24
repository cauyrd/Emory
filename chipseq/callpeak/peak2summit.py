#input is hpeak.out file
import sys
import numpy as np
ifp = open(sys.argv[1])
ofp = open(sys.argv[1]+'.summit.bed','w')
for i,line in enumerate(ifp):
	item = line.rstrip().split()
	chr = item[0]
	start = int(item[1])+int(np.round(float(item[4])))
	end = start+1
	name = 'peak'+str(i)
	strand = '+'
	print >> ofp, chr+'\t'+str(start)+'\t'+str(end)+'\t'+name+'\t'+strand
