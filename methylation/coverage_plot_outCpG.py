import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
import numpy as np
# for LNCaP data, 4 datasets read in +title+outputfile.
col = ['-'+x for x in ['o','^','s','+']]
name = ['5mc','5hmc','M-NGS','M-NGS2']
total_CpG_hg18 = 53478200
total_CpG_hg19 = 53910264
for i in range(4):
	ifp = open(sys.argv[i+1])
	CG_list = []
	for line in ifp:
		CG_list.append(int(line.split()[-4]))
	ifp.close()
	cg_num = []
	for j in range(1,25):
		cg_num.append(sum([x>=j for x in CG_list]))
	if i<2:
		pylab.plot(range(1,25),np.array(cg_num)/float(total_CpG_hg19)*100,col[i],label=name[i])	
	else:
		pylab.plot(range(1,25),np.array(cg_num)/float(total_CpG_hg18)*100,col[i],label=name[i])	

pylab.legend()
pylab.xlabel('read coverage (times)')
pylab.ylabel('% of CpGs')
pylab.title(sys.argv[5])
pylab.savefig(sys.argv[6],format='png')

