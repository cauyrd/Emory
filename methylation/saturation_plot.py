import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
import numpy as np
import pybedtools
# for LNCaP data, 4 datasets read in.
col = ['-'+x for x in ['o','^','s','+']]
name = ['5mc','5hmc','M-NGS','M-NGS2']
for i in range(4):
	print 'i=',i
	read_file = pybedtools.BedTool(sys.argv[i+1])
	cg_num = []
	if i == 0 or i==1:
		cg = pybedtools.BedTool('/scratch/bioinfo2/Jindan/job052812/coverage_plot/hg19_cg_position.bed')
	else:
		cg = pybedtools.BedTool('/scratch/bioinfo2/Jindan/job052812/coverage_plot/hg18_cg_position.bed')
	for rd_num in range(1,26):
		rd = pybedtools.BedTool(read_file[:rd_num*1e6])
		cg_rd = rd.coverage(cg) 
		total = sum([int(x[4])>=1 for x in cg_rd])
		print total
		cg_num.append(total)
	pylab.plot(range(1,26),np.log10(np.array(cg_num)),col[i],label=name[i])	
pylab.legend()
pylab.xlabel('number of millions reads')
pylab.ylabel('number of CG (log10)')
pylab.title('LNCaP')
pylab.savefig('LNCaP.coverage.png',format='png')

