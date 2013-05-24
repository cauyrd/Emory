'''
File: barplot.py
Author: Rendong Yang
Description: program bin_cg 5mc 5hmc M-NGS
'''
import sys, matplotlib, itertools
import numpy as np
matplotlib.use('Agg')
from matplotlib import pyplot
def define_GC_content(num):
	"""classfy GC_content into low, medium and high based on number of CpGs"""
	if num>=1 and num<10:
		return 'low'
	elif num>=10 and num<18:
		return 'medium'
	elif num>=18:
		return 'high'
	else:
		return 'none'
ifp1 = open(sys.argv[1])
ifp2 = open(sys.argv[2])
ifp3 = open(sys.argv[3])
ifp4 = open(sys.argv[4])
cg_dict = {'low':[0.0,0.0,0.0,0], 'medium':[0.0,0.0,0.0,0], 'high':[0.0,0.0,0.0,0], 'none':[0.0,0.0,0.0,0]} 
for line1, line2, line3, line4 in itertools.izip(ifp1,ifp2,ifp3,ifp4):
	label = define_GC_content(int(line1.split()[3])/2)
	cg_dict[label][0] += int(line2.split()[3])
	cg_dict[label][1] += int(line3.split()[3])
	cg_dict[label][2] += int(line4.split()[3])
	cg_dict[label][3] += 1
rpm = []
print cg_dict
for i in range(3):
	rpm.append([])
	for label in ['low','medium','high']:
		rpm[i].append(cg_dict[label][i]/cg_dict[label][3])
ind = np.arange(1,6,2)
width = 0.35
fig = pyplot.figure()
ax = fig.add_subplot(111)
for i,j,k in itertools.izip(range(3),['r','g','b'],['5mc','5hmc','M-NGS']):
	ax.bar(ind+width*i, rpm[i], width, color = j, label=k)
	ax.legend(prop={'size':10})
ax.set_ylabel('average rpm')
ax.set_xticks(ind+width)
ax.set_xticklabels(('low','medium','high'))
pyplot.savefig('cpg.png',format='png')
