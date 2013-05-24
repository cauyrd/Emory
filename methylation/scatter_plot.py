'''
File: select_bin.py
Author: Rendong Yang
Description: this porgram is used to generate two vectors in which each item >= threshold.
Usage: select_bin.py x_axis_file y_axis_file threshold output_figure_name
'''
import sys, itertools
from matplotlib import pylab
from scipy.stats.stats import pearsonr
import numpy as np
ifp1 = open(sys.argv[1])
ifp2 = open(sys.argv[2])
threshold = int(sys.argv[3])
x=[]
y=[]
for line1, line2 in itertools.izip(ifp1,ifp2):
	items1 = line1.split()
	items2 = line2.split()
	if int(items1[3])>=threshold and int(items2[3])>=threshold:
		x.append(int(items1[3]))
		y.append(int(items2[3]))
x_array = np.log10(np.array(x))
y_array = np.log10(np.array(y))
r_square = np.square(pearsonr(x_array,y_array)[0])
pylab.plot(x_array,y_array,'o')
pylab.xlabel('LNCaP (M-NGS) read counts log10')
pylab.ylabel('LNCaP-2 (M-NGS) read counts log10')
pylab.title('R square = '+str(r_square))
pylab.savefig(sys.argv[4],format='png')
