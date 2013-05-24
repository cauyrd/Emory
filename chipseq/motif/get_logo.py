import sys
import random
import os
import numpy as np
from Bio import SeqIO
from Bio import Motif
from scipy.stats import fisher_exact
loc = '/compbio/data/motif/human-mouse/'
ifp = open('logolist.txt')
for line in ifp:
	name = line.rstrip()
	mat = np.loadtxt(loc+name)
	pfm = np.transpose(mat)
	np.savetxt(name,pfm,fmt='%d')
	mymotif = Motif.read(open(name),'jaspar-pfm')
	mymotif.weblogo(name+'.png')
