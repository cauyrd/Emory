#usage: plot.py contig_name
import sys, HTSeq
import matplotlib
#matplotlib.use('Agg')
#from matplotlib import pyplot
import pylab
import numpy as np
for read in HTSeq.FastaReader('HAIB_T47D_FoxA1_unused.fasta.contig.merged.fa.coverage_p.fa'):
	if sys.argv[1] in read.name:
		plus = map(int,read.seq.split())
for read in HTSeq.FastaReader('HAIB_T47D_FoxA1_unused.fasta.contig.merged.fa.coverage_p.fa.smoothed.fa'):
	if sys.argv[1] in read.name:
		plus_sm = map(float,read.seq.split())
for read in HTSeq.FastaReader('HAIB_T47D_FoxA1_unused.fasta.contig.merged.fa.coverage_m.fa'):
	if sys.argv[1] in read.name:
		minus = map(int,read.seq.split())
for read in HTSeq.FastaReader('HAIB_T47D_FoxA1_unused.fasta.contig.merged.fa.coverage_m.fa.smoothed.fa'):
	if sys.argv[1] in read.name:
		minus_sm = map(float,read.seq.split())
#all = np.array(minus)+np.array(plus) 
pylab.plot(plus, linestyle='-', color= 'r', label = 'plus')
pylab.plot(minus,linestyle='-', color= 'b', label = 'minus')
pylab.plot(plus_sm, linestyle='-', color= 'g', label = 'plus_sm')
pylab.plot(minus_sm, linestyle='-', color= 'c', label = 'minus_sm')
pylab.title(sys.argv[1])
pylab.legend()
pylab.show()
