import sys
import os
import numpy as np
from Bio import SeqIO
from Bio import Motif
from scipy.stats import fisher_exact
def search_motif(motiflist, seq1, seq2):
	"""search pwm for each motif in the motiflist form sequence"""
	count_all = np.array([[0,0],[0,0]])
	mf = Motif.read(open(motiflist),'jaspar-pfm')
	cutoff = 0.8*mf.max_score()
	for sequence in seq1:
		count_3 = [(pos,score) for pos,score in mf.search_pwm(sequence)]
		max_score3 = max([j for i,j in count_3]) if count_3 else -1e5
		if max_score3 > cutoff:
			count_all[0,0] += 1
		else: 
			count_all[1,0] += 1
	for sequence in seq2:
		count_3 = [(pos,score) for pos,score in mf.search_pwm(sequence)]
		max_score3 = max([j for i,j in count_3]) if count_3 else -1e5
		if max_score3 > cutoff:
			count_all[0,1] += 1
		else: 
			count_all[1,1] += 1
	return count_all

treat = [record.seq for record in SeqIO.parse(sys.argv[1],'fasta')]
control = [record.seq for record in SeqIO.parse(sys.argv[2], 'fasta')]
loc = '/compbio/data/motif/human-mouse/'
ifp = open('/scratch/bioinfo2/chipseq-denovo/single-end/rawdata/matrixlist.txt')
ofp = open(sys.argv[1]+'.motiflist2.txt','w')
for line in ifp:
	name = line.rstrip()
	mat = np.loadtxt(loc+name)
	pfm = np.transpose(mat)
	np.savetxt(name,pfm,fmt='%d')
	mat_all = search_motif(name,treat,control)
	odd, pvalue = fisher_exact(mat_all,alternative='greater')
	print >> ofp, name+'\t'+str(pvalue)
	os.remove(name)
