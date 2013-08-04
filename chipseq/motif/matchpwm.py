# usage: matchpwm.py sequence control(optional)
import sys
import os
import random
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
		count = [(pos,score) for pos,score in mf.search_pwm(sequence)]
		max_score = max([j for i,j in count]) if count else -1e5
		if max_score > cutoff:
			count_all[0,0] += 1
		else: 
			count_all[1,0] += 1
	for sequence in seq2:
		count = [(pos,score) for pos,score in mf.search_pwm(sequence)]
		max_score3 = max([j for i,j in count]) if count else -1e5
		if max_score > cutoff:
			count_all[0,1] += 1
		else: 
			count_all[1,1] += 1
	return count_all

def get_control(inputfile):
	ofp = open(inputfile+'.shuffled.fa','w')
	for record in SeqIO.parse(inputfile,'fasta'):
		print >> ofp, '>'+record.id
		sequence = list(record.seq)
		random.shuffle(sequence)
		new_seq = ''.join(sequence)
		print >> ofp, new_seq
	ofp.close()
	return

treat = [record.seq for record in SeqIO.parse(sys.argv[1],'fasta')]
if (len(sys.argv)<3):
	get_control(sys.argv[1])
	control = [record.seq for record in SeqIO.parse(sys.argv[1]+'.shuffled.fa', 'fasta')]
else:
	control = [record.seq for record in SeqIO.parse(sys.argv[2], 'fasta')]

loc = '/compbio/data/motif/human-mouse/'
ifp = open(loc+'matrixlist.txt')
ofp = open(sys.argv[1]+'.motiflist.txt','w')
for line in ifp:
	name = line.rstrip()
	mat = np.loadtxt(loc+name)
	pfm = np.transpose(mat)
	np.savetxt(sys.argv[1]+'.'+name,pfm,fmt='%d')
	mat_all = search_motif(sys.argv[1]+'.'+name,treat,control)
	odd, pvalue = fisher_exact(mat_all,alternative='greater')
	motifname = name.split('.')[0]
	try:
		tmp = open(loc+motifname+'.consensus')
		seq = tmp.readline().rstrip()
	except IOError:
		seq = 'unknown'
	print >> ofp, motifname+'\t'+seq+'\t'+str(pvalue)
	os.remove(sys.argv[1]+'.'+name)
	tmp.close()
