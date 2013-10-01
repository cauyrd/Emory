# usage: count_motif.py motif.pwm sequence.fa
import sys
import os
import random
import numpy as np
from Bio import SeqIO
from Bio import Motif
from scipy.stats import fisher_exact

def search_motif(mf,seq):
	"""search pwm for each motif in the motiflist form sequence"""
	cutoff = 0.8*mf.max_score()
	result = [(score,pos) for pos,score in mf.search_pwm(seq)]
	if not result:
		return None
	scores = [item[0] for item in result]
	pos = [item[1] for item in result]
	if max(scores) > cutoff:
		return pos[scores.index(max(scores))]
	else: 
		return None

count = 0
mf = Motif.read(open(sys.argv[1]),'jaspar-pfm')
ofp = open(sys.argv[2]+'.motif.fa','w')
for i,record in enumerate(SeqIO.parse(sys.argv[2],'fasta')):
	hit = search_motif(mf,record.seq) 
	if hit == None:
		continue
	else:
		record.id = record.id+'_'+str(hit)
		print >> ofp, '>'+record.id
		print >> ofp, record.seq
		count += 1
print str(i+1)+'\ttotal seq(s)'
print str(count)+'\tcontains motifs ('+str(count/float(i+1)*100)+'%)'
