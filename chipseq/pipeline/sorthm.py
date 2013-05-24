import sys
import os
import numpy as np
ifp = open(sys.argv[1], 'r')

PATH = os.getcwd() + '/'
k27 = open(PATH + 'K27.bed.hm', 'r')
wt = open(PATH + 'Eed.WT.RNF2.bed.hm', 'r')
ko = open(PATH + 'Eed.KO.RNF2.uniq.bed.hm', 'r')

D = {}
for line in ifp:
	#D[line.split()[1]] =str( -float(line.split()[0]))
	D[line.split()[1]] = line.split()[0]

h1 = k27.readline().rstrip()
h2 = wt.readline().rstrip()
h3 = ko.readline().rstrip()
ncol = len(h1.split())
if sys.argv[2] == 'sec1.data':
	ofp1= open('header', 'w')
#	print >> ofp1, 'pvalue' + '\t' + h1 + '\t' + '\t' + '\t'.join(h2.split()[2:ncol]) + '\t' + '\t' + '\t'.join(h3.split()[2:ncol])
	for i in range(35):
		print >> ofp1, 'pvalue' + '\t' + h1.split()[0] + '\t' + h1.split()[1] + '\t' + '\t' + '\t' + '\t' + '\t'.join(h1.split()[2:ncol]) + '\t' + '\t' + '\t' + '\t' + '\t'.join(h2.split()[2:ncol]) + '\t' + '\t' + '\t' + '\t'  + '\t'.join(h3.split()[2:ncol]) +'\t' + 'blank' + '\t' + 'blank'
	ofp1.close()

ofp2 = open(sys.argv[2] , 'w')
while True:
	h1 = k27.readline().rstrip()
	if not h1:
		break
	#ncol = len(h1.split())
	h2 = wt.readline().rstrip()
	h3 = ko.readline().rstrip()
	if h1.split()[1] in D:
		print >> ofp2, D[h1.split()[1]] + '\t' + h1.split()[0] + '\t' + h1.split()[1] + '\t' + '\t' + '\t' + '\t' + '\t'.join(map(str,np.array(map(int,h1.split()[2:ncol]))*(-1))) + '\t' + '\t' + '\t' + '\t' + '\t'.join(h2.split()[2:ncol]) + '\t' + '\t' + '\t' + '\t'  + '\t'.join(h3.split()[2:ncol]) + '\t' + '\t' + '\t' + '\t'
	#else:
	#	print >> ofp2, '0' + '\t' + h1 + '\t' + '\t' + '\t'.join(h2.split()[2:ncol]) + '\t' + '\t' + '\t'.join(h3.split()[2:ncol])
for i in range(35):
	print >> ofp2, '0' + '\t' + 'blank' + '\t' + 'tmp' + '\t' + '\t' + '\t' + '\t' + '\t'.join(['\t']*(ncol-2)) + '\t' + '\t' + '\t' + '\t' + '\t'.join(['\t']*(ncol-2)) + '\t' + '\t' + '\t' + '\t'  + '\t'.join(['\t']*(ncol-2)) +'\t' + '\t' + '\t' + '\t'
ofp2.close()

os.system('sort -n -k 1,1 ' + sys.argv[2] + ' > ' + sys.argv[2] + '.sorted')

if sys.argv[2] == 'sec1.data':
	os.system('cat header ' + sys.argv[2] + '.sorted > '+ sys.argv[2] + '.tmp')
	os.system('mv ' + sys.argv[2] +'.tmp ' + sys.argv[2] +'.sorted')

