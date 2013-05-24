'''
File: generate_bins.py
Author: Rendong Yang
Description: the program is used to divide human genome into bins, bin length is specified by user
Usage: generate_bins.py chromosomal_end_file bin_length
'''
import sys
nchr = 24
CHROMOFILE = open(sys.argv[1])
resolution = int(sys.argv[2])
chromosize = [0]*nchr
for i,line in enumerate(CHROMOFILE):
	chromosize[i] = int(line.rstrip())
ofp = open('genome_bin.bed','w')
for i in range(nchr):
	upperlimit = chromosize[i]/resolution
	for j in range(upperlimit):
		start = resolution*j
		end = start + resolution
		if (i == (nchr-2)):
			print >> ofp, 'chrX\t'+str(start)+'\t'+str(end)
		elif (i == (nchr-1)):
			print >> ofp, 'chrY\t'+str(start)+'\t'+str(end)
		else:
			print >> ofp, 'chr'+str(i+1)+'\t'+str(start)+'\t'+str(end)
