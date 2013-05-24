import os
namedict = {}
ifp = open('/scratch/bioinfo2/ryang24/irfan/rawdata/brd4/merge/readme.txt')
for line in ifp:
	item = line.split()
	namedict[item[0]] = item[1]
path = '/scratch/bioinfo2/ryang24/irfan/rawdata/brd4/merge/'
flist = sorted(os.listdir(path))
for i in range(6):
	si = flist[2*i].split('_')[2]
	output = namedict[si]
	os.system('epd_python /compbio/software/pipeline/MapAndDedup.py hg19 1 '+path+flist[2*i]+' '+path+flist[2*i+1]+' '+output)
