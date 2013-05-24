# usage: slopbed.py input.bed genome.txt upper_kb down_kb
import sys
upper = 5000
down = 5000
if len(sys.argv) > 3:
	upper = int(sys.argv[3])
	down = int(sys.argv[4])
chr_name = ['chr'+str(i) for i in range(1,23)]
chr_name = chr_name+['chrX','chrY']
genome = dict.fromkeys(chr_name)
ifp = open(sys.argv[2])
for line in ifp:
	item = line.rstrip().split()
	try:
		genome[item[0]] = int(item[1])
	except KeyError:
		pass
ifp.close()
ifp = open(sys.argv[1])
ofp = open(sys.argv[1]+'.treated.bed','w')
for line in ifp:
	item = line.rstrip().split()
	if item[0] not in chr_name:
		continue
	if item[-1] == '-':
		item[1] = int(item[1]) + upper
		if item[1] > genome[item[0]]:
			item[1] = genome[item[0]]
		item[2] = int(item[2]) - down
		if item[2] < 0:
			item[2] = 0
		new_line = item[0]+'\t'+str(item[2])+'\t'+str(item[1])+'\t'+item[3]+'\t'+item[4]
	else:
		item[1] = int(item[1]) - upper
		if item[1] < 0:
			item[1] = 0
		item[2] = int(item[2]) + down
		if item[2] > genome[item[0]]:
			item[2] = genome[item[0]]
		new_line = item[0]+'\t'+str(item[1])+'\t'+str(item[2])+'\t'+item[3]+'\t'+item[4]
	print >> ofp, new_line
