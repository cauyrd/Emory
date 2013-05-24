import sys
ifp = open(sys.argv[1])
for line in ifp:
	items = line.rstrip().split()
	if items[0] == '23':
		chr = 'chrX'
	elif items[0] == '24':
		chr = 'chrY'
	else:
		chr = 'chr'+items[0]
	start = items[1]
	end = items[2]
	print chr+'\t'+start+'\t'+end+'\t'+'\t'.join(items[3:])

