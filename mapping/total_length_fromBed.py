import sys
ifp = open(sys.argv[1])
total_len = 0
for line in ifp:
	items = line.split()
	length = int(items[2]) - int(items[1])
	total_len += length
total_len = float(total_len)/1e6
print "total length (Mb): ", total_len
print "percentage of human genome: ",total_len/3098.0 
