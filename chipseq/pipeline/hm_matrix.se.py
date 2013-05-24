#usage: python hm_matrix.py sample.bed tss_file
import HTSeq
import numpy
import sys 
import re
def get_hm_bin(coverage, promoter, totalbins):
	hm_list = []
	"""this is to cut promoter into totalbins bins and count the coverage in each bin"""
	bins = numpy.linspace(promoter.start, promoter.end, totalbins+1)
	for i in range(totalbins):
		bin_range = HTSeq.GenomicInterval(promoter.chrom, int(bins[i]), int(bins[i+1]), '.')
		hm_list.append(int(sum(numpy.fromiter(coverage[bin_range], dtype='i', count=(int(bins[i+1])-int(bins[i]))))/float(reads/1e6)))
#		hm_list.append(numpy.fromiter(coverage[bin_range], dtype='i', count=(int(bins[i+1])-int(bins[i]))).sum()/float(bin_range.length))
#	if sum(hm_list):
	#	return hm_list/sum(hm_list)
	#else:
	return hm_list

coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
bedfile = open(sys.argv[1])
reads = 0
for line in bedfile:
	items = line.rstrip().split()
	reads += 1
	chr = items[0]
	start = int(items[1])
	end = int(items[2])
	strand = items[-1]
	iv = HTSeq.GenomicInterval(chr,start,end,strand)
	coverage[ iv ] += 1
bedfile.close()


tss = set()
ifp= open(sys.argv[2])
ofp = open(sys.argv[1] + '.hm', 'w')
D={}

for line in ifp:
	items = line.rstrip().split()
	chr = items[0]
	pos = int(items[2])
	strand = items[1]
	gene = items[-1]
	tss.add(HTSeq.GenomicPosition(chr,pos,strand))
	id = str(chr) + ':' + str(pos) + ':' + strand
	D[id] = gene

tssdis = int(sys.argv[3])
width= int(sys.argv[4])
totalbins = tssdis*2/width+1
f = re.compile('(.*).bed')
file = f.match(sys.argv[1])
label= [file.group(1) + '.bin' + str(i) for i in range(-tssdis/width,tssdis/width+1)]
colname = 'UNIQID' + '\t' + 'NAME' + '\t' + '\t'.join(label)
print >> ofp, colname
#profile1 = numpy.zeros(2*tssdis, dtype = 'i')
#profile2 = numpy.zeros(2*tssdis, dtype = 'i')
#profile3 = numpy.zeros(2*tssdis, dtype = 'i')

for p in tss:
	promoter = HTSeq.GenomicInterval(p.chrom, p.pos-tssdis-1, p.pos+tssdis-1, '.')
	if promoter.start < 0:
		continue
	wincvg = get_hm_bin(coverage,promoter,totalbins)
	id = str(p.chrom) + ':' + str(p.pos) + ':' + p.strand
	if p.strand == '+':
		print >> ofp, id + '\t' + D[id] + '\t' + '\t'.join(map(str, wincvg))
	else:
		print >> ofp, id + '\t' + D[id] + '\t' + '\t'.join(map(str,wincvg[::-1]))
    
