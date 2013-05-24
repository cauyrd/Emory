#usage: python hm_matrix.py sample.bed tss_file halfextention_size bin_size genome_size fragment/read
import HTSeq
import numpy
import sys 
import re
from operator import itemgetter
def get_hm_bin(coverage, promoter, totalbins):
	hm_list = []
	"""this is to cut promoter into totalbins bins and count the coverage in each bin"""
	bins = numpy.linspace(promoter.start, promoter.end, totalbins+1)
	for i in range(totalbins):
		bin_range = HTSeq.GenomicInterval(promoter.chrom, int(bins[i]), int(bins[i+1]), '.')
		hm_list.append(int(sum(numpy.fromiter(coverage[bin_range], dtype='i', count=(int(bins[i+1])-int(bins[i]))))/(reads/1e6)))
	return hm_list

coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
bedfile = open(sys.argv[1])
reads = 0
if sys.argv[6] == 'fragment':
	while True:
		line1 = bedfile.readline().rstrip()
		if not line1:
			break
		reads += 1
		line2 = bedfile.readline().rstrip()
		items1 = line1.split()
		items2 = line2.split()
		if items1[0] == items2[0]:
			chr = items1[0]
			start = min(int(items1[1]), int(items2[1]))
			end = max(int(items1[2]), int(items2[2]))
			iv = HTSeq.GenomicInterval(chr,start,end,'.')
			coverage[ iv ] += 1
elif sys.argv[6] == 'read':
	for line in bedfile:
		reads += 1
		item = line.rstrip().split()
		chr = item[0]
		start = int(item[1])
		end = int(item[2])
		iv = HTSeq.GenomicInterval(chr,start,end,'.')
		coverage[ iv ] += 1
else:
	print "specify using fragment or read as counting measurment."
	sys.exit(0)
bedfile.close()

tss = []
ifp= open(sys.argv[2])
ofp = open(sys.argv[1] + '.' + sys.argv[6] + '.'+sys.argv[3]+'.common.nosort.cdt', 'w')
D={}

for line in ifp:
	items = line.rstrip().split()
	chr = items[0]
	pos = int(items[1])
	strand = items[4]
	gene = items[3]
	tss.append(HTSeq.GenomicPosition(chr,pos,strand))
	id = str(chr) + ':' + str(pos) + ':' + strand
	D[id] = gene
ifp.close()
tssdis = int(sys.argv[3])
width= int(sys.argv[4])
totalbins = tssdis*2/width+1
#f = re.compile('(.*).bed')
#file = f.match(sys.argv[1])
label= ['bin' + str(i) for i in range(-tssdis/width,tssdis/width+1)]
colname = 'UNIQID' + '\t' + 'NAME' + '\t' + '\t'.join(label)
print >> ofp, colname

genomesize = {}
ifp = open(sys.argv[5])
for line in ifp:
	item = line.rstrip().split()
	genomesize[item[0]] = int(item[1])
ifp.close()

for p in tss:
	id = str(p.chrom) + ':' + str(p.pos) + ':' + p.strand
	promoter = HTSeq.GenomicInterval(p.chrom, p.pos-tssdis, p.pos+tssdis+1, '.')
	if promoter.start < 0 or promoter.end > genomesize[p.chrom]:
		continue
	wincvg = get_hm_bin(coverage,promoter,totalbins)
	if p.strand == '+':
		print >> ofp, id + '\t' + D[id] + '\t' + '\t'.join(map(str, wincvg))
	else:
		print >> ofp, id + '\t' + D[id] + '\t' + '\t'.join(map(str,wincvg[::-1]))
