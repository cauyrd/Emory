#usage: python program.py sample.bed tss_file halfextention_size bin_size genome_size fragment/read
import HTSeq
import numpy
import sys 
import re
from operator import itemgetter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from scipy.interpolate import interp1d
import smooth

tssdis = int(sys.argv[2])
width= int(sys.argv[3])
totalbins = tssdis*2/width+1
genomesize = {}
ifp = open(sys.argv[4])
for line in ifp:
	item = line.rstrip().split()
	genomesize[item[0]] = int(item[1])
ifp.close()
x = numpy.arange( -tssdis, tssdis+1)

ifp1 = open(sys.argv[1])
#mydata = []
for myline in ifp1:
	myitem = myline.rstrip().split()
	coverage = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
	bedfile = open(myitem[0])
	if sys.argv[5] == 'fragment':
		while True:
			line1 = bedfile.readline().rstrip()
			if not line1:
				break
			line2 = bedfile.readline().rstrip()
			items1 = line1.split()
			items2 = line2.split()
			if items1[0] == items2[0]:
				chr = items1[0]
				start = min(int(items1[1]), int(items2[1]))
				end = max(int(items1[2]), int(items2[2]))
				iv = HTSeq.GenomicInterval(chr,start,end,'.')
				coverage[ iv ] += 1
	elif sys.argv[5] == 'read':
		for line in bedfile:
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

	tss = set()
	ifp2= open(myitem[1])
	for myline2 in ifp2:
		myitem2 = myline2.rstrip().split()
		chr = myitem2[0]
		pos = int(myitem2[1])
		strand = myitem2[4]
		tss.add(HTSeq.GenomicPosition(chr,pos,strand))
	ifp2.close()
	profile = numpy.zeros(2*tssdis+1, dtype = 'i')
	for p in tss:
		promoter = HTSeq.GenomicInterval(p.chrom, p.pos-tssdis, p.pos+tssdis+1, '.')
		if promoter.start < 0 or promoter.end > genomesize[p.chrom]:
			continue
		wincvg = numpy.fromiter( coverage[promoter], dtype='i', count=2*tssdis+1)
		if p.strand == '+':
			profile += wincvg
		else:
			profile += wincvg[::-1]
#	mydata.append(profile/1e6)
	y = smooth.smooth(profile/float(myitem[-1])*1e6, 1000)
	pyplot.plot(x, y, myitem[2]+'-', lw = 2, label = myitem[3])

	
ifp1.close()
pyplot.legend(prop={'size':8})
pyplot.xlim(-1510,1510)
pyplot.xticks([-1500,-1000,-500,0,500,1000,1500],('-1500','-1000','-500','0','500','1000','1500'))
pyplot.xlabel('Distance form center (bp)')
pyplot.ylabel('Average coverage (reads per million)')
pyplot.title(sys.argv[6])
#pyplot.show()
pyplot.savefig(sys.argv[6]+'.'+sys.argv[5]+'.smooth.png')
pyplot.savefig(sys.argv[6]+'.'+sys.argv[5]+'.smooth.eps')
