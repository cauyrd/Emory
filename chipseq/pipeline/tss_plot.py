#usage: python tss_plot.py groseq.plus.bed groseq.minus.bed methy.bed tss_file output_image_file_name
import HTSeq
import numpy
import sys 
from matplotlib import pyplot
coverage1 = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
bedfile = open(sys.argv[1])
for line in bedfile:
	items = line.split()
	chr = items[0]
	start = int(items[1])
	end = start+len(items[3])
	strand = items[5]
	iv = HTSeq.GenomicInterval(chr,start,end,strand)
	coverage1[ iv ] += 1
bedfile.close()
coverage2 = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
bedfile = open(sys.argv[2])
for line in bedfile:
	items = line.split()
	chr = items[0]
	start = int(items[1])
	end = start+len(items[3])
	strand = items[5]
	iv = HTSeq.GenomicInterval(chr,start,end,strand)
	coverage2[ iv ] += 1
bedfile.close()
coverage3 = HTSeq.GenomicArray('auto',stranded=False,typecode='i')
bedfile = open(sys.argv[3])
for line in bedfile:
	items = line.split()
	chr = items[0]
	start = int(items[1])
	end = int(items[2])
	iv = HTSeq.GenomicInterval(chr,start,end,'.')
	coverage3[ iv ] += int(float(items[4])*1e6)
bedfile.close()
tsspos = set()
ifp = open(sys.argv[4])
for line in ifp:
	items = line.split()
	chr = items[0]
	pos = (int(items[1])+int(items[2]))/2
	strand = items[3]
	tsspos.add(HTSeq.GenomicPosition(chr,pos,strand))
halfwinwidth = 5000
profile1 = numpy.zeros(2*halfwinwidth, dtype = 'i')
profile2 = numpy.zeros(2*halfwinwidth, dtype = 'i')
profile3 = numpy.zeros(2*halfwinwidth, dtype = 'i')
for p in tsspos:
	window = HTSeq.GenomicInterval(p.chrom, p.pos-halfwinwidth, p.pos+halfwinwidth, '.')
	if window.start < 0:
		continue
	wincvg1 =  numpy.fromiter( coverage1[window], dtype='i', count=2*halfwinwidth )
	wincvg2 =  numpy.fromiter( coverage2[window], dtype='i', count=2*halfwinwidth )
	wincvg3 =  numpy.fromiter( coverage3[window], dtype='i', count=2*halfwinwidth )
	if p.strand == '+':
		profile1 += wincvg1
		profile2 += wincvg2
		profile3 += wincvg3
	else:
		profile1 += wincvg2[::-1]
		profile2 += wincvg1[::-1]
		profile3 += wincvg3[::-1]
x = numpy.arange( -halfwinwidth, halfwinwidth)
pyplot.subplot(2,1,1)
pyplot.plot(x, profile1, 'r-', label = 'sense')
pyplot.plot(x, profile2, 'b-', label = 'antisense')
pyplot.legend()
pyplot.subplot(2,1,2)
pyplot.plot(x, profile3/1e6, 'g-', label = 'bisulfite' )
pyplot.legend()
#pyplot.show()
pyplot.savefig(sys.argv[5],format='png')
