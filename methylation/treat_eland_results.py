# the script is used to generate the matrix for constructing Bayesian network. 
# row is all Refseq gene and column is TF, histone mark and DNA methylation.
# usage: python treat_eland_results.py inputfile tss_file halfwinwidth 
import sys
import gzip
import re
import HTSeq
eland_file = gzip.open(sys.argv[1])
output_file = open(sys.argv[1]+'.output','w')
coverage = HTSeq.GenomicArray("auto",stranded=True, typecode="i")
for almnt in eland_file:
    if 'chr' in almnt:
        item = almnt.rstrip().split()
        if '>' in item[0]:
            # this is eland_resutls file column 1 (start at 0) is sequence, column 6 is chr_name, column 7 is start position, column 8 is strand
            if int(item[7]) >= 0:
                if item[8] == 'F':
                    strand = '+'
                else:
                    strand = '-'
                pattern = re.compile('.*(chr.*)\.fa')
                match = pattern.match(item[6])
                read_interval = HTSeq.GenomicInterval(match.group(1),int(item[7]),int(item[7])+len(item[1]),strand)
                coverage[ read_interval ] += 1
        elif 'Run' in item[0]:
            if len(item)<14 and int(item[8]) >= 0:
                strand = '+' if item[9] == 'F' else '-'
                pattern = re.compile('.*(chr.*)\.fa')
                match = pattern.match(item[7])
                read_interal = HTSeq.GenomicInterval(match.group(1),int(item[8]),int(item[8])+len(item[5]),strand)
                coverage[ read_interal ] += 1
        else:
            # this is eland_export file column 6 is sequence, column 8 is chr, column 9 is start position, column 10 is strand
            if len(item)<15 and int(item[9]) >= 0:
                if item[10] == 'F':
                    strand = '+'
                else:
                    strand = '-'
                pattern = re.compile('.*(chr.*)\.fa')
                match = pattern.match(item[8])
                read_interval = HTSeq.GenomicInterval(match.group(1),int(item[9]),int(item[9])+len(item[6]),strand)
                coverage[ read_interval ] += 1

print "complete coverage vector"
#input the refseq gene TSS region
TSS_file = open(sys.argv[2])
tsspos_dict = {} 
for line in TSS_file:
    items = line.rstrip().split()
    mykey = '\t'.join([items[0], items[1], items[2]])
    tsspos_dict[mykey] = items[3:]

print "complete tsspos vector" 
key_list = sorted(tsspos_dict.keys())
halfwinwidth = int(sys.argv[3])
for mykey in key_list:
    items = mykey.split()
    p = HTSeq.GenomicPosition('chr'+items[0],int(items[2]),items[1])
    window = HTSeq.GenomicInterval( p.chrom, p.pos-halfwinwidth, p.pos+halfwinwidth, p.strand)
    count = sum(list(coverage[window]))
    #for iv2,value in coverage[ window ].steps():
    #    count += value
    print >> output_file, float(count),'\t'+p.chrom+'\t'+str(p.start)+'\t'+p.strand+'\t'+'\t'.join(tsspos_dict[mykey])

print "finish all"
