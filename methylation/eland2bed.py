# usage: python eland2bed.py eland_file bed_file 
import sys
import gzip
import re
import HTSeq
eland_file = gzip.open(sys.argv[1])
ofp = open(sys.argv[2],'w')
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
                print >> ofp, match.group(1)+'\t'+str(int(item[7]))+'\t'+str(int(item[7])+len(item[1]))
        elif 'Run' in item[0]:
            if len(item)<14 and int(item[8]) >= 0:
                strand = '+' if item[9] == 'F' else '-'
                pattern = re.compile('.*(chr.*)\.fa')
                match = pattern.match(item[7])
                print >> ofp, match.group(1)+'\t'+str(int(item[8]))+'\t'+str(int(item[8])+len(item[5]))
        else:
            # this is eland_export file column 6 is sequence, column 8 is chr, column 9 is start position, column 10 is strand
            if len(item)<15 and int(item[9]) >= 0:
                if item[10] == 'F':
                    strand = '+'
                else:
                    strand = '-'
                pattern = re.compile('.*(chr.*)\.fa')
                match = pattern.match(item[8])
                print >> ofp, match.group(1)+'\t'+str(int(item[9]))+'\t'+str(int(item[9])+len(item[6]))
