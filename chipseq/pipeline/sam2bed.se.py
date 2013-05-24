import sys
import gzip as gz

ifp = open(sys.argv[1], 'r')
ofp = open(sys.argv[1]+ '.bed', 'w')

while True:
	line = ifp.readline().rstrip()
	if not line:
		break
	if 'XT:A:U' in line:
		if int(line.split()[1]) & 0x10:
			print >> ofp, line.split()[2] + '\t' + str(int(line.split()[3])-1) + '\t'+ str(int(line.split()[3])-1 + len(line.split()[9])) + '\t' + line.split()[0] + '\t' + line.split()[4] + '\t' + '-'
		else:
			print >> ofp, line.split()[2] + '\t' + str(int(line.split()[3])-1) + '\t'+ str(int(line.split()[3])-1 + len(line.split()[9])) + '\t' + line.split()[0] + '\t    ' + line.split()[4] + '\t' + '+'
