import sys
import gzip as gz

ifp = open(sys.argv[1], 'r')
ofp = open(sys.argv[1]+ '.bed', 'w')

while True:
	line1 = ifp.readline().rstrip()
	if not line1:
		break
	if line1.split()[0] == '@SQ':
		continue
	line2 = ifp.readline().rstrip()

	if 'XT:A:U' in line1 and 'XT:A:U'in line2:
		if int(line1.split()[1]) & 0x10:
			strand1 = '-'
		else:
			strand1 = '+'
		if int(line2.split()[1]) & 0x10:
			strand2 = '-'
		else:
			strand2 = '+'
		if int(line1.split()[1]) & 0x40:
			pair1 = '1'
			pair2 = '2'
		else:
			pair1 = '2'
			pair2 = '1'
		print >> ofp, line1.split()[2] + '\t' + str(int(line1.split()[3])-1) + '\t'+ str(int(line1.split()[3])-1 + len(line1.split()[9])) + '\t' + line1.split()[0] + '/' + pair1 + '\t' + line1.split()[4] + '\t' + strand1
		print >> ofp, line2.split()[2] + '\t' + str(int(line2.split()[3])-1) + '\t'+ str(int(line2.split()[3])-1 + len(line2.split()[9])) + '\t' + line2.split()[0] + '/' + pair2 + '\t' + line2.split()[4] + '\t' + strand2
