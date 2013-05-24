import sys
ifp = open(sys.argv[1])
orientation = {'RR':0, 'RF':0, 'FF':0, 'FR':0}
total = 0
for line in ifp:
	if '@' == line[0]:
		continue
	total += 1
	items = line.split()
	flag_num = int(items[1])
	if (flag_num & 0x10) and (flag_num & 0x20):
		orientation['RR'] += 1
	elif (flag_num & 0x10) and not (flag_num & 0x20):
		orientation['FF'] += 1
	elif ((flag_num & 0x10) and (flag_num & 0x40)) or (not(flag_num & 0x10) and (flag_num & 0x80)):
		orientation['RF'] += 1
	else:
		orientation['FR'] += 1
for each in orientation:
	print each+':'+str(orientation[each])+'\t'+'percentage:'+str(orientation[each]*100/float(total))
