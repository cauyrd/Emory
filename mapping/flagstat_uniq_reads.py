# usage: python flagstat_uniq_reads.py samfile_with_dup samfile_no_dup max_insert_size
import sys
import numpy as np
from scipy import stats
def get_uniq_pairs_list(ifp):
	unique_discordant_pairs = set()
	while True:
		line1 = ifp.readline().rstrip()
		if not line1:
			break
		if line1[0] == '@':
			continue
		line2 = ifp.readline().rstrip()
		items = line1.split()
		flag_num = int(items[1])
		if not (flag_num&0x4) and not(flag_num&0x8):
			if 'XT:A:U' in line1 and 'XT:A:U' in line2 and not flag_num&0x2:
				unique_discordant_pairs.add(items[0])
	return unique_discordant_pairs

if __name__ == '__main__':
	ifp = open(sys.argv[1])
	unique_discordant_pairs = get_uniq_pairs_list(ifp)
	ifp = open(sys.argv[2])
	summary = {} 
	isize_set = []
	for line in ifp:
		if line[0] == '@':
			continue
		items = line.rstrip().split()
		flag_num = int(items[1])
		# select all mapped pairs
		if  not (flag_num & 0x4) and not (flag_num & 0x8):
			try:
				summary['mapped_pairs'] += 1 # count the number of reads (if count the number of pairs, divide this number by 2)
			except:
				summary['mapped_pairs'] = 1 
			# select proper pairs
			if flag_num & 0x2:
				try:
					summary['proper_pairs'] += 1 
				except: 
					summary['proper_pairs'] = 1 
				if int(items[8])>0:
					isize_set.append(int(items[8]))
			# select discordant pairs with both ends uniquely mapped
			elif items[0] in unique_discordant_pairs: 
				if items[6] != '=': 
					try:
						summary['diffchr_pairs'] += 1 
					except:
						summary['diffchr_pairs'] = 1 
				elif (abs(int(items[8])) > int(sys.argv[3])):
					try:
						summary['longisize_pairs'] += 1 
					except:
						summary['longisize_pairs'] = 1 
				elif (flag_num & 0x10) and (flag_num & 0x20):
					try:
						summary['discordantRR_pairs'] += 1 
					except KeyError: 
						summary['discordantRR_pairs'] = 1 
				elif not (flag_num & 0x10) and not (flag_num & 0x20):
					try:
						summary['discordantFF_pairs'] += 1 
					except KeyError:
						summary['discordantFF_pairs'] = 1
				elif ((flag_num & 0x10) and (flag_num & 0x40)) or (not(flag_num & 0x10) and (flag_num & 0x80)):
					try:
						summary['discordantRF_pairs'] += 1
					except KeyError:
						summary['discordantRF_pairs'] = 1 
				else:
					try:
						summary['discordantFR_minus_pairs'] += 1 
					except KeyError:
						summary['discordantFR_minus_pairs'] = 1 
			else:
				try:
					summary['repetitive_pairs'] += 1
				except KeyError:
					summary['repetitive_pairs'] = 1
				if flag_num&0x40 and 'XT:A:U' in line:
					try:
						summary['uniquely_mapped_left_reads_remains'] += 1
					except KeyError:
						summary['uniquely_mapped_left_reads_remains'] = 1
				elif flag_num&0x80 and 'XT:A:U' in line:
					try:
						summary['uniquely_mapped_right_reads_remains'] += 1
					except KeyError:
						summary['uniquely_mapped_right_reads_remains'] = 1
		elif (flag_num & 0x8) and not (flag_num & 0x4):
			try:
				summary['singletons'] += 1 
			except KeyError:
				summary['singletons'] = 1
			if flag_num&0x40 and 'XT:A:U' in line:
				try:
					summary['uniquely_mapped_left_reads_remains'] += 1
				except KeyError:
					summary['uniquely_mapped_left_reads_remains'] = 1
			elif flag_num&0x80 and 'XT:A:U' in line:
				try:
					summary['uniquely_mapped_right_reads_remains'] += 1
				except KeyError:
					summary['uniquely_mapped_right_reads_remains'] = 1
		if not (flag_num & 0x4) and flag_num&0x40 and 'XT:A:U' in line:
			try:
				summary['uniquely_mapped_left_reads'] += 1
			except KeyError:
				summary['uniquely_mapped_left_reads'] = 1
		elif not (flag_num & 0x4) and flag_num&0x80 and 'XT:A:U' in line:
			try:
				summary['uniquely_mapped_right_reads'] += 1
			except KeyError:
				summary['uniquely_mapped_right_reads'] = 1
	isize_set = np.array(isize_set)
	bottom = stats.scoreatpercentile(isize_set,5)
	top = stats.scoreatpercentile(isize_set,95)
	isize_set = isize_set[(isize_set>bottom)*(isize_set<top)]
	mu = np.mean(isize_set)
	std = np.std(isize_set)
	for mykey in sorted(summary):
		print mykey+'\t'+str(summary[mykey])
	print 'mean:',mu
	print 'std:',std
