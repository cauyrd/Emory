##### Usage: python callpeak.py genome pe TreatmentData ControlData(if any) Prefix ReadLength 

import sys
import os

genome = sys.argv[1]
pe = sys.argv[2]
tr = sys.argv[3]
if len(sys.argv) == 7:
	cl = sys.argv[4]
	prefix = sys.argv[5]
	rl = sys.argv[6]

if len(sys.argv) == 6:
	cl = None
	prefix = sys.argv[4]
	rl = sys.argv[5]
	

if cl is not None:
	os.system('echo ' + cl + ' > ' + 'c.inp')
	os.system('echo "mkdir ' + prefix + '" > ' + prefix + '.hpeak.sh')
	os.system('echo "cp c.inp ' + prefix + '" >> ' + prefix + '.hpeak.sh')
	os.system('echo "cd ' + prefix + '" >> ' + prefix + '.hpeak.sh')
	os.system('echo "echo ' + tr + ' > chip.inp" >> ' + prefix + '.hpeak.sh')
	if pe == '1':
		os.system('echo "perl /home/prj/compbio/software/HPeak3/HPeak.pl -sp ' + genome + ' -format BED -t chip.inp -c c.inp -n ' + prefix + ' -pe -r ' + rl + ' -ann -wig" >> ' + prefix + '.hpeak.sh')
	else:
		os.system('echo "perl /home/prj/compbio/software/HPeak3/HPeak.pl -sp ' + genome + ' -format BED -t chip.inp -c c.inp -n ' + prefix + ' -r ' + rl + ' -ann -wig" >> ' +prefix + '.hpeak.sh')
else:
	os.system('echo "mkdir ' + prefix + '" > ' + prefix + '.hpeak.sh')
	os.system('echo "cd ' + prefix + '" >> ' + prefix + '.hpeak.sh')
	os.system('echo "echo ' + tr + ' > chip.inp" >> ' + prefix + '.hpeak.sh')
	if pe == '1':
		os.system('echo "perl /home/prj/compbio/software/HPeak3/HPeak.pl -sp ' + genome + ' -format BED -t chip.inp -n ' + prefix + ' -pe -r ' + rl + ' -ann -wig" >> ' + prefix + '.hpeak.sh')
	else:
		os.system('echo "perl /home/prj/compbio/software/HPeak3/HPeak.pl -sp ' + genome + ' -format BED -t chip.inp -n ' + prefix + ' -r ' + rl + ' -ann -wig" >> ' + prefix + '.hpeak.sh')

os.system('echo "perl /home/prj/compbio/software/utl/disTSS.pl human 10000 ' + prefix + '.hpeak.out /home/prj/compbio/data/annotations/' + genome + '/' + genome + '.refseq.unique.TSS ' + prefix + '.dis" >> ' + prefix + '.hpeak.sh')
os.system('qsub -V -S /bin/bash -cwd -N ' + prefix + ' ' + prefix + '.hpeak.sh')


