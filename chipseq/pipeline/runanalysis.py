import os
ifp = open('readme.txt')
pwd = os.getcwd().split('/')
path = '/'.join(pwd[0:(len(pwd)-1)]) + '/mapping/'
igg = '/scratch/bioinfo2/ryang24/irfan/mapping/VCaP_AR_IgG.nodup.bed'
#input = path+'VCaP_AR_Input.nodup.bed'
for line in ifp:
	sample = line.rstrip()
	treat = path+sample+'.nodup.bed'
	os.system('python /compbio/software/pipeline/callpeak.py  hg19 1 '+treat+' '+igg+' '+sample+'_igg 101')
	os.system('sleep 1')
#	os.system('python /compbio/software/pipeline/callpeak.py  hg19 1 '+treat+' '+input+' '+sample+'_input 101')
#	os.system('sleep 1')
