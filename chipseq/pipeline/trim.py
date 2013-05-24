#|/usr/bin/env python


import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import stats
import sys
import gzip as gz

ifp = open(sys.argv[1])
#line1 = ifp.readline()
fraglen=[]
#fraglensq=[]
chr = ['chr'+str(i) for i in range(1,23)] + ['chrX']+ ['chrY']
while True:
    line1 = ifp.readline().rstrip()
    if not line1:
        break
    #range1 = line1.split()[1:3]
    line2 = ifp.readline().rstrip()
    #range2 = line2.split()[1:3]
    chrname1=line1.split()[0]
    chrname2=line2.split()[0]
    if chrname1 in chr and chrname1 == chrname2:
        pos = [int(line1.split()[1]),int(line1.split()[2]),int(line2.split()[1]),int(line2.split()[2])]
        dis = max(pos)-min(pos)
        fraglen.append(dis)
    #fraglensq +=[range^2]
ifp.close()
ofp1 = open(sys.argv[1]+'.summary.txt','w')

fraglen=np.array(fraglen)
#fraglen.sort()
chop=float(sys.argv[2])     #### choose a chop percentage, will cut both tails by the chosen percentage
bottom = stats.scoreatpercentile(fraglen, chop)    ### 5th percentile of frag length
top = stats.scoreatpercentile(fraglen, (100-chop))   ### 95th percentile of frag length
mu1 = np.mean(fraglen)                           #### mean before chopping
sigma1 = np.mean(fraglen)		      ### std before chopping
fraglen=fraglen[(fraglen>bottom)*(fraglen<top)]  ### trim off top 5 and bottom 5 

mu = np.mean(fraglen)
sigma = np.std(fraglen)

print >> ofp1, 'Chop by',chop, 'percent top and bottom'
print >> ofp1, 'Before choppoing', 'mean is', mu1, 'and std is', sigma1
print >> ofp1, 'After choppoing' 'mean is', mu, 'and std is', sigma
ofp1.close()

#ofp2=open(sys.argv[1]+'fraglen.txt', 'w')

#for i in fraglen:
#    print >> ofp2, i
#ofp2.close()

#fig = plt.figure()
#plt.hist(fraglen,50)

#fig.savefig(sys.argv[1]+'.png')


#lb= mu-2*sigma
#ub= mu+ 2*sigma

ofp = open(sys.argv[1]+'.trimmed.bed','w')
ifp = open(sys.argv[1], 'r')

while True:
    line1 = ifp.readline().rstrip()
    if not line1:
        break
    line2 = ifp.readline().rstrip()
    chrname1=line1.split()[0]
    chrname2=line2.split()[0]
    if chrname1 in chr and chrname1 == chrname2:
        pos = [int(line1.split()[1]),int(line1.split()[2]),int(line2.split()[1]),int(line2.split()[2])]
        dis = max(pos)-min(pos)
        if dis >= bottom and dis <= top:
            print >> ofp,line1.rstrip()
            print >> ofp,line2.rstrip()


