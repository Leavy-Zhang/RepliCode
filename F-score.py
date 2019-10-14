#! /usr/bin/python
import re,sys,os
from math import log
import numpy as np
libraryListSize=[float(j) for j in sys.argv[1].split(',')]
#libraryListSize=[]
#for j in sys.argv[1].split(','):
#	x=os.popen('samtools flagstat '+j).read()
#	n=0
#	for i in x.split('\n'):
#		if n==4:
#			libraryListSize.append(float(i.split(' ')[0]))
#			break
#		else:
#			n+=1

#scaleFactor=[sum(libraryListSize)/i for i in libraryListSize]
scaleFactor=list(min(libraryListSize)/np.array(libraryListSize))
print scaleFactor

f=open(sys.argv[4],'w')
treat=open(sys.argv[2]).readlines()
Input=open(sys.argv[3]).readlines()
for i in range(len(treat)):
	tLine=treat[i].strip().split('\t')
	iLine=Input[i].strip().split('\t')
	feature=tLine[:-1]
	valueTreat=scaleFactor[0]*float(tLine[-1])+1
	valueInput=scaleFactor[1]*float(iLine[-1])+1
	if log(valueTreat/valueInput,2)>0:
		feature.append(str(log(valueTreat/valueInput,2)))
	else:
		feature.append('0')
	print >> f, '\t'.join(feature)
	
	
