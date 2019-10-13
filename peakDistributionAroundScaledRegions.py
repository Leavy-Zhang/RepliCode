#! /usr/bin/python
'''
    This script was used to counts of overlapped features from A around features B.
    Prerequisite: BedTools
    Author: Leavy Zhang
    Date: 2016, Oct, 30th.
'''

######## imported module ########

import re,sys,os
from collections import Counter
from optparse import OptionParser
import logging

#################################


use_msg='python script.py [options] arg1 arg2'
descr='This script was used to counts of overlapped query features around target features in bed format. IDs for query and target features were required (namely the fourth column of bed files).'

parser = OptionParser(usage=use_msg,description=descr,add_help_option=True)

parser.add_option("-q","--query-peaks",type="string", dest="qpeaks",help="File of query peak regions in bed format.")
parser.add_option("-r","--target-regions",type="string", dest="tpeaks",help="File of target regions in bed format.")
parser.add_option("-e","--extend-distance", dest="edist",type="int",help="Extended distance (bp) at the downstream and upstream of all target regions (Default: 10000)",default=10000)
parser.add_option("-s","--scale-length", dest="scalelength",type="int",help="scaled length of all target regions (Default: 1000)",default=10000)
parser.add_option("-p","--to-percentage", dest="toPercentage",action="store_true", help="Whether to compute percentage (Default: False)", default=False)
parser.add_option('-b','--bin-size',dest='binsize',type='int',help="the size of bin in which overlapped feature will be counted. (Default:1000)",default=1000)
parser.add_option("-o","--output", dest="outFile",type="string",help="Specify the file in which the result should be stored.",default="NA")
(options,args)=parser.parse_args()


# -----constants-----
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info
# -------------------


def extendRegionAndRun(qpeaks,tpeaks,edistance):
    etpeaksfile=''.join([random.sample(['z','y','x','w','v','u','t','s','r','q','p','o','n','m','l','k','j','i','h','g','f','e','d','c','b','a'], 5),'.bed'])
    etpeaks=open(etpeaksfile,'w')
    for feature in open(tpeaks):
        feature=feature.strip().split('\t')
        feature[1]=str(int(feature[1])-edistance+1 if feature[1]<= edistance else 1)
        feature[2]=str(int(feature[2])+edistance-1)
        print >> etpeaks,'\t'.join(feature)
    etpeaks.close()
    overlappedResult= os.popen('intersectBed -a %s -b %s -wb' % (etpeaksfile,tpeaks))
    return overlappedResult
        
def distributionStat(overlappedResult,qpeaks,tpeaks,scaleLength,edistance,binSize,toPercentage,output):
    num=0
    for i in open(qpeaks):
        num+=1
    rawRegions={}
    for line in open(tpeaks):
        line=line.strip().split('\t')
        rawRegions[line[3]]=[line[1],line[2]]
    designedPos=range(-1*edistance,0+binSize,binSize)
    designedPos.extend(range(0,scaleLength,binSize))
    designedPos.extend(range(scaleLength,scaleLength+edistance+binSize,binSize))
    designedDist={}
    for pos in designedPos:
        designedDist[pos]=0
    for m in overlappedResult.split('\n'):
        m=m.split('\t')
        matchedPos=sum(int(m[1]),int(m[2]))/2
        if rawRegions[m[3]][1] <= matchedPos:
            tpos=matchedPos/binSize*binSize
            designedDist[tpos]+=1
            try:
                designedDist[tpos+binSize]+=1
            except:
                pass
        elif rawRegions[m[3]][0] >= matchedPos:
            tpos=int(scaleLength*(matchedPos-rawRegions[m[3]][0]+1)/float(rawRegions[m[3]][1]-rawRegions[m[3]][0]))/binSize*binSize
            designedDist[tpos]+=1
            designedDist[tpos+binSize]+=1
        else:
            tpos=matchedPos/binSize*binSize
            designedDist[tpos]+=1
            try:
                designedDist[tpos-binSize]+=1
            except:
                pass
    outputFile=open(output,'w')
    if toPercentage:
        countDist=[str(designedDist[i]/float(num)*100) for i in sorted(designedDist.keys())]
    else:
        countDist=[str(designedDist[i]) for i in sorted(designedDist.keys())]
    print >> outputFile, '\t'.join([str(i) for i in sorted(designedDist.keys())])
    print >> outputFile, '\t'.join(countDist)    
    
        
    


def main():
    qpeakFile=options.qpeaks
    tpeaksFile=options.tpeaks
    extendedDistance=options.edist
    scaleLength=options.scaleLength
    binSize=options.binsize
    outfile=options.outFile
    toPercentage=options.toPercentage
    featureCompare=extendRegionAndRun(qpeakFile,tpeaksFile,extendedDistance)
    distributionStat(featureCompare,qpeakFile,tpeaksFile,scaleLength,binSize,toPercentage,outfile)
