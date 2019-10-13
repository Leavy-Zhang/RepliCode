#! /usr/bin/python

'''
	This script was used for produce wiggle file using sorted bam file.
    Prerequisite: pysam, numpy
	Author: Leavy Zhang
	Date: 2016, Sep, 16th.
'''
######## imported module ########
import re,sys,os
from collections import Counter
from optparse import OptionParser
import logging
import pysam
import numpy as np
from bisect import insort


################################


use_msg='python script.py [options] arg1 arg2'
descr='This script was used for produce wiggle file using the input bam file.'

parser = OptionParser(usage=use_msg,description=descr,add_help_option=True)

parser.add_option("-i","--input-bam", type="string", dest="inbam",help="The inputted bam file.")
parser.add_option("-w","--sliding-window", dest="slideWindow", type="int",help="The width (bp) of sliding window during tag counting (DEFAULT: 10000).",default="10000")
parser.add_option("-t","--interval-step", dest="intervalStep",type="int",help="Step length (bp) for each counting (DEFAULT: 1000).",default="1000")
parser.add_option("-f","--output-format", dest="outfmt",type="string",help="bg:bedgraph, wig: wiggle.",default="bg")
parser.add_option("-o","--output-wig", dest="outfile",type="string",help="Specify the file in which the result should be stored.")
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



######## predefined class and funcs ########
class extractInfo:
	'''
	This class works for extracting and dealing with information from given pysam object.
	'''
	def __init__(self,bamObj=None,stepLen=None,windowLen=None):
		self.bamObj=bamObj
		self.stepLen=stepLen
		self.windowLen=windowLen

	def fetchChromSize(self):
		headerInfo=pysam.AlignmentFile(self.bamObj,'rb').header
		chromInfo={}
		for contig in headerInfo['SQ']:
			chromInfo[contig['SN']]=contig['LN']
		return chromInfo

	def slidingCount(self):
		posTag={}
		num=0
		for align in pysam.AlignmentFile(self.bamObj,'rb'):
			num+=1
			if num!=0 and num%1000000==0:
				info("%d tags counted for sliding window" % num)
			tagpos=(align.reference_start+1+align.reference_start+1+align.query_alignment_length)/2
			startBinSuffix= 1 if tagpos<self.windowLen/2 else (tagpos-self.windowLen/2)/self.stepLen+1
			endBinSuffix= (tagpos+self.windowLen/2)/self.stepLen+1
			tagbinpos=np.arange(startBinSuffix,endBinSuffix)*self.stepLen
			for binpos in tagbinpos:
				if  align.reference_id != -1:
					try:
						posTag[align.reference_name,binpos]+=1
					except:
						posTag[align.reference_name,binpos]=1
		return posTag


def prepWrite(postag):
	chromPos={}
	for element in postag.keys():
		try:
			insort(chromPos[element[0]],element[1])
		except:
			chromPos[element[0]]=[element[1]]
	return chromPos

def writeFile(chrompos,postag,chrominfo,windowlength,steplength,outfilename,outformat):
	outfile=open(outfilename,'w')
	if outformat =='bg':
		for chrom in sorted(chrompos.keys()):
			info('write profile for chromosome %s' % chrom)
			for pos in chrompos[chrom]:
				if (pos-steplength/2) <= chrominfo[chrom]:
					if pos+steplength/2 <= chrominfo[chrom]:
						print >> outfile, '%s\t%d\t%d\t%d' % (chrom,pos-steplength/2,pos+steplength/2,postag[chrom,pos])
					else:
						print >> outfile, '%s\t%d\t%d\t%d' % (chrom,pos-steplength/2, chrominfo[chrom],postag[chrom,pos])
				else:
					break
	else:
		print >> outfile , 'track type=wiggle_0 name="%s" description="Extended tag pileup using %d bp sliding window by %d bp interval."' % (outfile,windowlength,steplength)
		for chrom in sorted(chrompos.keys()):
			info('write profile for chromosome %s' % chrom)
			print >> outfile, 'variableStep chrom=%s span=%d' % (chrom,windowlength)
			for pos in chrompos[chrom]:
				if pos <= chrominfo[chrom]:
					print >> outfile, '%d\t%d' % (pos,postag[chrom,pos])
	outfile.close()

###############################


######## main function ########
def main():
	inbam=options.inbam
	stepLength=options.intervalStep
	windowLength=options.slideWindow
	outFmt=options.outfmt
	outFileName=options.outfile
	info('Collecting input parameters...')
	info('Input bam file: %s' % inbam)
	info('interval step length during tag counting: %d' % stepLength)
	info('interval sliding window length during tag counting: %d' % windowLength)
	info('format during writing result as output: %s' % outFmt)
	info('Name of output file: %s' % outFileName)
	info('Retrieving genome information')
	allAlignments=extractInfo(inbam,stepLength,windowLength)
	chromInfo=allAlignments.fetchChromSize()
	info('tag counting for sliding widows starts...')
	posTag= allAlignments.slidingCount()
	info('tag counting finished...')
	info('prepare for output...')
	chromPos=prepWrite(posTag)
	writeFile(chromPos,posTag,chromInfo,windowLength,stepLength,outFileName,outFmt)
	info('Job finished! Bye! :-)')

##############################

if __name__=='__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! Bye!\n")
        sys.exit(0)
