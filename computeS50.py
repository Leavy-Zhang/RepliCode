#! /usr/bin/python

'''
    This script was used for computing S50 with repli-seq.
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
descr='This script was used for computing S50 score using repli-seq data.'

parser = OptionParser(usage=use_msg,description=descr,add_help_option=True)

parser.add_option("-i","--input-files",type="string", dest="infiles",help="The inputted wig/bedgraph files separated by comma.")
parser.add_option("-t","--file-type", dest="filetype",type="string",help="Format of input files and output file, wiggle or bedgraph.(Default: bedgraph)",default="bedgraph")
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



class dataIntegrate:
    '''
    This class works for extracting and dealing with information from the input file list.
    '''
    def __init__(self,dataBase=None,newFile=None,listLen=None,fileTurn=None):
        self.fileturn=fileTurn
        self.db=dataBase
        self.newfile=newFile
        self.listlen=listLen

    def wiggleParse(self):
        for i in open(self.newfile):
            if re.match('track',i):
                pass
            elif re.match('variable',i):
                chrom=i.split(' ')[1].split('=')[1]
            else:
                i=i.strip().split('\t')
                try:
                    self.db[(chrom,int(i[0]))][self.fileturn]=float(i[-1])
                except:
                    self.db[(chrom,int(i[0]))]=[0.0]*self.listlen
                    self.db[(chrom,int(i[0]))][self.fileturn]=float(i[-1])
        return self.db

    def bedGraphParse(self):
        for i in open(self.newfile):
            if re.match('track',i):
                pass
            else:
                i=i.strip().split('\t')
                try:
                    self.db[(chrom,int(i[0]),int(i[1]))][self.fileturn]=float(i[-1])
                except:
                    self.db[(chrom,int(i[0]),int(i[1]))]=[0.0]*self.listlen
                    self.db[(chrom,int(i[0]),int(i[1]))][self.fileturn]=float(i[-1])
        return self.db

class S50compute:
    '''
    This class works for extracting and dealing with information from the input file list.
    '''
    def __init__(self, dataBase=None,listLen=None,outFile=None):
        self.db=dataBase
        self.listlen=listLen
        self.outfile=outFile
    def wiggleS50(self):
        self.outfile.writelines('track type=wiggle_0\n')
        ChrSite={}
        for key in self.db.keys():
            try:
                ChrSite[key[0]].append(key[1])
            except:
                ChrSite[key[0]]=[key[1]]
        for chrom in ChrSite.keys():
            pos=list(set(ChrSite[chrom]))
            print >> self.outfile, 'variableStep chrom=%s span=%d' % (chrom,pos[1]-pos[0])
            for singlePos in pos:
                cumSum=[sum(self.db[chrom,singlePos])[:k+1] for k in range(self.listlen)]
                dist={}
                for k in cumSum:
                    dist[abs(k-sum(self.db[chrom,singlePos])/2)] = k
                tkey=sorted(dist.keys())
                if len(tkey)>=2:
                    print >> self.outfile, '%d\t%f' % (singlePos,sum([dist[tkey[0]],dist[tkey[1]]])/2/sum(self.db[chrom,singlePos]))
                else:
                    print >> self.outfile, '%d\t%f' % (singlePos,dist[tkey[0]]/sum(self.db[chrom,singlePos]))
                break
        
    def bedGraphS50(self):
        self.outfile.writelines('track type=bedGraph"\n')
        ChrSite={}
        for key in self.db.keys():
            try:
                ChrSite[key[0]].append([key[1],key[2]])
            except:
                ChrSite[key[0]]=[[key[1],key[2]]]
        for chrom in ChrSite.keys():
            sePos={}
            for se in ChrSite[chrom]:
                sePos[se[0]]=se[1]
            posKey=list(set(sePos.keys()))
            for key in posKey:
                cumSum=[sum(self.db[chrom,key,sePos[key]])[:k+1] for k in range(self.listlen)]
                dist={}
                for k in cumSum:
                    dist[abs(k-sum(self.db[chrom,key,sePos[key]])/2)] = k
                tkey=sorted(dist.keys())
                if len(tkey)>=2:
                    print >> self.outfile, '%s\t%d\t%d\t%f' % (chrom,key,sePos[key],sum([dist[tkey[0]],dist[tkey[1]]])/2/sum(chroms[chrom,key,sePos[key]]))
                else:
                    print >> self.outfile, '%s\t%d\t%d\t%f' % (chrom,key,sePos[key],dist[tkey[0]]/sum(chroms[chrom,key,sePos[key]]))
                break
    
def main():
    infiles=options.infiles
    filetype=options.filetype
    outfile=options.outfile
    info('Collecting input parameters...')
    info('Input file list: %s' % infiles)
    info('Format during writing result as output: %s' % filetype)
    info('Name of output file: %s' % outfile)
    for turn in xrange(len(infiles.split(','))):
        try:
            rawdb
        except:
            rawdb={}
        singleFile=infiles[turn]
        info('Read information from %s' % singleFile)
        readinfo=dataIntegrate(rawdb,singleFile,len(infiles.split(',')),turn)
        if filetype=='bedgraph':
            rawdb=readinfo.bedGraphParse()            
        elif filetype=='wiggle':
            rawdb=readinfo.wiggleParse()
        else:
            sys.exit('unknown file format, sorry! Bye~')
            
    info('Write to output file.')
    if filetype == 'bedgraph':
        fetchS50=S50compute(rawdb,len(infiles.split(',')),outfile)
        fetchS50.bedGraphS50()
    elif filetype=='wiggle':
        fetchS50=S50compute(rawdb,len(infiles.split(',')),outfile)
        fetchS50.wiggleS50()
    else:
        sys.exit('unknown file format, sorry! Bye~')
    info('Job finished!')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit('User interrupt me! Bye!\n')
