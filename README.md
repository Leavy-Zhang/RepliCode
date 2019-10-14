## Code for sequencing data analysis of DNA replication

### Major codes for paper "H2A.Z Facilitates Licensing and Activation of Early Replication" are deposited here. There might be some error during running them, please let me know if needed (E-mail: lwzhanghz@163.com)

---
### slidingCountFromBam.py
#### Description:
This script was used for produce wiggle file using the input bam file.
#### Prerequisite:
*Python module:* pysam, numpy
#### Usage:
python script.py [options] arg1 arg2

##### Options:
  -h, --help            show this help message and exit
  -i INBAM, --input-bam=INBAM
                        The inputted bam file.
  -w SLIDEWINDOW, --sliding-window=SLIDEWINDOW
                        The width (bp) of sliding window during tag counting
                        (DEFAULT: 10000).
  -t INTERVALSTEP, --interval-step=INTERVALSTEP
                        Step length (bp) for each counting (DEFAULT: 1000).
  -f OUTFMT, --output-format=OUTFMT
                        bg:bedgraph, wig: wiggle.
  -o OUTFILE, --output-wig=OUTFILE
                        Specify the file in which the result should be stored.
##### Example:
python slidingCountFromBam.py -i input.bam -w 50000 -t 1000 -f bg -o output.bedgraph

---
### computeS50.py
#### Description:
H2A.Z Facilitates Licensing and Activation of Early Replication. S50 algorithm was used to describe the DNA replication-timing of genomic regions of interest using reli-seq data [Gaetano Ivan Dellino and *et. al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530669/)

Input files are a series of read per million (RPM), or other algorithm, normalized bedgraph/wiggle formatted files of time-course repli-seq data, such as S1,S2,S3 and *et.al*.
1. Genomic bins of input files should be generated with the same parameters using **slidingCountFromBam.py** to ensure the consistency of genomic bin positions between files.
2. Using **intersectBed** to retrieve read count for each genomic bin at each time point and then normalized using RPM-like method. An example of files is shown below:

</br>chr&emsp;start&emsp;end&emsp;normValue<br/>
</br>chr1&emsp;1&emsp;5000&emsp;0.36<br/>
</br>chr1&emsp;5001&emsp;10000&emsp;0.37<br/>
</br>chr1&emsp;10001&emsp;15000&emsp;0.34<br/>

3. Run script to calculate S50 score for all genomic bins.

#### Usage:
python script.py [options] arg1 arg2

Options:
  -h, --help            show this help message and exit
  -i INFILES, --input-files=INFILES
                        The inputted wig/bedgraph files separated by comma.
  -t FILETYPE, --file-type=FILETYPE
                        Format of input files and output file, wiggle or
                        bedgraph.(Default: bedgraph)
  -o OUTFILE, --output-wig=OUTFILE
                        Specify the file in which the result should be stored.

#### Example:
python computeS50.py -i S1.bedgraph,S2.bedgraph,S3.bedgraph,S4.bedgraph -t bedgraph -o S50.bedgraph

---
### F-score.py
#### Description:
F-score (firing score) was used to assess the effective firing efficiency of genomic regions using nascent strand (NS)-seq data, each sample contains two parallel data, RNase untreated and treated sequencing data.

#### Usage:
1. Region files can be prepared similarly as that of **computeS50.py** or a set of pre-predicted peak regions in bed format.
2. Using **intersectBed** to **intersectBed** to retrieve read count for region/peaks of interest.
3. Run script to compute F-score:
   python F-score.py total_reads_of_untreated,total_reads_of_treated untreated.bed treated.bed output.F-score
   
#### Example:
   python F-score.py 17321232,12987432 untreated.bed treated.bed my.F-score
   
