## Code for sequencing data analysis of DNA replication

### Major codes for paper "H2A.Z Facilitates Licensing and Activation of Early Replication" are deposited here. There might be some error during running them, please let me know if needed (E-mail: lwzhanghz@163.com)

---
### slidingCountFromBam.py
#### Description:
This script was used for produce genomic profiles in wiggle/bedgraph format using the input bam file.
#### Prerequisite:
*Python module:* pysam, numpy
#### Usage:
python script.py [options] arg1 arg2

##### Options:
  -h, --help&emsp;show this help message and exit</br>
  -i INBAM, --input-bam=INBAM</br>
  &emsp;&emsp;&emsp;The inputted bam file.</br>
  -w SLIDEWINDOW, --sliding-window=SLIDEWINDOW</br>
&emsp;&emsp;&emsp;The width (bp) of sliding window during tag counting</br>
&emsp;&emsp;&emsp;(DEFAULT: 10000).</br>
  -t INTERVALSTEP, --interval-step=INTERVALSTEP</br>
&emsp;&emsp;&emsp;Step length (bp) for each counting (DEFAULT: 1000).</br>
  -f OUTFMT, --output-format=OUTFMT</br>
&emsp;&emsp;&emsp;bg:bedgraph, wig: wiggle.</br>
  -o OUTFILE, --output-wig=OUTFILE</br>
&emsp;&emsp;&emsp;Specify the file in which the result should be stored.</br>
##### Example:</br>
&emsp;&emsp;python slidingCountFromBam.py -i input.bam -w 50000 -t 1000 -f bg -o output.bedgraph</br>

---
### computeS50.py</br>
#### Description:</br>
H2A.Z Facilitates Licensing and Activation of Early Replication. S50 algorithm was used to describe the DNA replication-timing of genomic regions of interest using reli-seq data [Gaetano Ivan Dellino and *et. al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530669/)

Input files are a series of read per million (RPM), or other algorithm, normalized bedgraph/wiggle formatted files of time-course repli-seq data, such as S1,S2,S3 and *et.al*.</br>
1. Genomic bins of input files should be generated with the same parameters using **slidingCountFromBam.py** to ensure the consistency of genomic bin positions between files.</br>
2. Using **intersectBed** to retrieve read count for each genomic bin at each time point and then normalized using RPM-like method. An example of files is shown below:</br>

</br>chr&emsp;start&emsp;end&emsp;normValue<br/>
chr1&emsp;1&emsp;5000&emsp;0.36<br/>
chr1&emsp;5001&emsp;10000&emsp;0.37<br/>
chr1&emsp;10001&emsp;15000&emsp;0.34<br/>

3. Run script to calculate S50 score for all genomic bins.</br>

#### Usage:
python script.py [options] arg1 arg2

Options:</br>
  -h, --help&emsp;show this help message and exit</br>
  -i INFILES, --input-files=INFILES</br>
&emsp;&emsp;&emsp;The inputted wig/bedgraph files separated by comma.</br>
  -t FILETYPE, --file-type=FILETYPE</br>
&emsp;&emsp;&emsp;Format of input files and output file, wiggle or</br>
&emsp;&emsp;&emsp;bedgraph.(Default: bedgraph)</br>
  -o OUTFILE, --output-wig=OUTFILE</br>
&emsp;&emsp;&emsp;Specify the file in which the result should be stored.</br></br>

#### Example:
&emsp;&emsp;python computeS50.py -i S1.bedgraph,S2.bedgraph,S3.bedgraph,S4.bedgraph -t bedgraph -o S50.bedgraph

---
### F-score.py
#### Description:
F-score (firing score) was used to assess the effective firing efficiency of genomic regions using nascent strand (NS)-seq data, each sample contains two parallel data, RNase untreated and treated sequencing data.

#### Usage:
1. Region files can be prepared similarly as that of **computeS50.py** or a set of pre-predicted peak regions in bed format.</br>
2. Using **intersectBed** to **intersectBed** to retrieve read count for region/peaks of interest.</br>
3. Run script to compute F-score:</br>
&emsp;&emsp;python F-score.py total_reads_of_untreated,total_reads_of_treated untreated.bed treated.bed output.F-score</br>
   
#### Example:
&emsp;python F-score.py 17321232,12987432 untreated.bed treated.bed my.F-score
   
