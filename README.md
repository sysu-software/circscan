# circscan
circScan: Discovering the Circular RNA Interactions with RNA-binding Proteins and Ribosomes

/*****************************************************************************************
 *	circScan - An algorithm for discovering circRNA interactions with RNA-binding proteins 
 *  (RBPs) and ribosomes from short-read (18-50nt) CLIP-seq and ribo-seq datasets
 *
 *	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>
 * 
 *	RNA Information Center, School of Life Sciences, Sun Yat-Sen University
 *	
 *  Create date: 09/18/2016
 *  
 ****************************************************************************************/

Overview:
---------
circScan is an algorithm for identifying circRNA interactions with RBPs and ribosomes from 
short-read(~18-50nt) CLIP-seq and ribosome profiling (ribo-seq) data.
Using circScan, one can easy and efficient to detect circRNAs from various CLIP-seq and ribo-seq data. 

Usage:
---------
Usage:  circScan [options] <genome file> <genome fai> <annotation file, bed12> <mapped alignments, bam> <unmapped read, fa><BR>
File format for annotation must be UCSC bed12<BR>
Mapped alignments must be sorted BAM format<BR>

[options]<BR>
-v/--verbose                : verbose information<BR>
-V/--version                : circScan version
-o/--outfile <string>       : output file
-p/--parclip                : par-clip datasets
-f/--fdr                    : false discovery rate(<0.05)
-m/--max-num <int>          : maximum locations for each read [default=1]
-s/--seed-len <int>         : the length of seed for searching [default>=10]
-x/--mismatch <int>         : maximum mismatch number for each read [default<=1]
-e/--error-rate <double>    : maximum error rate [default<=0.05], =mismatch/readLen
-i/--min-score <int>        : minimum score [default>=20], =matches - mismatches

Installation:
---------
Download circScan-0.1.tar.gz from http://rna.sysu.edu.cn/circscan/Download.php; unpack it, and make:
tar -xzvf circScan-0.1.tar.gz 
cd circScan-0.1
make
The newly compiled binary (circScan) is in the circScan /bin directory. You can run it from there, as in this example:
bin/circScan ./test_data/testGenome.fa ./test_data/testGenome.fa.fai ./test_data/testAnnoFile.bed12 ./test_data/testGenome.sorted.bam ./test_data/testGenome.unmapped.fa >./test_data/test_circScan_candidates.txt
You can check the output file test_circScan_candidates.txt in test_data directory.

System requirements:
---------
Operating system: circScan is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.
Libraries and other installation requirements: circScan includes one software library: the BamTools library package. All will automatically compile during circScan installation process. 
By default, circScan does not require any additional libraries to be installed by you.

Prerequisites
---------
Dependencies: The input of circScan is sorted BAM file. So you need the short read mapper bowtie2 or other mappers (e.g. STAR) and samtools up and running.
You can get the most fresh versions: 
(1)	Samtools: http://www.htslib.org/download/;
(2)	Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

At this point you should have everything to run a built-in test data set
cd test_data
./run_test.sh
If you get error messages here, first make sure the Prerequisites are
really installed correctly and run on their own.

If everything goes well you can get started with your real data! :)   
You need to have the reference genome, annotation file, and  bowtie2 indexes for genome and annotation.
Option 1: You can download these datasets in our circScan download web page: http://rna.sysu.edu.cn/circscan/Download.php
Option 2: You can constructed these datasets by yourself using following steps:
As an example, let's assume you use human genome (version hg38).
(1)	Genome:
mkdir genome
wget -c 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
gzip -d hg38.fa.gz
cd ..
(2)	Annotation:
You can use Table Browser to get the bed12 and fasta files for genome annotation(e.g. GENCODE)
http://genome.ucsc.edu/cgi-bin/hgTables
e.g. You can save the output files in genome directory: hg38.gencode.bed12 and hg38.gencode.fa
(3) Build the genome index:
bowtie2-build genome/hg38.fa genome/hg38
(4)build the fai index:
samtools faidx genome/hg38.fa
(5)build the transcriptome index:
bowtie2-build genome/hg38.gencode.fa genome/hg38.gencode

How to get unmapped reads:
---------
1. aligned the short CLIP-seq reads to genome and transcriptome
(1) bowtie2 -p 16 -t -k 4 --no-unal -D 200 -R 3 -N 0 -L 15 -i S,1,0.5 \
    --score-min=C,-16,0 -f -U <your_reads.fa> -x ./genome/hg38 --un ./genome/genome.unmapped.fa \
    >./genome/genome.alignments.bwt2 2>./genome/bowtie2.genome.log
(2) cd genome;
(3) samtools view -bS genome.alignments.bwt2 >genome.alignments.bam
(4) samtools sort genome.alignments.bam genome.alignments.sorted
(5) samtools index genome.alignments.sorted.bam
(6) aligned the unmapped reads to transcriptome
    bowtie2 -p 16 -t --no-unal -D 200 -R 3 -N 0 -L 15 -i S,1,0.5 \
    --score-min=C,-16,0 -f -U genome.unmapped.fa -x ./genome/hg38.gencode --un ./genome/genome.transcriptome.unmapped.fa \
    >./genome/transcriptome.alignments.bwt2 2>./genome/bowtie2.transcriptome.log

Now, you can get the unmapped read file "genome.transcriptome.unmapped.fa" in genome directory for the circScan

run circScan:
---------
bin/circScan ./genome/hg38.fa ./genome/hg38.fa.fai ./genome/hg38.gencode.bed12 ./genome/genome.alignments.sorted.bam ./genome/genome.transcriptome.unmapped.fa \
>./genome/circScan_candidate_circRNAs.txt

Output:
---------
#chrom	chromStart	chromEnd	name	score	strand	mismatchNum	t2cNum	acceptReadNum	donorReadNum	FDR	junctionPosInRead	readName	readSeq
testSeq	50	1535	ENST00000625883.1|CDR1-AS|ENSG00000281508.1|CDR1-AS-001|antisense	52	+	0	0	3	5	0.00000	22	1665211-1	TTCCAGTGTCTGCAATATCCAGGGTTTCCGATGGCACCTGTGTCAAGGTCTT
Note: # is comment line

Acknowledgements:
---------
Thanks a lot to everyone who contributed to the public code (e.g. BamTools) used by circScan.

Contact :
---------
JianHua, Yang<yangjh7@mail.sysu.edu.cn or lsp03yjh@gmail.com>
last modified time: /01/2017/
