# MapSplice3 

__MapSplice3__ is an spliced aligner. It maps RNA-seq onto reference genome, and can be used to detect splice junctions, gene fusions, and circular RNAs.

## System Requirements
MapSplice3 has been tested on Linux platforms with the following system settings.
  * Red Hat 4.4.5-6
  * gcc version 4.9.1

# Manual

There are two subfolders in the package: "indexing" and "mapping".

Code in subfolder "indexing" is used to build index.
Code in subfolder "mapping" is used to perform the spliced alignment.

How to run MapSplice3?

## 1. Build index for MapSplice3.

Note: before building index, you need to put all the sequence files of reference genome into a directory (like /PATH/hg19/). And all the sequence files are required to be in the following format:

(1) In "FASTA" format, with ".fa" extension.
(2) One chromosome per sequence file.
(3) Chromosome name in the header line (">" not included) is the same as the sequence file base name, and does not contain any blank space. 
	E.g. If the header line is ">chr1", then the sequence file name should be "chr1.fa".
(5) No other files in the same folder.


## Installation:
Go to directory "indexing", and run "make buildIndex".

## Running: 

Command Line:
     
```
./buildIndex <path-to-chromosomes-folder> <path-to-index-folder>
```

Example:
 
```
./buildIndex /PATH/hg19/ /PATH/index/
```

## 2. Map reads onto reference genome and detect co-linear splice junctions, back-splice junctions (circular RNA), and gene fusions.

Installation:

Go to directory "mapping", and run "make mps".

Running: There are four modes in total. All the parameters are required for each mode.

### i). map co-linear transcript reads only to detect co-linear splice junctions
Command Line:

```
bin/mps_regular -G <path-to-index> -1 <read-end1> -2 <read-end2> -T <threads-num> -O <output-folder>
```
    
Example:
    
```
bin/mps_regular -G /PATH/index/ -1 /PATH/read_end1.fa -2 /PATH/read_end2.fa -T 16 -O /PATH/mps3_results/
```

### ii). map co-linear transcript and circualr RNA reads to detect both co-linear splice junctions and back-splice junctions
Command Line:

```
bin/mps_regular_circRNA -G <path-to-index> -1 <read-end1> -2 <read-end2> -T <threads-num> -O <output-folder>
```
    
Example:
    
```
bin/mps_regular_circRNA -G /PATH/index/ -1 /PATH/read_end1.fa -2 /PATH/read_end2.fa -T 16 -O /PATH/mps3_results/
```

### iii). map co-linear transcript and gene fusion reads to detect co-linear splice junctions and gene fusions
Command Line:

```
# Step 1
bin/mps_regular_fusion -G <path-to-index> -1 <read-end1> -2 <read-end2> -T <threads-num> -O <output-folder>
```
```
# Step 2
bin/mps3_fusion_post_new \
<path-to-index> \
<output-folder-step1> \
<path-to-gtf-file> \
<threads-num> \
<fusion-junc-read-support-min> \
<path-to-paralog-gene-file> \
<output-folder>
```
    
Example:
    
```
# Step 1
bin/mps_regular_fusion -G /PATH/index/ -1 /PATH/read_end1.fa -2 /PATH/read_end2.fa -T 16 -O /PATH/mps3_results/
```
```
# Step 2
bin/mps3_fusion_post_new \
/PATH/index/ \
/PATH/step1_results/ \
/PATH/gtf_file \
16 \
5 \
/PATH/paralog_gene_file \
/PATH/fusion_post_output/
```

### iv). map co-linear transcript, circualr RNA and gene fusion reads to detect co-linear splice junctions, back-splice junctions, and gene fusions
Command Line:

```
# Step 1
bin/mps_regular_circRNA_fusion -G <path-to-index> -1 <read-end1> -2 <read-end2> -T <threads-num> -O <output-folder>
```
```
# Step 2
bin/mps3_fusion_post_new \
<path-to-index> \
<output-folder-step1> \
<path-to-gtf-file> \
<threads-num> \
<fusion-junc-read-support-min> \
<path-to-paralog-gene-file> \
<output-folder>
```

Example:
    
```
# Step 1
bin/mps_regular_circRNA_fusion -G /PATH/index/ -1 /PATH/read_end1.fa -2 /PATH/read_end2.fa -T 16 -O /PATH/mps3_results/
```
```
# Step 2
bin/mps3_fusion_post_new \
/PATH/index/ \
/PATH/step1_results/ \
/PATH/gtf_file \
16 \
5 \
/PATH/paralog_gene_file \
/PATH/fusion_post_output/
```
<!---
# How to interpret MapSplice3 results:

Two result files are generated:

"output.sam" records alignments (in SAM format) for all the reads;
"output.junc" records detected junctions in following format:

Note: all the positions are 1-based.
A. Each line records one splice junction;
B. Detailed description of all the columns:
	column 1: the chromosome name involved in the junction,
	column 2: last base of the upstream exon,
	column 3: first base of the downstream exon,
	column 4: number of reads aligned to the junction.
For example:
chr1	12227	12613	10
describes a splice junction whose:
upstream exon ends at position chr1: 12227;
downstream exon starts at position chr1: 12613;
involoved intronic region is chr1: 12226 ~ 12612;
supporting read number is 10. -->

### License
Please refer to LICENSE.txt
