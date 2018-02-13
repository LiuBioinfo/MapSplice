MapSplice3 Manual

There are two subfolders in the package: "indexing" and "mapping".

Code in subfolder "indexing" is used to build index.
Code in subfolder "mapping" is used to perform the spliced alignment.


/**********************************************************************/
/**********************************************************************/
The current version MapSplice3 is tested in the following environment:
OS: Linux Red Hat 4.4.6-3
Compiler: g++ 4.4.6


/**********************************************************************/
/**********************************************************************/
How to run MapSplice3?

1. Build index for MapSplice3.

Note: before building index, you need to put all the sequence files of reference genome into a directory (like /PATH/hg19/). And all the sequence files are required to be in the following format:

(1) In "FASTA" format, with ".fa" extension.
(2) One chromosome per sequence file.
(3) Chromosome name in the header line (">" not included) is the same as the sequence file base name, and does not contain any blank space. 
	E.g. If the header line is ">chr1", then the sequence file name should be "chr1.fa".
(5) No other files in the same folder.


Installation:
     Go to directory "indexing", and run "make buildIndex".

Running: 

     Command Line: 
         ./buildIndex <path-to-chromosomes-folder> <path-to-index-folder>

     Example: 
         ./buildIndex /PATH/hg19/ /PATH/index/


2. Map reads onto reference genome.

Installation:
     Go to directory "mapping", and run "make mps".

Running: all the parameters are required.

     Command Line: 
         ./mps -J -G <path-to-index> -1 <read-end1> -2 <read-end2> -T <threads-num> -O <output-folder>

     Example: 
         ./mps -J -G /PATH/index/ -1 /PATH/read_end1.fa -2 /PATH/read_end2.fa -T 16 -O /PATH/mps3_results/


/**********************************************************************/
/**********************************************************************/
How to interpret MapSplice3 results:

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
supporting read number is 10. 