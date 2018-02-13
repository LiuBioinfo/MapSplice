// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 outputScript" << endl;
		exit(1);
	}
	string outputScript = argv[1];
	
	string cmd_cd_hisat2dir = "cd /home/lcph222/Xinan/tools/hisat2-2.0.5/";
	string cmd_chrSeqList_1 = "/home/xli262/chrom_Index/hg19_noRandom/chr1.fa,/home/xli262/chrom_Index/hg19_noRandom/chr2.fa,\\";
	string cmd_chrSeqList_2 = "/home/xli262/chrom_Index/hg19_noRandom/chr3.fa,/home/xli262/chrom_Index/hg19_noRandom/chr4.fa,\\";
	string cmd_chrSeqList_3 = "/home/xli262/chrom_Index/hg19_noRandom/chr5.fa,/home/xli262/chrom_Index/hg19_noRandom/chr6.fa,\\";
	string cmd_chrSeqList_4 = "/home/xli262/chrom_Index/hg19_noRandom/chr7.fa,/home/xli262/chrom_Index/hg19_noRandom/chr8.fa,\\";
	string cmd_chrSeqList_5 = "/home/xli262/chrom_Index/hg19_noRandom/chr9.fa,/home/xli262/chrom_Index/hg19_noRandom/chr10.fa,\\";
	string cmd_chrSeqList_6 = "/home/xli262/chrom_Index/hg19_noRandom/chr11.fa,/home/xli262/chrom_Index/hg19_noRandom/chr12.fa,\\";
	string cmd_chrSeqList_7 = "/home/xli262/chrom_Index/hg19_noRandom/chr13.fa,/home/xli262/chrom_Index/hg19_noRandom/chr14.fa,\\";
	string cmd_chrSeqList_8 = "/home/xli262/chrom_Index/hg19_noRandom/chr15.fa,/home/xli262/chrom_Index/hg19_noRandom/chr16.fa,\\";
	string cmd_chrSeqList_9 = "/home/xli262/chrom_Index/hg19_noRandom/chr17.fa,/home/xli262/chrom_Index/hg19_noRandom/chr18.fa,\\";
	string cmd_chrSeqList_10 = "/home/xli262/chrom_Index/hg19_noRandom/chr19.fa,/home/xli262/chrom_Index/hg19_noRandom/chr20.fa,\\";
	string cmd_chrSeqList_11 = "/home/xli262/chrom_Index/hg19_noRandom/chr21.fa,/home/xli262/chrom_Index/hg19_noRandom/chr22.fa,\\";
	string cmd_chrSeqList_12 = "/home/xli262/chrom_Index/hg19_noRandom/chrM.fa,/home/xli262/chrom_Index/hg19_noRandom/chrX.fa,\\";
	string cmd_chrSeqList_13 = "/home/xli262/chrom_Index/hg19_noRandom/chrY.fa, \\";
	string cmd_outputDir = "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/HISAT2_index/";

	int dataSetNum = 74;
	ofstream script_ofs(outputScript.c_str());
	script_ofs << cmd_cd_hisat2dir << endl;
	for(int tmpDataNO = 1; tmpDataNO <= dataSetNum; tmpDataNO ++)
	{
		script_ofs << "mkdir -p " << cmd_outputDir << tmpDataNO << endl;

		script_ofs << "./hisat2-build -p 16 --snp /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/";
		script_ofs << tmpDataNO << "/twoHap_SNP.hisat2.txt \\" << endl;
		script_ofs << cmd_chrSeqList_1 << endl;
		script_ofs << cmd_chrSeqList_2 << endl;
		script_ofs << cmd_chrSeqList_3 << endl;
		script_ofs << cmd_chrSeqList_4 << endl;
		script_ofs << cmd_chrSeqList_5 << endl;
		script_ofs << cmd_chrSeqList_6 << endl;
		script_ofs << cmd_chrSeqList_7 << endl;
		script_ofs << cmd_chrSeqList_8 << endl;
		script_ofs << cmd_chrSeqList_9 << endl;
		script_ofs << cmd_chrSeqList_10 << endl;
		script_ofs << cmd_chrSeqList_11 << endl;
		script_ofs << cmd_chrSeqList_12 << endl;
		script_ofs << cmd_chrSeqList_13 << endl;
		script_ofs << cmd_outputDir << tmpDataNO << "/KGP_" << tmpDataNO << endl << endl;
	}
	script_ofs.close();
	return 0;
}