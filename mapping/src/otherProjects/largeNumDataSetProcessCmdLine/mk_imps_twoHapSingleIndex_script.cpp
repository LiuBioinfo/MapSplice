#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>
#include <sstream>

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
	ofstream script_ofs(outputScript.c_str());
	int dataNO_1st = 6;
	int dataNO_last = 74;
	script_ofs << "cd /home/lcph222/Xinan/mps/imps_github/mps-master/code/" << endl;
	for(int tmp = dataNO_1st; tmp <= dataNO_last; tmp++)
	{
		script_ofs << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output/repeat_region" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output/oneEndUnmapped" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output/incomplete" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output/completePair" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase1_output/bothEndsUnmapped" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/phase2_output" << endl;

		script_ofs << "./iMapSplice \\" << endl << "-P /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmp << "/twoHap_SNP.txt \\" << endl
			<< "-Q /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmp << "/twoHap_SNPmer_index/ \\" << endl
			<< "-G /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416 \\" << endl
			<< "-L /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416/2ndLevelIndex/ \\" << endl
			<< "-1 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmp << "/1_paired.fastq \\" << endl
			<< "-2 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmp << "/2_paired.fastq \\" << endl
			<< "-T 16 -O /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmp << "_impsTwoHapSingleIndex/" << endl;
	}
	script_ofs.close();
	return 0;
}