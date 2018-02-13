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
		script_ofs << "cat \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" 
			<< tmp << "/unique_1_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" 
			<< tmp << "/unique_2_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" 
			<< tmp << "/shared_SNP.txt > \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" 
			<< tmp << "/twoHap_SNP.txt" << endl << endl;		
		script_ofs << "./generateSyntheticSNPgenomSeq /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416/ \\"
			<< endl << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/genCode.reissuedID.gaf \\" << endl
			<< "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmp << "/twoHap_SNP.txt \\" << endl
			<< "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmp << "/twoHap_SNPmer 201" << endl << endl;
		script_ofs << "./buildWholeGenome_mergedFa_noPreIndex \\" << endl << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" 
			<< tmp << "/twoHap_SNPmer/SNPinAnn.fa \\" << endl << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/"
			<< tmp << "/twoHap_SNPmer_index" << endl;
	}
	script_ofs.close();
	return 0;
}