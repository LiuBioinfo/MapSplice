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

void getIdVec(string& dataListFile, vector<string>& idVec)
{
	ifstream id_ifs(dataListFile.c_str());
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		if(tabLoc == string::npos)
			idVec.push_back(tmpStr);
		else
			idVec.push_back(tmpStr.substr(0, tabLoc));
	}
	id_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 dataList" << endl;
		cout << "#2 outputScript" << endl;
		exit(1);
	}
	string dataList = argv[1];
	string outputScript = argv[2];

	vector<string> idVec;
	getIdVec(dataList, idVec);

	ofstream script_ofs(outputScript.c_str());
	script_ofs << "cd /home/lcph222/Xinan/mps/imps_github/mps-master/code/" << endl;
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		script_ofs << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output/repeat_region" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output/oneEndUnmapped" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output/incomplete" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output/completePair" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase1_output/bothEndsUnmapped" << endl;
		script_ofs << "mkdir /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/phase2_output" << endl;

		script_ofs << "./iMapSplice \\" << endl << "-P /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNP.txt \\" << endl
			<< "-Q /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNPmer_index/ \\" << endl
			<< "-G /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416 \\" << endl
			<< "-L /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416/2ndLevelIndex/ \\" << endl
			<< "-1 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/1_paired.fastq \\" << endl
			<< "-2 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/2_paired.fastq \\" << endl
			<< "-T 16 -O /scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" << tmpId << "_impsTwoHapSingleIndex/" << endl;
	}
	script_ofs.close();
	return 0;
}