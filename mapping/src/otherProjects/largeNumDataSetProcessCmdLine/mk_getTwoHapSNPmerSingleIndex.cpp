// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
// cd /home/lcph222/Xinan/MPS3/mapping/ ./cmpTwoHaplotypeSNP
// cat > twoHap_SNP.txt
// cd /home/lcph222/Xinan/mps/imps_github/mps-master/code ./generateSyntheticSNPgenomSeq
// ./buildWholeGenome_mergedFa_noPreIndex
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
		cout << "#1 outputScript" << endl;
		exit(1);
	}
	string dataList = argv[1];
	string outputScript = argv[2];
	vector<string> idVec;
	getIdVec(dataList, idVec);

	ofstream script_ofs(outputScript.c_str());
	script_ofs << "cd /home/lcph222/Xinan/MPS3/mapping/" << endl;
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		script_ofs << endl;
		script_ofs << "./cmpTwoHaplotypeSNP /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416 \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/formattedSNP_file/" << tmpId << "_hap1_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/formattedSNP_file/" << tmpId << "_hap2_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << endl << endl;
		script_ofs << "cat /scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/shared_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/unique_1_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/unique_2_SNP.txt > \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNP.txt" << endl;
	}

	script_ofs << endl << "cd /home/lcph222/Xinan/mps/imps_github/mps-master/code" << endl;
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		script_ofs << endl;
		script_ofs << "./generateSyntheticSNPgenomSeq /home/xli262/chrom_Index/adSA_table/hg19_noRandom_table/index_0416/ \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/genCode.reissuedID.gaf \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNP.txt \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNPmer 201" << endl << endl;
		script_ofs << "./buildWholeGenome_mergedFa_noPreIndex \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNPmer/SNPinAnn.fa \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/SNP_twoHaplotype/" << tmpId << "/twoHap_SNPmer_index" << endl;
	}

	script_ofs.close();
	return 0;
}