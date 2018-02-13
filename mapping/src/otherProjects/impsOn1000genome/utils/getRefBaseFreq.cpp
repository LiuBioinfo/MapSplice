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
#include <hash_map>
#include <map>
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 alleleReadCountFile" << endl;
		cout << "#2 outputResults" << endl;
		cout << "#3 minReadCount" << endl;
		exit(1);
	}
	string alleleReadCountFile = argv[1];
	string outputResults = argv[2];
	string minReadCountStr = argv[3];
	int minReadCount = atoi(minReadCountStr.c_str());

	// int totalCount_mps = 0;
	// int totalCount_imps = 0;
	// int totalCount_hisat2genome = 0;
	// int totalCount_hisat2iSNP = 0;
	// int totalCount_hisat2pSNP = 0;
	// int totalCount_star = 0;
	// vector<int> refReadCountVec_mps;
	// vector<int> refReadCountVec_imps;
	// vector<int> refReadCountVec_hisat2genome;
	// vector<int> refReadCountVec_hisat2iSNP;
	// vector<int> refReadCountVec_hisat2pSNP;
	// vector<int> refReadCountVec_star;
	// for(int tmp = 0; tmp <= (100/binSize); tmp++)
	// {
	// 	refReadCountVec_mps.push_back(0);
	// 	refReadCountVec_imps.push_back(0);
	// 	refReadCountVec_hisat2genome.push_back(0);
	// 	refReadCountVec_hisat2iSNP.push_back(0);
	// 	refReadCountVec_hisat2pSNP.push_back(0);
	// 	refReadCountVec_star.push_back(0);
	// }
	string outputResults_mps = outputResults + ".mps";
	string outputResults_imps = outputResults + ".imps";
	string outputResults_mps_mask = outputResults + ".mps_mask";
	string outputResults_hisat2genome = outputResults + ".hisat2genome";
	string outputResults_hisat2iSNP = outputResults + ".hisat2iSNP";
	string outputResults_hisat2pSNP = outputResults + ".hisat2pSNP";
	string outputResults_hisat2mask = outputResults + ".hisat2mask";
	string outputResults_star = outputResults + ".star";
	string outputResults_star_mask = outputResults + ".star_mask";
	ifstream alleleReadCount_ifs(alleleReadCountFile.c_str());
	ofstream results_ofs_mps(outputResults_mps.c_str());
	ofstream results_ofs_imps(outputResults_imps.c_str());
	ofstream results_ofs_mps_mask(outputResults_mps_mask.c_str());
	ofstream results_ofs_hisat2genome(outputResults_hisat2genome.c_str());
	ofstream results_ofs_hisat2iSNP(outputResults_hisat2iSNP.c_str());
	ofstream results_ofs_hisat2pSNP(outputResults_hisat2pSNP.c_str());
	ofstream results_ofs_hisat2mask(outputResults_hisat2mask.c_str());
	ofstream results_ofs_star(outputResults_star.c_str());
	ofstream results_ofs_star_mask(outputResults_star_mask.c_str());
	while(!alleleReadCount_ifs.eof())
	{
		string tmpStr;
		getline(alleleReadCount_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
		int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
		int tabLoc_8 = tmpStr.find("\t", tabLoc_7 + 1);
		int tabLoc_9 = tmpStr.find("\t", tabLoc_8 + 1);
		int tabLoc_10 = tmpStr.find("\t", tabLoc_9 + 1);
		int tabLoc_11 = tmpStr.find("\t", tabLoc_10 + 1);
		int tabLoc_12 = tmpStr.find("\t", tabLoc_11 + 1);
		int tabLoc_13 = tmpStr.find("\t", tabLoc_12 + 1);
		int tabLoc_14 = tmpStr.find("\t", tabLoc_13 + 1);
		int tabLoc_15 = tmpStr.find("\t", tabLoc_14 + 1);
		int tabLoc_16 = tmpStr.find("\t", tabLoc_15 + 1);
		int tabLoc_17 = tmpStr.find("\t", tabLoc_16 + 1);						
		string refCount_mps_str = tmpStr.substr(0, tabLoc_1);
		string altCount_mps_str = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string refCount_imps_str = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string altCount_imps_str = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string refCount_mps_mask_str = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string altCount_mps_mask_str = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);		
		string refCount_hisat2genome_str = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
		string altCount_hisat2genome_str = tmpStr.substr(tabLoc_7 + 1, tabLoc_8 - tabLoc_7 - 1);
		string refCount_hisat2iSNP_str = tmpStr.substr(tabLoc_8 + 1, tabLoc_9 - tabLoc_8 - 1);
		string altCount_hisat2iSNP_str = tmpStr.substr(tabLoc_9 + 1, tabLoc_10 - tabLoc_9 - 1);
		string refCount_hisat2pSNP_str = tmpStr.substr(tabLoc_10 + 1, tabLoc_11 - tabLoc_10 - 1);
		string altCount_hisat2pSNP_str = tmpStr.substr(tabLoc_11 + 1, tabLoc_12 - tabLoc_11 - 1);
		string refCount_hisat2mask_str = tmpStr.substr(tabLoc_12 + 1, tabLoc_13 - tabLoc_12 - 1);
		string altCount_hisat2mask_str = tmpStr.substr(tabLoc_13 + 1, tabLoc_14 - tabLoc_13 - 1);		
		string refCount_star_str = tmpStr.substr(tabLoc_14 + 1, tabLoc_15 - tabLoc_14 - 1);
		string altCount_star_str = tmpStr.substr(tabLoc_15 + 1, tabLoc_16 - tabLoc_15 - 1);
		string refCount_star_mask_str = tmpStr.substr(tabLoc_16 + 1, tabLoc_17 - tabLoc_16 - 1);
		string altCount_star_mask_str = tmpStr.substr(tabLoc_17  + 1);

		int refCount_mps = atoi(refCount_mps_str.c_str());
		int altCount_mps = atoi(altCount_mps_str.c_str());
		
		int refCount_imps = atoi(refCount_imps_str.c_str());
		int altCount_imps = atoi(altCount_imps_str.c_str());

		int refCount_mps_mask = atoi(refCount_mps_mask_str.c_str());
		int altCount_mps_mask = atoi(altCount_mps_mask_str.c_str());
			
		int refCount_hisat2genome = atoi(refCount_hisat2genome_str.c_str());
		int altCount_hisat2genome = atoi(altCount_hisat2genome_str.c_str());
		
		int refCount_hisat2iSNP = atoi(refCount_hisat2iSNP_str.c_str());
		int altCount_hisat2iSNP = atoi(altCount_hisat2iSNP_str.c_str());
		
		int refCount_hisat2pSNP = atoi(refCount_hisat2pSNP_str.c_str());
		int altCount_hisat2pSNP = atoi(altCount_hisat2pSNP_str.c_str());

		int refCount_hisat2mask = atoi(refCount_hisat2mask_str.c_str());
		int altCount_hisat2mask = atoi(altCount_hisat2mask_str.c_str());	

		int refCount_star = atoi(refCount_star_str.c_str());
		int altCount_star = atoi(altCount_star_str.c_str());

		int refCount_star_mask = atoi(refCount_star_mask_str.c_str());
		int altCount_star_mask = atoi(altCount_star_mask_str.c_str());				

		int count_mps = refCount_mps + altCount_mps;
		int count_imps = refCount_imps + altCount_imps;	
		int count_mps_mask = refCount_mps_mask + altCount_mps_mask;
		int count_hisat2genome = refCount_hisat2genome + altCount_hisat2genome;
		int count_hisat2iSNP = refCount_hisat2iSNP + altCount_hisat2iSNP;
		int count_hisat2pSNP = refCount_hisat2pSNP + altCount_hisat2pSNP;
		int count_hisat2mask = refCount_hisat2mask + altCount_hisat2mask;
		int count_star = refCount_star + altCount_star;		
		int count_star_mask = refCount_star_mask + altCount_star_mask;		

		int refFreq_mps, refFreq_imps, refFreq_mps_mask, refFreq_hisat2genome, 
			refFreq_hisat2iSNP, refFreq_hisat2pSNP, refFreq_hisat2mask, refFreq_star, refFreq_star_mask;
		if(count_mps >= minReadCount)
		{
			refFreq_mps = 100 * (double)refCount_mps / (double)count_mps;
			results_ofs_mps << refFreq_mps << endl;
		}
		if(count_imps >= minReadCount)
		{
			refFreq_imps = 100 * (double)refCount_imps / (double)count_imps;
			results_ofs_imps << refFreq_imps << endl;
		}
		if(count_mps_mask >= minReadCount)
		{
			refFreq_mps_mask = 100 * (double)refCount_mps_mask / (double)count_mps_mask;
			results_ofs_mps_mask << refFreq_mps_mask << endl;
		}

		if(count_hisat2genome >= minReadCount)
		{
			refFreq_hisat2genome = 100 * (double)refCount_hisat2genome / (double)count_hisat2genome;
			results_ofs_hisat2genome << refFreq_hisat2genome << endl;
		}
		if(count_hisat2iSNP >= minReadCount)
		{
			refFreq_hisat2iSNP = 100 * (double)refCount_hisat2iSNP / (double)count_hisat2iSNP;
			results_ofs_hisat2iSNP << refFreq_hisat2iSNP << endl;
		}
		if(count_hisat2pSNP >= minReadCount)
		{
			refFreq_hisat2pSNP = 100 * (double)refCount_hisat2pSNP / (double)count_hisat2pSNP;
			results_ofs_hisat2pSNP << refFreq_hisat2pSNP << endl;
		}
		if(count_hisat2mask >= minReadCount)
		{
			refFreq_hisat2mask = 100 * (double)refCount_hisat2mask / (double)count_hisat2mask;
			results_ofs_hisat2mask << refFreq_hisat2mask << endl;
		}


		if(count_star >= minReadCount)
		{
			refFreq_star = 100 * (double)refCount_star / (double)count_star;
			results_ofs_star << refFreq_star << endl;
		}
		if(count_star_mask >= minReadCount)
		{
			refFreq_star_mask = 100 * (double)refCount_star_mask / (double)count_star_mask;
			results_ofs_star_mask << refFreq_star_mask << endl;
		}		
	}
	results_ofs_mps.close();
	results_ofs_imps.close();
	results_ofs_mps_mask.close();
	results_ofs_hisat2genome.close();
	results_ofs_hisat2iSNP.close();
	results_ofs_hisat2pSNP.close();
	results_ofs_hisat2mask.close();
	results_ofs_star.close();
	results_ofs_star_mask.close();
	alleleReadCount_ifs.close();
	return 0;
}