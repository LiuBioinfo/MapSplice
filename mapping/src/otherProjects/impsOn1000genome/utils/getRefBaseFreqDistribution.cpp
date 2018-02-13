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
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 alleleReadCountFile" << endl;
		cout << "#2 outputResults" << endl;
		cout << "#3 minReadCount" << endl;
		cout << "#4 binSize" << endl;
		exit(1);
	}
	string alleleReadCountFile = argv[1];
	string outputResults = argv[2];
	string outputResults_count = outputResults + ".count";
	string outputResults_stats = outputResults + ".stats";
	string minReadCountStr = argv[3];
	int minReadCount = atoi(minReadCountStr.c_str());
	string binSizeStr = argv[4];
	int binSize = atoi(binSizeStr.c_str());

	int totalCount_mps = 0;
	int totalCount_imps = 0;
	int totalCount_mps_mask = 0;
	int totalCount_hisat2genome = 0;
	int totalCount_hisat2iSNP = 0;
	int totalCount_hisat2pSNP = 0;
	int totalCount_hisat2mask = 0;
	int totalCount_star = 0;
	int totalCount_star_mask = 0;
	vector<int> refReadCountVec_mps;
	vector<int> refReadCountVec_imps;
	vector<int> refReadCountVec_mps_mask;
	vector<int> refReadCountVec_hisat2genome;
	vector<int> refReadCountVec_hisat2iSNP;
	vector<int> refReadCountVec_hisat2pSNP;
	vector<int> refReadCountVec_hisat2mask;
	vector<int> refReadCountVec_star;
	vector<int> refReadCountVec_star_mask;
	for(int tmp = 0; tmp <= (100/binSize); tmp++)
	{
		refReadCountVec_mps.push_back(0);
		refReadCountVec_imps.push_back(0);
		refReadCountVec_mps_mask.push_back(0);
		refReadCountVec_hisat2genome.push_back(0);
		refReadCountVec_hisat2iSNP.push_back(0);
		refReadCountVec_hisat2pSNP.push_back(0);
		refReadCountVec_hisat2mask.push_back(0);
		refReadCountVec_star.push_back(0);
		refReadCountVec_star_mask.push_back(0);
	}

	ifstream alleleReadCount_ifs(alleleReadCountFile.c_str());
	ofstream results_ofs(outputResults.c_str());
	ofstream results_ofs_count(outputResults_count.c_str());
	ofstream stats_ofs(outputResults_stats.c_str());
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
			refFreq_mps = (100/binSize) * refCount_mps / count_mps;
			refReadCountVec_mps[refFreq_mps] ++;
			totalCount_mps ++;
		}
		if(count_imps >= minReadCount)
		{
			refFreq_imps = (100/binSize) * refCount_imps / count_imps;
			refReadCountVec_imps[refFreq_imps] ++;
			totalCount_imps ++;
		}
		if(count_mps_mask >= minReadCount)
		{
			refFreq_mps_mask = (100/binSize) * refCount_mps_mask / count_mps_mask;
			refReadCountVec_mps_mask[refFreq_mps_mask] ++;
			totalCount_mps_mask ++;
		}

		if(count_hisat2genome >= minReadCount)
		{
			refFreq_hisat2genome = (100/binSize) * refCount_hisat2genome / count_hisat2genome;
			refReadCountVec_hisat2genome[refFreq_hisat2genome] ++;
			totalCount_hisat2genome ++;
		}
		if(count_hisat2iSNP >= minReadCount)
		{
			refFreq_hisat2iSNP = (100/binSize) * refCount_hisat2iSNP / count_hisat2iSNP;
			refReadCountVec_hisat2iSNP[refFreq_hisat2iSNP] ++;
			totalCount_hisat2iSNP ++;
		}
		if(count_hisat2pSNP >= minReadCount)
		{
			refFreq_hisat2pSNP = (100/binSize) * refCount_hisat2pSNP / count_hisat2pSNP;
			refReadCountVec_hisat2pSNP[refFreq_hisat2pSNP] ++;
			totalCount_hisat2pSNP ++;
		}
		if(count_hisat2mask >= minReadCount)
		{
			refFreq_hisat2mask = (100/binSize) * refCount_hisat2mask / count_hisat2mask;
			refReadCountVec_hisat2mask[refFreq_hisat2mask] ++;
			totalCount_hisat2mask ++;
		}

		if(count_star >= minReadCount)
		{
			refFreq_star = (100/binSize) * refCount_star / count_star;
			refReadCountVec_star[refFreq_star] ++;
			totalCount_star ++;
		}
		if(count_star_mask >= minReadCount)
		{
			refFreq_star_mask = (100/binSize) * refCount_star_mask / count_star_mask;
			refReadCountVec_star_mask[refFreq_star_mask] ++;
			totalCount_star_mask ++;
		}

	}

	results_ofs << "Aligner\tMapSplice\tiMapSplice\tMapSplice Mask\tHISAT2\tHISAT2 SNP\tHISAT2 POP\tHISAT2 MASK\tSTAR\tSTAR MASK" << endl;
	//results_ofs << "Total\t1\t1\t1\t1\t1\t1\t1\t1\t1" << endl;
	for(int tmp = 0; tmp <= (100/binSize); tmp++)
		results_ofs << (double)tmp/(double)binSize << "\t" 
			<< ((double)refReadCountVec_mps[tmp]/(double)totalCount_mps)*100 << "\t"
			<< ((double)refReadCountVec_imps[tmp]/(double)totalCount_imps)*100 << "\t"
			<< ((double)refReadCountVec_mps_mask[tmp]/(double)totalCount_mps_mask)*100 << "\t"
			<< ((double)refReadCountVec_hisat2genome[tmp]/(double)totalCount_hisat2genome)*100 << "\t"
			<< ((double)refReadCountVec_hisat2iSNP[tmp]/(double)totalCount_hisat2iSNP)*100 << "\t"
			<< ((double)refReadCountVec_hisat2pSNP[tmp]/(double)totalCount_hisat2pSNP)*100 << "\t"
			<< ((double)refReadCountVec_hisat2mask[tmp]/(double)totalCount_hisat2mask)*100 << "\t"
			<< ((double)refReadCountVec_star[tmp]/(double)totalCount_star)*100 << "\t"
			<< ((double)refReadCountVec_star_mask[tmp]/(double)totalCount_star_mask)*100 << endl;

	results_ofs_count << "Aligner\tMapSplice\tiMapSplice\tMapSplice Mask\tHISAT2\tHISAT2 SNP\tHISAT2 POP\tHISAT2 MASK\tSTAR\tSTAR MASK" << endl;
	for(int tmp = 0; tmp <= (100/binSize); tmp++)
		results_ofs_count << (double)tmp*((double)binSize/100) << "\t" 
			<< refReadCountVec_mps[tmp] << "\t" << refReadCountVec_imps[tmp] << "\t" << refReadCountVec_mps_mask[tmp] << "\t"
			<< refReadCountVec_hisat2genome[tmp] << "\t" << refReadCountVec_hisat2iSNP[tmp] << "\t"
			<< refReadCountVec_hisat2pSNP[tmp] << "\t" << refReadCountVec_hisat2mask[tmp] << "\t"
			<< refReadCountVec_star[tmp] << "\t" << refReadCountVec_star_mask[tmp] << endl;

	stats_ofs << "Aligner\tMapSplice\tiMapSplice\tMapSplice Mask\tHISAT2\tHISAT2 SNP\tHISAT2 POP\tHISAT2 MASK\tSTAR\tSTAR MASK" << endl;
	stats_ofs << "Total\t" << totalCount_mps << "\t" << totalCount_imps << "\t" << totalCount_mps_mask << "\t"
		<< totalCount_hisat2genome << "\t" << totalCount_hisat2iSNP << "\t" << totalCount_hisat2pSNP << "\t" << totalCount_hisat2mask 
		<< "\t" << totalCount_star << "\t" << totalCount_star_mask << endl;

	results_ofs.close();
	results_ofs_count.close();
	stats_ofs.close();
	alleleReadCount_ifs.close();
	return 0;
}