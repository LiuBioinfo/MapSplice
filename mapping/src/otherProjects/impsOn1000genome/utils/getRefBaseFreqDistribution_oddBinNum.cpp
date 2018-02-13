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

string returnBinLabelFromBinIndex(int tmpBinIndex)
{
	switch(tmpBinIndex)
	{
		case 0:
			return "[0,0.025)";
		case 1:
			return "[0.025,0.075)"; 
		case 2:
			return "[0.075,0.125)"; 
		case 3:
			return "[0.125,0.175)"; 
		case 4:
			return "[0.175,0.225)";
		case 5:
			return "[0.225,0.275)"; 
		case 6:
			return "[0.275,0.325)"; 
		case 7:
			return "[0.325,0.375)"; 
		case 8:
			return "[0.375,0.425)";
		case 9:
			return "[0.425,0.475)"; 
		case 10:
			return "[0.475,0.525]";
		case 11:
			return "(0.525,0.575]";
		case 12:
			return "(0.575,0.625]";
		case 13:
			return "(0.625,0.675]";
		case 14:
			return "(0.675,0.725]";
		case 15:
			return "(0.725,0.775]";
		case 16:
			return "(0.775,0.825]";
		case 17:
			return "(0.825,0.875]";
		case 18:
			return "(0.875,0.925]";
		case 19:
			return "(0.925,0.975]";
		case 20:
			return "(0.975,1]";
		default:
			cout << "invalid binIndex: " << tmpBinIndex << endl;
			exit(1);																									
	}
}

int returnBinIndexFromRefFreq_bin21(double tmpRefFreq)
{
	if(tmpRefFreq < 0.025)
		return 0;
	else if((tmpRefFreq >= 0.025)&&(tmpRefFreq < 0.075))
		return 1;
	else if((tmpRefFreq >= 0.075)&&(tmpRefFreq < 0.125))
		return 2;
	else if((tmpRefFreq >= 0.125)&&(tmpRefFreq < 0.175))
		return 3;
	else if((tmpRefFreq >= 0.175)&&(tmpRefFreq < 0.225))
		return 4;
	else if((tmpRefFreq >= 0.225)&&(tmpRefFreq < 0.275))
		return 5;
	else if((tmpRefFreq >= 0.275)&&(tmpRefFreq < 0.325))
		return 6;
	else if((tmpRefFreq >= 0.325)&&(tmpRefFreq < 0.375))
		return 7;
	else if((tmpRefFreq >= 0.375)&&(tmpRefFreq < 0.425))
		return 8;
	else if((tmpRefFreq >= 0.425)&&(tmpRefFreq < 0.475))
		return 9;
	else if((tmpRefFreq >= 0.475)&&(tmpRefFreq <= 0.525))
		return 10;
	else if((tmpRefFreq > 0.525)&&(tmpRefFreq <= 0.575))
		return 11;
	else if((tmpRefFreq > 0.575)&&(tmpRefFreq <= 0.625))
		return 12;
	else if((tmpRefFreq > 0.625)&&(tmpRefFreq <= 0.675))
		return 13;
	else if((tmpRefFreq > 0.675)&&(tmpRefFreq <= 0.725))
		return 14;
	else if((tmpRefFreq > 0.725)&&(tmpRefFreq <= 0.775))
		return 15;
	else if((tmpRefFreq > 0.775)&&(tmpRefFreq <= 0.825))
		return 16;
	else if((tmpRefFreq > 0.825)&&(tmpRefFreq <= 0.875))
		return 17;
	else if((tmpRefFreq > 0.875)&&(tmpRefFreq <= 0.925))
		return 18;
	else if((tmpRefFreq > 0.925)&&(tmpRefFreq <= 0.975))
		return 19;
	else if(tmpRefFreq > 0.975)
		return 20;
	else
	{
		cout << "invalid tmpRefFreq: " << tmpRefFreq << endl;
		exit(1);
	}																															
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 alleleReadCountFile" << endl;
		cout << "#2 outputResults" << endl;
		cout << "#3 minReadCount" << endl;
		//cout << "#4 binSize" << endl;
		exit(1);
	}
	string alleleReadCountFile = argv[1];
	string outputResults = argv[2];
	string outputResults_count = outputResults + ".count";
	string outputResults_stats = outputResults + ".stats";
	string minReadCountStr = argv[3];
	int minReadCount = atoi(minReadCountStr.c_str());
	//string binSizeStr = argv[4];
	//int binSize = atoi(binSizeStr.c_str());

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
	for(int tmp = 0; tmp <= 20; tmp++)
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

		//int refFreq_mps, refFreq_imps, refFreq_mps_mask, refFreq_hisat2genome, 
		//	refFreq_hisat2iSNP, refFreq_hisat2pSNP, refFreq_hisat2mask, refFreq_star, refFreq_star_mask;

		if(count_mps >= minReadCount)
		{
			double refFreq_mps = (double)refCount_mps/(double)count_mps;
			//cout << "refFreq_mps: " << refFreq_mps << endl;
			//cout << "refCount_mps: " << refCount_mps << endl;
			//cout << "count_mps: " << count_mps << endl;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_mps);
			refReadCountVec_mps[tmpIndex] ++;
			totalCount_mps ++;
		}
		if(count_imps >= minReadCount)
		{
			double refFreq_imps = (double)refCount_imps/(double)count_imps;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_imps);
			refReadCountVec_imps[tmpIndex] ++;
			totalCount_imps ++;
		}
		if(count_mps_mask >= minReadCount)
		{
			double refFreq_mps_mask = (double)refCount_mps_mask/(double)count_mps_mask;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_mps_mask);
			refReadCountVec_mps_mask[tmpIndex] ++;
			totalCount_mps_mask ++;
		}

		if(count_hisat2genome >= minReadCount)
		{
			double refFreq_hisat2genome = (double)refCount_hisat2genome/(double)count_hisat2genome;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_hisat2genome);
			refReadCountVec_hisat2genome[tmpIndex] ++;
			totalCount_hisat2genome ++;
		}
		if(count_hisat2iSNP >= minReadCount)
		{
			double refFreq_hisat2iSNP = (double)refCount_hisat2iSNP/(double)count_hisat2iSNP;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_hisat2iSNP);
			refReadCountVec_hisat2iSNP[tmpIndex] ++;
			totalCount_hisat2iSNP ++;
		}
		if(count_hisat2pSNP >= minReadCount)
		{
			double refFreq_hisat2pSNP = (double)refCount_hisat2pSNP/(double)count_hisat2pSNP;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_hisat2pSNP);
			refReadCountVec_hisat2pSNP[tmpIndex] ++;
			totalCount_hisat2pSNP ++;
		}
		if(count_hisat2mask >= minReadCount)
		{
			double refFreq_hisat2mask = (double)refCount_hisat2mask/(double)count_hisat2mask;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_hisat2mask);
			refReadCountVec_hisat2mask[tmpIndex] ++;
			totalCount_hisat2mask ++;
		}

		if(count_star >= minReadCount)
		{
			double refFreq_star = (double)refCount_star / (double)count_star;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_star);
			refReadCountVec_star[tmpIndex] ++;
			totalCount_star ++;
		}
		if(count_star_mask >= minReadCount)
		{
			double refFreq_star_mask = (double)refCount_star_mask / (double)count_star_mask;
			int tmpIndex = returnBinIndexFromRefFreq_bin21(refFreq_star_mask);
			refReadCountVec_star_mask[tmpIndex] ++;
			totalCount_star_mask ++;
		}

	}

	results_ofs << "Aligner\tMapSplice\tiMapSplice\tMapSplice Mask\tHISAT2\tHISAT2 SNP\tHISAT2 POP\tHISAT2 MASK\tSTAR\tSTAR MASK" << endl;
	//results_ofs << "Total\t1\t1\t1\t1\t1\t1\t1\t1\t1" << endl;
	for(int tmp = 0; tmp <= 20; tmp++)
		results_ofs << returnBinLabelFromBinIndex(tmp) << "\t" 
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
	for(int tmp = 0; tmp <= 20; tmp++)
		results_ofs_count << returnBinLabelFromBinIndex(tmp) << "\t" 
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