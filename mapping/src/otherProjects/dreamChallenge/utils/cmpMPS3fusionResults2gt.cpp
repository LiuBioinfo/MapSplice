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
#include "../../../general/index_info.h"
#include "../general/geneAnnEntryHash.h"

using namespace std;

void parseGtFusionStr(string& tmpGtFusionStr, int& fusionChrNameInt_1,
	int& fusionChrNameInt_2, int& fusionPos_1, int& fusionPos_2,
	string& fusionId_1, string& fusionId_2, string& fusionStrand_1,
	string& fusionStrand_2, Index_Info* indexInfo)
{
	vector<string> tmpFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 7; tmp++)
	{
		int tabLoc = tmpGtFusionStr.find("\t", startLoc);
		string tmpField = tmpGtFusionStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpGtFusionStr.substr(startLoc));
	string tmpChrName_1 = tmpFieldVec[0];
	string tmpChrPosStr_1 = tmpFieldVec[1];
	string tmpChrName_2 = tmpFieldVec[2];
	string tmpChrPosStr_2 = tmpFieldVec[3];
	fusionChrNameInt_1 = indexInfo->convertStringToInt(tmpChrName_1);
	fusionChrNameInt_2 = indexInfo->convertStringToInt(tmpChrName_2);	
	fusionPos_1 = atoi(tmpChrPosStr_1.c_str());
	fusionPos_2 = atoi(tmpChrPosStr_2.c_str());
	fusionId_1 = tmpFieldVec[4]; 
	fusionId_2 = tmpFieldVec[5];
	fusionStrand_1 = tmpFieldVec[6];
	fusionStrand_2 = tmpFieldVec[7];
}

void parseGtFusionFile2pairVec(string& fusion_gt_file, 
	vector< pair<int,int> >& fusionChrNameIntPairVec_gt, 
	vector< pair<int,int> >& fusionPosPairVec_gt, 
	vector< pair<string, string> >&	fusionIdPairVec_gt, 
	vector< pair<string, string> >& fusionStrandPairVec_gt, 
	Index_Info* indexInfo)
{
	ifstream fusion_gt_ifs(fusion_gt_file.c_str());
	while(!fusion_gt_ifs.eof())
	{
		string tmpStr;
		getline(fusion_gt_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpFusionChrNameInt_1, tmpFusionChrNameInt_2, tmpFusionPos_1, tmpFusionPos_2;
		string tmpFusionId_1, tmpFusionId_2, tmpFusionStrand_1, tmpFusionStrand_2;
		parseGtFusionStr(tmpStr, tmpFusionChrNameInt_1, tmpFusionChrNameInt_2, 
			tmpFusionPos_1, tmpFusionPos_2, tmpFusionId_1, tmpFusionId_2, 
			tmpFusionStrand_1, tmpFusionStrand_2, indexInfo);
		fusionChrNameIntPairVec_gt.push_back(pair<int,int>(tmpFusionChrNameInt_1, tmpFusionChrNameInt_2));
		fusionPosPairVec_gt.push_back(pair<int,int>(tmpFusionPos_1, tmpFusionPos_2));
		fusionIdPairVec_gt.push_back(pair<string,string>(tmpFusionId_1, tmpFusionId_2));
		fusionStrandPairVec_gt.push_back(pair<string,string>(tmpFusionStrand_1, tmpFusionStrand_2));
	}
	fusion_gt_ifs.close();
}

void parseMPS3fusion(string& tmpStr, int& tmpChrNameInt_1, int& tmpPos_1, 
	int& tmpChrNameInt_2, int& tmpPos_2, Index_Info* indexInfo)
{
	vector<string> tmpFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 4; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	string tmpChrName_1 = tmpFieldVec[0];
	string tmpChrName_2 = tmpFieldVec[1];
	string tmpPosStr_1 = tmpFieldVec[2];
	string tmpPosStr_2 = tmpFieldVec[3];
	tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrName_1);
	tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrName_2);
	tmpPos_1 = atoi(tmpPosStr_1.c_str());
	tmpPos_2 = atoi(tmpPosStr_2.c_str());
}

void generateUniqGeneIdVec(vector<string>& tmpEntryGeneIdVec, vector<string>& tmpUniqGeneIdVec)
{
	int tmpEntryGeneIdVecSize = tmpEntryGeneIdVec.size();
	for(int tmp1 = 0; tmp1 < tmpEntryGeneIdVecSize; tmp1++)
	{
		string tmpEntryGeneId = tmpEntryGeneIdVec[tmp1];
		int tmpCurrentUniqGeneIdVecSize = tmpUniqGeneIdVec.size();
		bool alreadyInUniqGeneIdVec_bool = false;
		for(int tmp2 = 0; tmp2 < tmpCurrentUniqGeneIdVecSize; tmp2++)
		{
			string tmpUniqGeneId = tmpUniqGeneIdVec[tmp2];
			if(tmpUniqGeneId == tmpEntryGeneId)
			{
				alreadyInUniqGeneIdVec_bool = true;
				break;
			}
		}
		if(!alreadyInUniqGeneIdVec_bool)
			tmpUniqGeneIdVec.push_back(tmpEntryGeneId);
	}
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexPath inputGTF fusion_gt fusion_mps3 outputFolder" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	cout << "initiate indexInfo ..." << endl;	
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();
	cout << "end of initiating indexInfo" << endl;

	cout << "start to initiate GeneAnnEntryInfoVec and GeneAnnEntryArea2infoIndexMapVec" << endl;
	string inputGTFentryPrefix = argv[2];
	string inputGTFentryPrefix_index = inputGTFentryPrefix + ".index";
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiateGeneAnnEntryInfoVec(inputGTFentryPrefix, indexInfo);
	tmpGeneAnnHashInfo.initiateGeneAnnEntryArea2infoIndexMapVec(inputGTFentryPrefix_index, indexInfo);
	cout << "end of initiating GeneAnnEntryInfoVec and GeneAnnEntryArea2infoIndexMapVec" << endl;

	cout << "start to load fusion ground truth" << endl;
	string fusion_gt_file = argv[3];
	vector< pair<int,int> > fusionChrNameIntPairVec_gt;
	vector< pair<int,int> > fusionPosPairVec_gt;
	vector< pair<string, string> > fusionIdPairVec_gt;
	vector< pair<string, string> > fusionStrandPairVec_gt;
	parseGtFusionFile2pairVec(fusion_gt_file, fusionChrNameIntPairVec_gt, fusionPosPairVec_gt, 
		fusionIdPairVec_gt, fusionStrandPairVec_gt, indexInfo);
	int gtFusionGenePairNum = fusionChrNameIntPairVec_gt.size();
	vector<bool> gtFusionGenePairDetectionBoolVec;
	vector<bool> gtFusionPosPairDetectionBoolVec;
	for(int tmp = 0; tmp < gtFusionGenePairNum; tmp++)
	{
		gtFusionGenePairDetectionBoolVec.push_back(false);
		gtFusionPosPairDetectionBoolVec.push_back(false);
	}
	cout << "end of loading fusion ground truth" << endl;

	cout << "start to create output folder ..." << endl;
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	string mps3_fusion_cmp2gt_file_gene = outputFolderStr + "gt_cmp_gene.txt";
	string mps3_fusion_cmp2gt_file_pos = outputFolderStr + "gt_cmp_pos.txt";
	string mps3_fusion_trueGenePair_file = outputFolderStr + "trueGene.txt";
	string mps3_fusion_trueGenePair_falsePosPair_file = outputFolderStr + "trueGene_falsePos.txt";
	string mps3_fusion_falseGenePair_file = outputFolderStr + "falseGene.txt";
	ofstream mps3_fusion_cmp2gt_ofs_gene(mps3_fusion_cmp2gt_file_gene.c_str());
	ofstream mps3_fusion_cmp2gt_ofs_pos(mps3_fusion_cmp2gt_file_pos.c_str());
	ofstream mps3_fusion_trueGenePair_ofs(mps3_fusion_trueGenePair_file.c_str());
	ofstream mps3_fusion_trueGenePair_falsePosPair_ofs(mps3_fusion_trueGenePair_falsePosPair_file.c_str());
	ofstream mps3_fusion_falseGenePair_ofs(mps3_fusion_falseGenePair_file.c_str());
	cout << "start to load fusions detected by mps3, and compare 2 ground truth" << endl;
	int totalPosPairNum = 0, truePosPairNum = 0,  falsePosPairNum_falseGene = 0;
	string fusion_mps3_file = argv[4];
	ifstream fusion_mps3_ifs(fusion_mps3_file.c_str());
	while(!fusion_mps3_ifs.eof())
	{
		string tmpStr;
		getline(fusion_mps3_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpChrNameInt_1, tmpPos_1, tmpChrNameInt_2, tmpPos_2;
		parseMPS3fusion(tmpStr, tmpChrNameInt_1, tmpPos_1, tmpChrNameInt_2, tmpPos_2, indexInfo);
		vector<string> tmpEntryGeneIdVec_1, tmpEntryGeneIdVec_2;
		tmpGeneAnnHashInfo.searchAndReturnGeneAnnEntryGeneIdVec(tmpEntryGeneIdVec_1, tmpChrNameInt_1, tmpPos_1, indexInfo);
		tmpGeneAnnHashInfo.searchAndReturnGeneAnnEntryGeneIdVec(tmpEntryGeneIdVec_2, tmpChrNameInt_2, tmpPos_2, indexInfo);		
		vector<string> tmpUniqGeneIdVec_1, tmpUniqGeneIdVec_2;
		generateUniqGeneIdVec(tmpEntryGeneIdVec_1, tmpUniqGeneIdVec_1);
		generateUniqGeneIdVec(tmpEntryGeneIdVec_2, tmpUniqGeneIdVec_2);
		bool tmpFusionGenePair_true_or_false_bool = false;
		bool tmpFusionPosPair_true_or_false_bool = false;
		int tmpUniqGeneIdVec_1_size = tmpUniqGeneIdVec_1.size();
		int tmpUniqGeneIdVec_2_size = tmpUniqGeneIdVec_2.size();
		for(int tmp1 = 0; tmp1 < tmpUniqGeneIdVec_1_size; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < tmpUniqGeneIdVec_2_size; tmp2++)
			{
				string tmpGeneId_1 = tmpUniqGeneIdVec_1[tmp1];
				string tmpGeneId_2 = tmpUniqGeneIdVec_2[tmp2];
				for(int tmp3 = 0; tmp3 < gtFusionGenePairNum; tmp3++)
				{
					if(((tmpGeneId_1 == fusionIdPairVec_gt[tmp3].first)&&(tmpGeneId_2 == fusionIdPairVec_gt[tmp3].second))
						||((tmpGeneId_2 == fusionIdPairVec_gt[tmp3].first)&&(tmpGeneId_1 == fusionIdPairVec_gt[tmp3].second)))
					{
						tmpFusionGenePair_true_or_false_bool = true;
						gtFusionGenePairDetectionBoolVec[tmp3] = true;
						if(((tmpPos_1 == fusionPosPairVec_gt[tmp3].first)&&(tmpPos_2 == fusionPosPairVec_gt[tmp3].second))
							||((tmpPos_2 == fusionPosPairVec_gt[tmp3].first)&&(tmpPos_1 == fusionPosPairVec_gt[tmp3].second)))
						{
							tmpFusionPosPair_true_or_false_bool = true;
							gtFusionPosPairDetectionBoolVec[tmp3] = true;
						}
					}
				}
			}
		}
		if(tmpFusionGenePair_true_or_false_bool)
		{
			if(tmpFusionPosPair_true_or_false_bool)
				mps3_fusion_trueGenePair_ofs << tmpStr << endl;
			else
				mps3_fusion_trueGenePair_falsePosPair_ofs << tmpStr << endl;
		}
		else
			mps3_fusion_falseGenePair_ofs << tmpStr << endl;
	}
	int detectedGtFusionGenePairNum = 0;
	int detectedGtFusionPosPairNum = 0;
	for(int tmp = 0; tmp < gtFusionGenePairNum; tmp++)
	{
		mps3_fusion_cmp2gt_ofs_gene << fusionIdPairVec_gt[tmp].first << "-" << fusionIdPairVec_gt[tmp].second << "\t";
		mps3_fusion_cmp2gt_ofs_pos << fusionIdPairVec_gt[tmp].first << "-" << fusionIdPairVec_gt[tmp].second << "\t";			
		if(gtFusionGenePairDetectionBoolVec[tmp])
		{
			detectedGtFusionGenePairNum ++;
			mps3_fusion_cmp2gt_ofs_gene << "Y" << endl;
		}
		else
			mps3_fusion_cmp2gt_ofs_gene << "N" << endl;

		if(gtFusionPosPairDetectionBoolVec[tmp])
		{
			detectedGtFusionPosPairNum ++;
			mps3_fusion_cmp2gt_ofs_pos << "Y" << endl;
		}
		else
			mps3_fusion_cmp2gt_ofs_pos << "N" << endl;
	}

	double detectedGtFusionGenePairPerc = 100 * (double)detectedGtFusionGenePairNum/(double)gtFusionGenePairNum; 
	double missedGtFusionGenePairPerc = 100 - detectedGtFusionGenePairPerc;
	log_ofs << "Total fusion gene pair: " << gtFusionGenePairNum << endl;
	log_ofs << "Detected fusion gene pair: " << detectedGtFusionGenePairNum << " -- " << detectedGtFusionGenePairPerc << "%" << endl;
	log_ofs << "Missed fusion gene pair: " << gtFusionGenePairNum - detectedGtFusionGenePairNum
		<< " -- " << missedGtFusionGenePairPerc << "%" << endl;
	fusion_mps3_ifs.close();
	cout << "end of loading fusions detected by mps3, and compare 2 ground truth" << endl;
	mps3_fusion_cmp2gt_ofs_gene.close();
	mps3_fusion_cmp2gt_ofs_pos.close();
	mps3_fusion_trueGenePair_ofs.close();
	mps3_fusion_trueGenePair_falsePosPair_ofs.close();
	mps3_fusion_falseGenePair_ofs.close();
	log_ofs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}