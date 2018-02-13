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

void parseGtFusionFile2geneNamePairVec(string& fusion_gt_file, 
	//vector< pair<int,int> >& fusionChrNameIntPairVec_gt, 
	//vector< pair<int,int> >& fusionPosPairVec_gt, 
	vector< pair<string, string> >&	fusionGeneNamePairVec_gt, 
	//vector< pair<string, string> >& fusionStrandPairVec_gt,
	GeneAnnEntry_Hash_Info& tmpGeneAnnEntryInfo,
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
		string tmpGeneName_1 = tmpGeneAnnEntryInfo.searchAndReturnGeneName_withGeneId(tmpFusionId_1);
		string tmpGeneName_2 = tmpGeneAnnEntryInfo.searchAndReturnGeneName_withGeneId(tmpFusionId_2);
		fusionGeneNamePairVec_gt.push_back(pair<string,string>(tmpGeneName_1, tmpGeneName_2));
		//fusionChrNameIntPairVec_gt.push_back(pair<int,int>(tmpFusionChrNameInt_1, tmpFusionChrNameInt_2));
		//fusionPosPairVec_gt.push_back(pair<int,int>(tmpFusionPos_1, tmpFusionPos_2));
		//fusionIdPairVec_gt.push_back(pair<string,string>(tmpFusionId_1, tmpFusionId_2));
		//fusionStrandPairVec_gt.push_back(pair<string,string>(tmpFusionStrand_1, tmpFusionStrand_2));
	}
	fusion_gt_ifs.close();
}

bool geneVecPairCmp2gt(vector<string>& gene1Vec_mps3, vector<string>& gene2Vec_mps3,
	string& gene1_gt, string& gene2_gt)
{
	vector< pair<string, string> > genePairVec_mps3;
	for(int tmp1 = 0; tmp1 < gene1Vec_mps3.size(); tmp1++)
	{
		string tmpGene1 = gene1Vec_mps3[tmp1];
		for(int tmp2 = 0; tmp2 < gene2Vec_mps3.size(); tmp2++)
		{
			string tmpGene2 = gene2Vec_mps3[tmp2];
			genePairVec_mps3.push_back(pair<string,string>(tmpGene1, tmpGene2));
		}
	}
	for(int tmp1 = 0; tmp1 < genePairVec_mps3.size(); tmp1++)
	{
		string tmpGene1_mps3 = genePairVec_mps3[tmp1].first;
		string tmpGene2_mps3 = genePairVec_mps3[tmp1].second;
		if(((tmpGene1_mps3 == gene1_gt)&&(tmpGene2_mps3 == gene2_gt))
			||((tmpGene1_mps3 == gene2_gt)&&(tmpGene2_mps3 == gene1_gt)))
			return true;
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputGTF fusion_gt fusion_mps3 outputFolder" << endl;
		exit(1);
	}

	cout << "start to initiate GeneAnnEntryInfoVec and GeneAnnEntryArea2infoIndexMapVec" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();

	string inputGTFentryPrefix = argv[2];
	string inputGTFentryPrefix_index = inputGTFentryPrefix + ".index";
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiateGeneAnnEntryInfoVec(inputGTFentryPrefix, indexInfo);
	tmpGeneAnnHashInfo.initiateGeneAnnEntryArea2infoIndexMapVec(inputGTFentryPrefix_index, indexInfo);
	cout << "end of initiating GeneAnnEntryInfoVec and GeneAnnEntryArea2infoIndexMapVec" << endl;

	cout << "start to load fusion ground truth" << endl;
	string fusion_gt_file = argv[3];
	vector< pair<string, string> > fusionGeneNamePairVec_gt;
	parseGtFusionFile2geneNamePairVec(fusion_gt_file, 
		fusionGeneNamePairVec_gt, tmpGeneAnnHashInfo, indexInfo);
	int gtFusionGenePairNum = fusionGeneNamePairVec_gt.size();
	vector<bool> gtFusionGenePairDetectionBoolVec;
	for(int tmp = 0; tmp < gtFusionGenePairNum; tmp++)
		gtFusionGenePairDetectionBoolVec.push_back(false);

	cout << "end of loading fusion ground truth" << endl;

	cout << "start to create output folder ..." << endl;
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_cmp_file = outputFolderStr + "cmp.txt";
	ofstream cmp_ofs(output_cmp_file.c_str());
	string output_trueFusionGenePair_file = outputFolderStr + "trueFusionGenePair.txt";
	string output_falseFusionGenePair_file = outputFolderStr + "falseFusionGenePair.txt";
	string output_detectedFusionGenePair_file = outputFolderStr + "detectedFusionGenePair.txt";
	string output_missedFusionGenePair_file = outputFolderStr + "missedFusionGenePair.txt";
	ofstream true_ofs(output_trueFusionGenePair_file.c_str());
	ofstream false_ofs(output_falseFusionGenePair_file.c_str());
	ofstream detected_ofs(output_detectedFusionGenePair_file.c_str());
	ofstream missed_ofs(output_missedFusionGenePair_file.c_str());

	int true_num = 0, false_num = 0, detected_num = 0, missed_num;
	string fusion_mps3_file = argv[4];
	ifstream fusion_mps3_ifs(fusion_mps3_file.c_str());
	while(!fusion_mps3_ifs.eof())
	{
		string tmpStr;
		getline(fusion_mps3_ifs, tmpStr);
		if(tmpStr == "")
			break;	
		////////////////////////////
		// start to parse geneVecStr
		////////////////////////////
		vector<string> tmpGeneVec_1, tmpGeneVec_2;		
		int tabLoc = tmpStr.find("\t");
		string tmpGeneVecStr_1 = tmpStr.substr(0, tabLoc);
		string tmpGeneVecStr_2 = tmpStr.substr(tabLoc + 1);
		int startCommaLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabCommaLoc = tmpGeneVecStr_1.find(",", startCommaLoc);
			if(tabCommaLoc == string::npos)
				break;
			else
			{
				string tmpField = tmpGeneVecStr_1.substr(startCommaLoc, tabCommaLoc - startCommaLoc);
				//cout << "tmpGeneId_1: " << tmpField << endl;
				tmpGeneVec_1.push_back(tmpField);
			}
			startCommaLoc = tabCommaLoc + 1;
		}
		startCommaLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabCommaLoc = tmpGeneVecStr_2.find(",", startCommaLoc);
			if(tabCommaLoc == string::npos)
				break;
			else
			{
				string tmpField = tmpGeneVecStr_2.substr(startCommaLoc, tabCommaLoc - startCommaLoc);
				//cout << "tmpGeneId_2: " << tmpField << endl;
				tmpGeneVec_2.push_back(tmpField);
			}
			startCommaLoc = tabCommaLoc + 1;
		}
		////////////////////////////
		// start to cmp 2 ground truth
		////////////////////////////
		bool geneVecPairCmp2gt_bool = false;
		for(int tmp = 0; tmp < fusionGeneNamePairVec_gt.size(); tmp ++)
		{	
			bool tmp_geneVecPairCmp2gt_bool = geneVecPairCmp2gt(tmpGeneVec_1, tmpGeneVec_2, 
				fusionGeneNamePairVec_gt[tmp].first, fusionGeneNamePairVec_gt[tmp].second);
			if(tmp_geneVecPairCmp2gt_bool)
			{	
				gtFusionGenePairDetectionBoolVec[tmp] = true;
				geneVecPairCmp2gt_bool = true;
			}
		}
		if(geneVecPairCmp2gt_bool)
		{
			true_ofs << tmpStr << endl;
			true_num ++;
		}
		else
		{	
			false_ofs << tmpStr << endl;
			false_num ++;
		}
	}
	fusion_mps3_ifs.close();

	for(int tmp = 0; tmp < gtFusionGenePairDetectionBoolVec.size(); tmp++)
	{
		if(gtFusionGenePairDetectionBoolVec[tmp])
		{
			detected_ofs << fusionGeneNamePairVec_gt[tmp].first << "\t" << fusionGeneNamePairVec_gt[tmp].second << endl;
			detected_num ++;
		}
		else
		{
			missed_ofs << fusionGeneNamePairVec_gt[tmp].first << "\t" << fusionGeneNamePairVec_gt[tmp].second << endl;
			missed_num ++;
		}
	}
	
	cmp_ofs << "True fusion gene pair     : " << true_num << endl;
	cmp_ofs << "False fusion gene pair    : " << false_num << endl;
	cmp_ofs << "Detected fusion gene pair : " << detected_num << endl;
	cmp_ofs << "Missed fusion gene pair   : " << missed_num << endl;
	cmp_ofs.close();	
	true_ofs.close();
	false_ofs.close();
	detected_ofs.close();
	missed_ofs.close();
	return 0;
}

