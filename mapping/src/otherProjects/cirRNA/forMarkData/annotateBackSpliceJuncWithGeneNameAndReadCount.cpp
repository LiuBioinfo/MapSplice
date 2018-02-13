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
#include "../../dreamChallenge/general/geneAnnEntryHash.h"
using namespace std;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

void parseGeneNameGeneCountVecFromFeatureCountPastedStr(string& tmpFeatureCountPastedStr, string& tmpGeneName, vector<int>& tmpGeneCountVec)
{
	vector<string> tmpFieldVec;
	parseStr2fieldVec(tmpFieldVec, tmpFeatureCountPastedStr);
	tmpGeneName = tmpFieldVec[0];
	int countVecSize = tmpFieldVec.size() / 2;
	for(int tmp = 0; tmp < countVecSize; tmp++)
	{
		int tmpIndex = tmp * 2 + 1;
		int tmpGeneCount = atoi(tmpFieldVec[tmpIndex].c_str());
		tmpGeneCountVec.push_back(tmpGeneCount);
	}
}

void copyGeneCountVec2target_withGeneName(string& tmpQueryGeneName, vector< pair<string, vector<int> > >& geneNameCountVecVec, vector<int>& targetGeneCountVec)
{
	for(int tmp1 = 0; tmp1 < geneNameCountVecVec.size(); tmp1++)
	{
		string tmpGeneName = geneNameCountVecVec[tmp1].first;
		if(tmpQueryGeneName == tmpGeneName)
		{
			int tmpGeneCountVecSize = (geneNameCountVecVec[tmp1].second).size();
			for(int tmp2 = 0; tmp2 < tmpGeneCountVecSize; tmp2++)
				targetGeneCountVec.push_back((geneNameCountVecVec[tmp1].second)[tmp2]);
			return;
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputGeneAnnEntry inputFeatureCountPastedFile inputBackSJlist outputFile" << endl;
		exit(1);
	}
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	// load gene ann entry file
	string inputGeneAnnEntryFile = argv[2];
	cout << "start to load simplified gene annotation file" << endl;
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	cout << "start to initiate GeneAnnEntryArea2infoIndexMapVec" << endl;
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);

	// load feature count pasted file
	vector< pair<string, vector<int> > > geneNameCountVecVec;
	string inputFeatureCountPastedFile = argv[3];
	ifstream fc_ifs(inputFeatureCountPastedFile.c_str());
	while(!fc_ifs.eof())
	{
		string tmpStr;
		getline(fc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpGeneName;
		vector<int> tmpGeneCountVec;
		parseGeneNameGeneCountVecFromFeatureCountPastedStr(tmpStr, tmpGeneName, tmpGeneCountVec);
		geneNameCountVecVec.push_back(pair<string, vector<int> >(tmpGeneName, tmpGeneCountVec));
	}
	fc_ifs.close();

	string outputFile = argv[5];
	ofstream annotated_ofs(outputFile.c_str());
	string outputFile_invalid = outputFile + ".invalid";
	ofstream invalid_ofs(outputFile_invalid.c_str());

	string inputBackSJlistFile = argv[4];
	ifstream backSJ_ifs(inputBackSJlistFile.c_str());
	while(!backSJ_ifs.eof())
	{
		string tmpStr;
		getline(backSJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpBackSJ_startPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpBackSJ_endPosStr;
		if(tabLoc_3 == string::npos)
			tmpBackSJ_endPosStr = tmpStr.substr(tabLoc_2 + 1);
		else
			tmpBackSJ_endPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpBackSJ_chrName = tmpStr.substr(0, tabLoc_1);
		int tmpBackSJ_chrNameInt = indexInfo->convertStringToInt(tmpBackSJ_chrName);
		int tmpBackSJ_startPos = atoi(tmpBackSJ_startPosStr.c_str());
		int tmpBackSJ_endPos = atoi(tmpBackSJ_endPosStr.c_str());
		string tmpBackSJ_geneName_startPos = tmpGeneAnnHashInfo.searchEndPosAndReturn1stGeneName(tmpBackSJ_chrNameInt, tmpBackSJ_startPos);
		string tmpBackSJ_geneName_endPos = tmpGeneAnnHashInfo.searchStartPosAndReturn1stGeneName(tmpBackSJ_chrNameInt, tmpBackSJ_endPos);
		if((tmpBackSJ_geneName_startPos != tmpBackSJ_geneName_endPos)||(tmpBackSJ_geneName_startPos == "NULL"))
		{
			invalid_ofs << tmpStr << "\t" << tmpBackSJ_geneName_startPos << "\t" << tmpBackSJ_geneName_endPos << endl;
			continue;
		}
		string tmpBackSJ_geneName = tmpBackSJ_geneName_startPos;
		vector<int> tmpBackSJ_geneCountVec;
		copyGeneCountVec2target_withGeneName(tmpBackSJ_geneName, geneNameCountVecVec, tmpBackSJ_geneCountVec);
		annotated_ofs << tmpStr << "\t" << tmpBackSJ_geneName;
		for(int tmp = 0; tmp < tmpBackSJ_geneCountVec.size(); tmp++)
			annotated_ofs << "\t" << tmpBackSJ_geneCountVec[tmp];
		annotated_ofs << endl;
	}
	backSJ_ifs.close();
	annotated_ofs.close();
	invalid_ofs.close();
	return 0;
}