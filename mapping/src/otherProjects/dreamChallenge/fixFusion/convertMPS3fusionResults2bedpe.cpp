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
#include "../general/geneAnnEntryHash.h"

using namespace std;

time_t nowtime;
struct tm *local;

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

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath validGeneAnnEntryFile inputMPS3fusionResultsFile outputBedPeFile" << endl;
		exit(1);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load indexInfo" << endl;
	string indexFolderPath = argv[1];
	cout << "initiate indexInfo ..." << endl;	
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading indexInfo" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load simplified gene annotation file" << endl;	
	string inputGeneAnnEntryFile = argv[2];	
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading simplified gene annotation file" << endl;

	string inputMPS3fusionResultsFile = argv[3];
	string outputBedPeFile = argv[4];

	ifstream mps3fusionResults_ifs(inputMPS3fusionResultsFile.c_str());
	ofstream bedpe_ofs(outputBedPeFile.c_str());
	while(!mps3fusionResults_ifs.eof())
	{
		string tmpStr;
		getline(mps3fusionResults_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpFusion_chrName_1 = tmpFieldVec[0];
		int tmpFusion_chrNameInt_1 = indexInfo->convertStringToInt(tmpFusion_chrName_1);
		string tmpFusion_chrName_2 = tmpFieldVec[1];
		int tmpFusion_chrNameInt_2 = indexInfo->convertStringToInt(tmpFusion_chrName_2);
		string tmpFusion_posStr_1 = tmpFieldVec[2];
		int tmpFusion_pos_1 = atoi(tmpFusion_posStr_1.c_str());
		string tmpFusion_posStr_2 = tmpFieldVec[3];
		int tmpFusion_pos_2 = atoi(tmpFusion_posStr_2.c_str());
		string tmpFusion_strand_1 = tmpFieldVec[4];
		// if(tmpFusion_strand_1 == "N")
		// 	tmpFusion_strand_1 = ".";
		string tmpFusion_strand_2 = tmpFieldVec[5];
		// if(tmpFusion_strand_2 == "N")
		// 	tmpFusion_strand_2 = ".";
		if((tmpFusion_strand_1 == "N")||(tmpFusion_strand_2 == "N"))
			continue;
		if(((tmpFusion_strand_1 != "+")&&(tmpFusion_strand_1 != "-"))
			||((tmpFusion_strand_2 != "+")&&(tmpFusion_strand_2 != "-")))
			continue;
		string tmpFusion_geneNameVecStr_1 = tmpFieldVec[tmpFieldVec.size() - 2];
		string tmpFusion_geneNameVecStr_2 = tmpFieldVec[tmpFieldVec.size() - 1];
		string tmpFusion_geneId_1 = tmpGeneAnnHashInfo.searchAndReturnSingleGeneId_geneNameVec_strand(
			tmpFusion_geneNameVecStr_1, tmpFusion_strand_1, tmpFusion_chrNameInt_1, tmpFusion_pos_1, indexInfo);
		string tmpFusion_geneId_2 = tmpGeneAnnHashInfo.searchAndReturnSingleGeneId_geneNameVec_strand(
			tmpFusion_geneNameVecStr_2, tmpFusion_strand_2, tmpFusion_chrNameInt_2, tmpFusion_pos_2, indexInfo);
		bedpe_ofs << tmpFusion_chrName_1.substr(3) << "\t" << tmpFusion_pos_1 - 1 << "\t" << tmpFusion_pos_1 
			<< "\t" << tmpFusion_chrName_2.substr(3) << "\t" << tmpFusion_pos_2 - 1 << "\t" << tmpFusion_pos_2
			<< "\t" << tmpFusion_geneId_1 << "-" << tmpFusion_geneId_2 << "\t0\t" << tmpFusion_strand_1 
			<< "\t" << tmpFusion_strand_2 << endl;  
	}
	bedpe_ofs.close();
	mps3fusionResults_ifs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}