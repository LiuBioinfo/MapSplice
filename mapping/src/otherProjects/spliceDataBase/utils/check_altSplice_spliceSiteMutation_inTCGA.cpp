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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

typedef map<int,int> JuncEndPos2supNumMap;
typedef map<int, JuncEndPos2supNumMap> JuncStartPos2nextLevelMap;

void extractSampleNameInTCGAjuncFile1stLine(string& tmp1stLineStr, int& totalSampleNum, vector<string>& sampleNameVec)
{
	vector<string> sampleStrFieldVec;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmp1stLineStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpSampleNameField = tmp1stLineStr.substr(startLoc, tabLoc - startLoc);
		sampleStrFieldVec.push_back(tmpSampleNameField);
		startLoc = tabLoc + 1;
	}
	sampleStrFieldVec.push_back(tmp1stLineStr.substr(startLoc));
	for(int tmp = 1; tmp < sampleStrFieldVec.size(); tmp++)
		sampleNameVec.push_back(sampleStrFieldVec[tmp]);
	totalSampleNum = sampleStrFieldVec.size() - 1;
}

bool extractJuncChrNamePosSupNumVecInTCGAjuncFileSpecificLine(string& tmpStr, int& tmpJuncChrNameInt, 
	int& tmpJuncStartPos, int& tmpJuncEndPos, string& tmpJuncStrand, vector<int>& tmpJuncSupNumVec, Index_Info* indexInfo)
{
	vector<string> juncNumInSampleStrFieldVec;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpJuncNumInSampleStrField = tmpStr.substr(startLoc, tabLoc - startLoc);
		juncNumInSampleStrFieldVec.push_back(tmpJuncNumInSampleStrField);
		startLoc = tabLoc + 1;
	}
	juncNumInSampleStrFieldVec.push_back(tmpStr.substr(startLoc));
	string tmpJuncChrNamePosStrand = juncNumInSampleStrFieldVec[0];
	int colonLoc_1 = tmpJuncChrNamePosStrand.find(":");
	int colonLoc_2 = tmpJuncChrNamePosStrand.find(":", colonLoc_1 + 1);
	int comma_1 = tmpJuncChrNamePosStrand.find(",", colonLoc_2 + 1);
	int colonLoc_3 = tmpJuncChrNamePosStrand.find(":", comma_1 + 1);
	int colonLoc_4 = tmpJuncChrNamePosStrand.find(":", colonLoc_3 + 1);
	string tmpJuncChrNameStr = tmpStr.substr(0, colonLoc_1);
	string tmpJuncStartPosStr = tmpStr.substr(colonLoc_1 + 1, colonLoc_2 - colonLoc_1 - 1);
	string tmpJuncEndPosStr = tmpStr.substr(colonLoc_3 + 1, colonLoc_4 - colonLoc_3 - 1);
	tmpJuncChrNameInt = indexInfo->convertStringToInt(tmpJuncChrNameStr);
	if(tmpJuncChrNameInt < 0)
		return false;

	tmpJuncStartPos = atoi(tmpJuncStartPosStr.c_str());
	tmpJuncEndPos = atoi(tmpJuncEndPosStr.c_str());
	tmpJuncStrand = tmpStr.substr(colonLoc_2 + 1, 1);
	for(int tmp = 1; tmp < juncNumInSampleStrFieldVec.size(); tmp++)
	{
		string tmpJuncNumStr = juncNumInSampleStrFieldVec[tmp];
		int tmpJuncNum = atoi(tmpJuncNumStr.c_str());
		tmpJuncSupNumVec.push_back(tmpJuncNum);
	}
	return true;
}

void parseTCGAjuncFile(string& inputJuncFile_TCGA, int& totalSampleNum, int& totalJuncNum, int& validJuncNum, 
	vector<string>& sampleNameVec, vector<int>& juncChrNameIntVec, vector<int>& juncStartPosVec, 
	vector<int>& juncEndPosVec, vector<string>& juncStrandVec, vector< vector<int> >& juncNumInSampleVecVec, Index_Info* indexInfo)
{
	cout << "parseTCGAjuncFile starts ..." << endl;
	totalJuncNum = 0;
	validJuncNum = 0;
	ifstream junc_ifs(inputJuncFile_TCGA.c_str());
	string tmp1stLineStr;
	getline(junc_ifs, tmp1stLineStr);
	extractSampleNameInTCGAjuncFile1stLine(tmp1stLineStr, totalSampleNum, sampleNameVec);
	cout << "totalSampleNum: " << totalSampleNum << endl;
	cout << "sampleNameVec.size(): " << sampleNameVec.size() << endl;
	string tmp2ndLineStr;
	getline(junc_ifs, tmp2ndLineStr);
	while(!junc_ifs.eof())
	{
		string tmpStr;
		getline(junc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpJuncChrNameInt, tmpJuncStartPos, tmpJuncEndPos;
		string tmpJuncStrand;
		vector<int> tmpJuncSupNumVec;
		bool parse_juncStr_bool = extractJuncChrNamePosSupNumVecInTCGAjuncFileSpecificLine(tmpStr, 
			tmpJuncChrNameInt, tmpJuncStartPos, tmpJuncEndPos, tmpJuncStrand, tmpJuncSupNumVec, indexInfo);
		//cout << "parse_juncStr_bool: " << parse_juncStr_bool << endl;
		totalJuncNum ++;
		if(!parse_juncStr_bool)
			continue;
		validJuncNum ++;
		int tmpJuncSupNumVecSize = tmpJuncSupNumVec.size();
		if(tmpJuncSupNumVecSize != totalSampleNum)
		{
			cout << "error: tmpJuncSupNumVecSize: " << tmpJuncSupNumVecSize << endl 
				<< "totalSampleNum: " << totalSampleNum << endl;
			exit(1);
		}
		juncChrNameIntVec.push_back(tmpJuncChrNameInt);
		juncStartPosVec.push_back(tmpJuncStartPos);
		juncEndPosVec.push_back(tmpJuncEndPos);
		juncStrandVec.push_back(tmpJuncStrand);
		juncNumInSampleVecVec.push_back(tmpJuncSupNumVec);
	}
	junc_ifs.close();
}


int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputSNPfile_TCGA inputJuncFile_TCGA outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	log_ofs << "inputIndexPath: " << argv[1] << endl;
	log_ofs << "inputSNPfile_TCGA: " << argv[2] << endl;
	log_ofs << "inputJuncFile_TCGA: " << argv[3] << endl;
	log_ofs << "outputFolder: " << argv[4] << endl;

	string indexFolderPath = argv[1];
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	Index_Info* indexInfo = new Index_Info(indexParameterFileStr);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;
	cout << endl << "[" << asctime(local) << "... start to parse junc files" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to parse junc files" << endl;		

	string inputSNPfile_TCGA = argv[2];
	string inputJuncFile_TCGA = argv[3];

	int totalSampleNum, totalJuncNum, validJuncNum;
	vector<string> sampleNameVec;
	vector<int> juncChrNameIntVec;
	vector<int> juncStartPosVec;
	vector<int> juncEndPosVec;
	vector<string> juncStrandVec;
	vector< vector<int> > juncNumInSampleVecVec;
	parseTCGAjuncFile(inputJuncFile_TCGA, totalSampleNum, totalJuncNum, validJuncNum, sampleNameVec, 
		juncChrNameIntVec, juncStartPosVec, juncEndPosVec, juncStrandVec, juncNumInSampleVecVec, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... start to generate junc file for each sample" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to generate junc file for each sample" << endl;
	cout << "totalSampleNum: " << totalSampleNum << endl;
	cout << "totalJuncNum: " << totalJuncNum << endl;
	cout << "validJuncNum: " << validJuncNum << endl;
	for(int tmpSampleIndex = 0; tmpSampleIndex < totalSampleNum; tmpSampleIndex ++)
	{
		string tmpSampleName = sampleNameVec[tmpSampleIndex];
		string tmpSample_juncFile = outputFolderStr + tmpSampleName + "_junc.txt";
		ofstream tmpSample_junc_ofs(tmpSample_juncFile.c_str());
		for(int tmpJuncIndex = 0; tmpJuncIndex < validJuncNum; tmpJuncIndex ++)
		{
			int tmpJuncChrNameInt = juncChrNameIntVec[tmpJuncIndex];
			string tmpJuncChrNameStr = indexInfo->returnChrNameStr(tmpJuncChrNameInt);
			int tmpJuncStartPos = juncStartPosVec[tmpJuncIndex];
			int tmpJuncEndPos = juncEndPosVec[tmpJuncIndex];
			string tmpJuncStrand = juncStrandVec[tmpJuncIndex];
			int tmpJuncSupNum = (juncNumInSampleVecVec[tmpJuncIndex])[tmpSampleIndex];
			tmpSample_junc_ofs << tmpJuncChrNameStr << "\t" << tmpJuncStartPos << "\t"
				<< tmpJuncEndPos << "\t" << tmpJuncStrand << "\t" << tmpJuncSupNum << endl;
		}
		tmpSample_junc_ofs.close();
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of generating junc file for each sample" << endl;
	log_ofs << endl << "[" << asctime(local) << "... end of generating junc file for each sample" << endl;	
	delete  indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	log_ofs.close();
	return 0;
}