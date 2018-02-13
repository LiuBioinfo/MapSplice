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
//#include <omp.h>
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "../../general/index_info.h"

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

void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpJumpCodeLength;
	string tmpJumpCodeType;
	int jumpCodeStartPosInCigarStr = 0;
	int jumpCodeEndPosInCigarStr;
	string candidateJumpCodeType = "SMNIDXH";
	while(1)
	{
		jumpCodeEndPosInCigarStr = 
			jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
		if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
			{break;}
		else
		{
			tmpJumpCodeLength = 
				atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
			tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
			if(tmpJumpCodeType != "H")
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
			jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
		}
	}
}

int getEndLocInReadOfSpecificJumpCode(
	vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	if(jumpCodeIndex < 0)
		return 0;
	int tmpEndLocInRead = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "M")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
			tmpEndLocInRead += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "D")
		{}
		else if(tmpJumpCodeType == "N")
		{}
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return tmpEndLocInRead;
}	

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int tmpEndPos = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{}
		else if(tmpJumpCodeType == "M")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
		{}
		else if(tmpJumpCodeType == "D")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "N")
			tmpEndPos += tmpJumpCodeLength;
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return (tmpEndPos + startPos-1);
}

int getEndMapPos(int startPos, string& jumpCodeStr)
{
	//cout << "start to do getEndMapPos " << endl;
	vector<Jump_Code> cigarStringJumpCodeVec;
	//cout << "start to do cigarString2jumpCodeVec" << endl;
	cigarString2jumpCodeVec(jumpCodeStr, cigarStringJumpCodeVec);
	int jumpCodeIndex_lastMatch = -1;
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
	{
		if(cigarStringJumpCodeVec[tmp].type == "M")
			jumpCodeIndex_lastMatch = tmp;
	}
	//cout << "jumpCodeIndex_lastMatch: " << jumpCodeIndex_lastMatch << endl;
	if(jumpCodeIndex_lastMatch < 0)
	{
		cout << "error, jumpCodeIndex_lastMatch < 0, ==: " << jumpCodeIndex_lastMatch << endl;
		cout << "in checkAlignmentAccuracy.cpp" << endl;
		exit(1);
	}
	return getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec, jumpCodeIndex_lastMatch);
}

void parseBeersSam(string& tmpSamStr_1, string& tmpSamStr_2, string& tmpName_1, string& tmpName_2, 
	int& tmp_PE_chrNameInt, int& tmp_PE_startPos, int& tmp_PE_endPos, string& tmpSeq_1, string& tmpSeq_2, Index_Info* indexInfo)
{
	vector<string> samFieldVec_1;
	vector<string> samFieldVec_2;
	parseStr2fieldVec(samFieldVec_1, tmpSamStr_1);
	parseStr2fieldVec(samFieldVec_2, tmpSamStr_2);
	tmpName_1 = samFieldVec_1[0];
	tmpName_2 = samFieldVec_2[0];
	string tmp_PE_chrNameStr = samFieldVec_1[2];
	tmp_PE_chrNameInt = indexInfo->convertStringToInt(tmp_PE_chrNameStr);
	tmpSeq_1 = samFieldVec_1[9];
	tmpSeq_2 = samFieldVec_2[9];

	string tmpStartPosStr_1 = samFieldVec_1[3];
	string tmpStartPosStr_2 = samFieldVec_2[3];
	string tmpCigarString_1 = samFieldVec_1[5];
	string tmpCigarString_2 = samFieldVec_2[5];

	int tmpStartPos_1 = atoi(tmpStartPosStr_1.c_str());
	int tmpStartPos_2 = atoi(tmpStartPosStr_2.c_str());
	int tmpEndPos_1 = getEndMapPos(tmpStartPos_1, tmpCigarString_1);
	int tmpEndPos_2 = getEndMapPos(tmpStartPos_2, tmpCigarString_2);
	if(tmpStartPos_1 < tmpStartPos_2)
		tmp_PE_startPos = tmpStartPos_1;
	else
		tmp_PE_startPos = tmpStartPos_2;

	if(tmpEndPos_1 < tmpEndPos_2)
		tmp_PE_endPos = tmpEndPos_2;
	else
		tmp_PE_endPos = tmpEndPos_1;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputSam outputFaFilePrefix class_num" << endl;
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

	string class_num_str = argv[4];
	int genomeRegionNum = atoi(class_num_str.c_str()) - 2;
	cout << "genomeSubRegionNum: " << genomeRegionNum << endl;
	unsigned int genomeSubRegionSize = (indexInfo->returnIndexSize())/genomeRegionNum;
	cout << "genomeSubRegionSize: " << genomeSubRegionSize << endl;

	string inputSam = argv[2];
	ifstream sam_ifs(inputSam.c_str());
	string outputFaFilePrefix = argv[3];
	string outputFaFile_end1 = outputFaFilePrefix + ".1.fa";
	string outputFaFile_end2 = outputFaFilePrefix + ".2.fa";
	ofstream end1_fa_ofs(outputFaFile_end1.c_str());
	ofstream end2_fa_ofs(outputFaFile_end2.c_str());
	while(!sam_ifs.eof())
	{
		string tmpSamStr_1;
		getline(sam_ifs, tmpSamStr_1);
		if(tmpSamStr_1 == "")
			break;
		string tmpSamStr_2;
		getline(sam_ifs, tmpSamStr_2);
	
		string tmpName_1, tmpName_2, tmpSeq_1, tmpSeq_2;
		int tmp_PE_chrNameInt, tmp_PE_startPos, tmp_PE_endPos;
		parseBeersSam(tmpSamStr_1, tmpSamStr_2, tmpName_1, tmpName_2, tmp_PE_chrNameInt, 
			tmp_PE_startPos, tmp_PE_endPos, tmpSeq_1, tmpSeq_2, indexInfo);
		unsigned int tmpStartPosInGenome = indexInfo->getWholeGenomeLocation(tmp_PE_chrNameInt, tmp_PE_startPos);
		unsigned int tmpEndPosInGenome = indexInfo->getWholeGenomeLocation(tmp_PE_chrNameInt, tmp_PE_endPos);
		
		int tmpRegionId_startPos = (tmpStartPosInGenome/genomeSubRegionSize) + 1;
		int tmpRegionId_endPos = (tmpEndPosInGenome/genomeSubRegionSize) + 1;
		
		end1_fa_ofs << tmpName_1 << "_" << indexInfo->returnChrNameStr(tmp_PE_chrNameInt) << "_" << tmp_PE_startPos << "_" << tmp_PE_endPos
			<< "_" << tmpRegionId_startPos << "_" << tmpRegionId_endPos << endl << tmpSeq_1 << endl;
		end2_fa_ofs << tmpName_2 << "_" << indexInfo->returnChrNameStr(tmp_PE_chrNameInt) << "_" << tmp_PE_startPos << "_" << tmp_PE_endPos
			<< "_" << tmpRegionId_startPos << "_" << tmpRegionId_endPos << endl << tmpSeq_2 << endl;
	}
	end1_fa_ofs.close();
	end2_fa_ofs.close();
	sam_ifs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}