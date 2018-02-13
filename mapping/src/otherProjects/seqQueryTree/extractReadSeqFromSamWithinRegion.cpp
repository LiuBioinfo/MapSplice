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

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputTotalSam" << endl;
		cout << "#2 chrName" << endl;
		cout << "#3 startPos" << endl;
		cout << "#4 endPos" << endl;
		cout << "#5 outputReadSeqFromSamWithinRegion" << endl;
		exit(1);
	}
	string inputTotalSam = argv[1];
	string chrName = argv[2];
	string startPosStr = argv[3];
	int startPos = atoi(startPosStr.c_str());
	string endPosStr = argv[4];
	int endPos = atoi(endPosStr.c_str());
	string outputReadSeqFromSamWithinRegion = argv[5];

	ifstream rawSam_ifs(inputTotalSam.c_str());
	ofstream seq_ofs(outputReadSeqFromSamWithinRegion.c_str());
	int tmpLineNO = 0;
	int tmpSeqNO = 0;
	while(!rawSam_ifs.eof())
	{
		string samStr;
		getline(rawSam_ifs, samStr);
		if(samStr == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 1000000;
		if(tmpLineNO == tmpThousandIndex * 1000000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		if(samStr.at(0) == '@')
			continue;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 10; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}	
		string tmpReadId = samFieldVec[0];
		string tmpChrName = samFieldVec[2];
		if(tmpChrName != chrName)
			continue;
		string tmpStartPosStr = samFieldVec[3];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		if(tmpStartPos > endPos)
			continue;
		else if(tmpStartPos >= startPos)
		{
			string tmpSeq = samFieldVec[9];
			tmpSeqNO ++;
			seq_ofs << ">" << tmpReadId << "_" << tmpSeqNO << endl << tmpSeq << endl;
			//continue;
		}
		else // tmpStartPos < startPos
		{
			string tmpCigarString = samFieldVec[5];
			int tmpEndPos = getEndMapPos(tmpStartPos, tmpCigarString);
			if(tmpEndPos > endPos)
			{
				string tmpSeq = samFieldVec[9];
				tmpSeqNO ++;
				seq_ofs << ">" << tmpReadId << "_" << tmpSeqNO << endl << tmpSeq << endl;
			}
		}
	}
	seq_ofs.close();;
	rawSam_ifs.close();
	return 0;
}