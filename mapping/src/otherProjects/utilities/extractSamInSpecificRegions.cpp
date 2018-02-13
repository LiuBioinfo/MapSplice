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
#include <sstream>
#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	int jumpCodeIndex)
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

void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpJumpCodeLength;
	string tmpJumpCodeType;
	int jumpCodeStartPosInCigarStr = 0;
	int jumpCodeEndPosInCigarStr;		
	string candidateJumpCodeType = "SMNIDX";
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
			cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
			jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
		}
	}
}

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
	if(argc < 6)
	{
		cout << "Executable inputSam outputSamInSpecficRegions chrName_1 startPos_1 endPos_1 (chrName_2 startPos_2 endPos_2 ...)" << endl;
		exit(1);
	}
	int chrNamePosTriple = (argc - 3)/3;
	if(chrNamePosTriple * 3 + 3 != argc)
	{
		cout << "Error ! chrNamePosTriple * 3 + 3 != argc" << endl;
		exit(1);
	}
	vector<string> chrNameVec_region;
	vector< pair<int,int> > chrPosPairVec_region;
	for(int tmp = 0; tmp < chrNamePosTriple; tmp++)
	{
		string tmpChrName = argv[3 + tmp * 3];
		string tmpChrStartPosStr = argv[4 + tmp * 3];
		int tmpChrStartPos = atoi(tmpChrStartPosStr.c_str());
		string tmpChrEndPosStr = argv[5 + tmp * 3];
		int tmpChrEndPos = atoi(tmpChrEndPosStr.c_str());
		chrNameVec_region.push_back(tmpChrName);
		chrPosPairVec_region.push_back(pair<int,int>(tmpChrStartPos, tmpChrEndPos));
	}

	string inputSam = argv[1];
	string outputSamInSpecficRegions = argv[2];
	ifstream oriSam_ifs(inputSam.c_str());
	ofstream specificRegionSam_ofs(outputSamInSpecficRegions.c_str());
	while(!oriSam_ifs.eof())
	{
		string tmpStr;
		getline(oriSam_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.at(0) == '@')
		{
			specificRegionSam_ofs << tmpStr << endl;
			continue;
		}
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpSam_chrName = tmpFieldVec[2]; // chrName
		string tmpSam_startPosStr = tmpFieldVec[3];
		int tmpSam_startPos = atoi(tmpSam_startPosStr.c_str()); // startPos
		string tmpSam_cigarString = tmpFieldVec[5]; // cigar string
		vector<Jump_Code> tmpSam_jumpCodeVec; // jumpCode vec
		cigarString2jumpCodeVec(tmpSam_cigarString, tmpSam_jumpCodeVec); 
		int tmpSam_endPos = getEndPosOfSpecificJumpCode(tmpSam_startPos, 
			tmpSam_jumpCodeVec, tmpSam_jumpCodeVec.size() - 1); // endPos
		bool tmpSam_inSomeSpecificRegion = false;
		for(int tmpRegionIndex = 0; tmpRegionIndex < chrNameVec_region.size(); tmpRegionIndex ++)
		{
			string tmpRegion_chrName = chrNameVec_region[tmpRegionIndex];
			int tmpRegion_startPos = chrPosPairVec_region[tmpRegionIndex].first;
			int tmpRegion_endPos = chrPosPairVec_region[tmpRegionIndex].second;
			if(tmpRegion_chrName == tmpSam_chrName)
			{
				if(!((tmpRegion_startPos > tmpSam_endPos)||(tmpRegion_endPos < tmpSam_startPos)))
				{
					tmpSam_inSomeSpecificRegion = true;
					break;
				}	
			}
		}
		if(tmpSam_inSomeSpecificRegion)
			specificRegionSam_ofs << tmpStr << endl;
	}
	specificRegionSam_ofs.close();
	oriSam_ifs.close();
	return 0;
}