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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../general/learnedCandiSNPhash_info.h"

using namespace std;

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
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "M")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "I")
		{
			tmpEndLocInRead += tmpJumpCodeLength;
		}
		else if(tmpJumpCodeType == "D")
		{
		}
		else if(tmpJumpCodeType == "N")
		{
		}
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

void generateExonLocInReadPosInChr(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
	vector<int>& endLocVecInRead, vector<int>& endPosVecInChr, vector<int>& lenVec)
{
	for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp ++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
		int tmpJumpCodeIndex = tmp;
		if(tmpJumpCodeType == "M")
		{
			int tmpEndLocInRead = getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmpJumpCodeIndex);
			int tmpEndPosInChr = getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec, tmpJumpCodeIndex);
			endLocVecInRead.push_back(tmpEndLocInRead);
			endPosVecInChr.push_back(tmpEndPosInChr);
			lenVec.push_back(tmpJumpCodeLength);
		}
	}
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

bool parseSam2chrNamePosCigarString(string& tmpSamStr, int& chrNameInt, int& chrMapPos, 
	string& cigarString, string& readSeq, int& mismatchNum, int& candiMapOptionNum, 
	Index_Info* indexInfo, bool BeersSamOrNot)
{
	int tabLoc_1 = tmpSamStr.find("\t");
	int tabLoc_2 = tmpSamStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSamStr.find("\t", tabLoc_2 + 1);		
	int tabLoc_4 = tmpSamStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpSamStr.find("\t", tabLoc_4 + 1);
	int tabLoc_6 = tmpSamStr.find("\t", tabLoc_5 + 1);
	int tabLoc_7 = tmpSamStr.find("\t", tabLoc_6 + 1);		
	int tabLoc_8 = tmpSamStr.find("\t", tabLoc_7 + 1);
	int tabLoc_9 = tmpSamStr.find("\t", tabLoc_8 + 1);
	int tabLoc_10 = tmpSamStr.find("\t", tabLoc_9 + 1);
	int tabLoc_11 = tmpSamStr.find("\t", tabLoc_10 + 1);
	int tabLoc_12 = tmpSamStr.find("\t", tabLoc_11 + 1);
	int tabLoc_13 = tmpSamStr.find("\t", tabLoc_12 + 1);
	string tmpChrName = tmpSamStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	chrNameInt = indexInfo->convertStringToInt(tmpChrName);
	if(chrNameInt < 0)
		return false;
	string tmpChrPosStr = tmpSamStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	chrMapPos = atoi(tmpChrPosStr.c_str());
	if(chrMapPos < 0)
		return false;
	cigarString = tmpSamStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
	if(cigarString == "*")
		return false;
	readSeq = tmpSamStr.substr(tabLoc_9 + 1, tabLoc_10 - tabLoc_9 - 1);
	if(BeersSamOrNot)
	{
		string tmpFlagStr = tmpSamStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		if(tmpFlagStr == "0")
		{}
		else if(tmpFlagStr == "16")
		{
			readSeq = getRcmSeq(readSeq);
		}
		else
		{
			cout << "error ! flag in Beers Sam: " << tmpFlagStr << endl;
			cout << "exiting ......" << endl;
			exit(1);
		}
	}
	string tmpNMfieldStr = tmpSamStr.substr(tabLoc_11 + 1, tabLoc_12 - tabLoc_11 - 1);
	string tmpIHfieldStr = tmpSamStr.substr(tabLoc_12 + 1, tabLoc_13 - tabLoc_12 - 1);
	string tmpNMstr = tmpNMfieldStr.substr(5);
	string tmpIHstr = tmpIHfieldStr.substr(5);
	mismatchNum = atoi(tmpNMstr.c_str());
	candiMapOptionNum = atoi(tmpIHstr.c_str());
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputSam outputFile BeersSamOrNot" << endl;
		exit(1);
	}
	int distance2exonBoundary_min = 10;
	int segLen_min = distance2exonBoundary_min * 2 + 1;
	int mismatchNum_max = 2;
	int candiMapOptionNum_max = 1;

	bool BeersSamOrNot = false;
	string BeersSamOrNotStr = argv[4];
	if((BeersSamOrNotStr == "true")||(BeersSamOrNotStr == "TRUE")||(BeersSamOrNotStr == "True"))
		BeersSamOrNot = true;
	else if((BeersSamOrNotStr == "false")||(BeersSamOrNotStr == "FALSE")||(BeersSamOrNotStr == "False"))
		BeersSamOrNot = false;
	else
	{
		cout << "BeersSamOrNot should be set as true or false" << endl;
		cout << "exiting ......" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;	

	cout << "start to initiate learnedCandiSNPhashInfo ....." << endl;
	LearnedCandiSNPhash_Info learnedCandiSNPhashInfo;
	learnedCandiSNPhashInfo.initiate(chromNum);

	cout << "start to extract mismatch from SAM and add 2 learnedCandiSNPhashInfo" << endl;
	string inputSAMpath = argv[2];
	ifstream SAM_ifs(inputSAMpath.c_str());
	int tmpReadNum = 0;
	while(!SAM_ifs.eof())
	{
		string tmpSamStr;
		getline(SAM_ifs, tmpSamStr);
		//cout << "tmpSamStr: " << tmpSamStr << endl;
		if(tmpSamStr == "")
			break;
		if(tmpSamStr.substr(0,1) == "@")
			continue;
		tmpReadNum ++;
		int tmpThousandIndex = tmpReadNum / 100000;
		if(tmpReadNum == tmpThousandIndex * 100000)
			cout << "Processed Read #: " << tmpReadNum << endl;
		int chrNameInt, chrMapPos;
		string cigarString, readSeq;
		int mismatchNum, candiMapOptionNum;
		//cout << "start to parse " << endl;
		bool parseBool = parseSam2chrNamePosCigarString(tmpSamStr, chrNameInt, chrMapPos, cigarString, readSeq, 
			mismatchNum, candiMapOptionNum, indexInfo, BeersSamOrNot);
		//cout << "parseBool: " << parseBool << endl;
		if(!parseBool)
			continue;
		if(mismatchNum > mismatchNum_max)
			continue;
		if(candiMapOptionNum > candiMapOptionNum_max)
			continue;
		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		vector<int> endLocVecInRead;
		vector<int> endPosVecInChr;
		vector<int> lenVec;
		generateExonLocInReadPosInChr(chrMapPos, cigarStringJumpCodeVec, endLocVecInRead, endPosVecInChr, lenVec);		
		int exonNum = lenVec.size();
		//cout << "exonNum: " << exonNum << endl;
		for(int tmp = 0; tmp < exonNum; tmp++)
		{
			int tmpEndLocInRead = endLocVecInRead[tmp];
			int tmpEndPosInChr = endPosVecInChr[tmp];
			int tmpLen = lenVec[tmp];
			//cout << "tmpEndLocInRead: " << tmpEndLocInRead << endl;
			//cout << "tmpEndPosInChr: " << tmpEndPosInChr << endl;
			//cout << "tmpLen: " << tmpLen << endl;
			int tmpStartLocInRead = tmpEndLocInRead - tmpLen + 1;
			int tmpStartPosInChr = tmpEndPosInChr - tmpLen + 1;
			string tmpSeqInRead = readSeq.substr(tmpStartLocInRead - 1, tmpLen);
			string tmpSeqInChr = indexInfo->returnChromStrSubstr(chrNameInt, tmpStartPosInChr, tmpLen);
			//cout << "tmpSeqInRead: " << tmpSeqInRead << endl;
			//cout << "tmpSeqInChr: " << tmpSeqInChr << endl;
			if(tmpLen >= segLen_min)
			{	
				for(int tmpBase = 0 + distance2exonBoundary_min; 
					tmpBase < tmpLen - distance2exonBoundary_min; tmpBase ++)
				{
					if(tmpSeqInRead[tmpBase] != tmpSeqInChr[tmpBase])
					{
						//cout << "tmpBase: " << tmpBase << endl;
						string tmpBaseInRead = tmpSeqInRead.substr(tmpBase, 1);
						// mismatch found
						learnedCandiSNPhashInfo.addMismatchInSam(chrNameInt, tmpStartPosInChr + tmpBase, tmpBaseInRead, indexInfo);
					}
				}
			}
		}
	}

	string outputFile = argv[3];
	learnedCandiSNPhashInfo.outputCandiSNP(outputFile, indexInfo);
	//ofstream candiSNP_ofs(outputFile.c_str());


	//candiSNP_ofs.close();
	free(chrom);
	delete indexInfo;
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	SAM_ifs.close();
	return 0;
}