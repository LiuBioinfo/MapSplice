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
#include "../../../general/splice_info.h"
#include "general/alu_total_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

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

bool parseSam2readNameChrNamePosCigarString(string& tmpSamStr, string& tmpReadName, int& chrNameInt, int& chrMapPos, 
	string& cigarString, Index_Info* indexInfo)
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
	tmpReadName = tmpSamStr.substr(0, tabLoc_1);
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
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{	
		cout << "Executable inputIndexFolderPath inputAluInfoFile inputAluIndexMapSamFile outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "initiate indexInfo ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "end of initiating indexInfo; start to do aluInfo initiating" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of initiating indexInfo; start to do aluInfo initiating" << endl;

	string inputAluInfoFile = argv[2];
	string outputInvalidAluInfoFile = outputFolderStr + "invalid_aluInfo.txt";
	string outputValidAluInfoFile = outputFolderStr + "valid_aluInfo.txt";
	Alu_Total_Info tmpAluTotalInfo;
	tmpAluTotalInfo.initiate_withAluInfoFile(inputAluInfoFile, indexInfo, outputInvalidAluInfoFile, outputValidAluInfoFile);

	vector<int> wholeGenomePos2aluelementVecIndexVec;
	unsigned int wholeGenomeLength = indexInfo->returnGenomeLength();
	for(unsigned int tmp = 0; tmp < wholeGenomeLength; tmp++)
		wholeGenomePos2aluelementVecIndexVec.push_back(-1);
	tmpAluTotalInfo.initiate_wholeGenomePos2aluelementVecIndexVec(wholeGenomePos2aluelementVecIndexVec, indexInfo);

	string readName2aluElementId_path = outputFolderStr + "/read2aluElement.txt";
	ofstream readName2aluElementId_ofs(readName2aluElementId_path.c_str());
	string inputSAMpath = argv[3];
	ifstream SAM_ifs(inputSAMpath.c_str());
	int tmpReadNum = 0;
	while(!SAM_ifs.eof())
	{
		string tmpSamStr;
		getline(SAM_ifs, tmpSamStr);
		if(tmpSamStr == "")
			break;
		if(tmpSamStr.substr(0,1) == "@")
			continue;
		//cout << "tmpSamStr: " << tmpSamStr << endl;
		tmpReadNum ++;
		int tmpThousandIndex = tmpReadNum / 100000;
		if(tmpReadNum == tmpThousandIndex * 100000)
			cout << "Processed Read #: " << tmpReadNum << endl;
		string tmpReadName;
		int chrNameInt, chrMapPos;
		string cigarString;
		bool parseBool = parseSam2readNameChrNamePosCigarString(tmpReadName, tmpSamStr, chrNameInt, chrMapPos, cigarString, indexInfo);
		if(!parseBool)
			continue;
		vector<Jump_Code> cigarStringJumpCodeVec;
		cigarString2jumpCodeVec(cigarString, cigarStringJumpCodeVec);
		int cigarStringJumpCodeVecSize = cigarStringJumpCodeVec.size();
		int chrMapPos_end = getEndPosOfSpecificJumpCode(chrMapPos, cigarStringJumpCodeVec, cigarStringJumpCodeVecSize - 1);
		vector<int> endLocVecInRead;
		vector<int> endPosVecInChr;
		vector<int> lenVec;
		//cout << "start to do generateExonLocInReadPosInChr " << endl;
		set<int> tmpSamAluElementSet;
		generateExonLocInReadPosInChr(chrMapPos, cigarStringJumpCodeVec, endLocVecInRead, endPosVecInChr, lenVec);
		for(int tmp = 0; tmp < endLocVecInRead.size(); tmp++)
		{
			int tmpEndLocInRead = endLocVecInRead[tmp];
			int tmpEndPosInChr = endPosVecInChr[tmp];
			int tmpLen = lenVec[tmp];
			int tmpStartPosInChr = tmpEndPosInChr - tmpLen + 1;
			for(int tmpBase = 0; tmpBase < tmpLen; tmpBase ++)
			{
				int tmpPosInChr = tmpStartPosInChr + tmpBase;
				unsigned int tmpPosInWholeGenome = indexInfo->getWholeGenomeLocation(chrNameInt, tmpPosInChr);
				int tmpAluElementId = wholeGenomePos2aluelementVecIndexVec[tmpPosInWholeGenome];
				tmpSamAluElementSet.insert(tmpAluElementId);
			}
		}
		int tmpSamAluElementSetSize = tmpSamAluElementSet.size();
		if(tmpSamAluElementSetSize > 0)
		{
			readName2aluElementId_ofs << tmpReadName << "\t";
			for(set<int>::iterator tmpSetIter = tmpSamAluElementSet.begin();
				tmpSetIter != tmpSamAluElementSet.end(); tmpSetIter)
			{
				int tmpAluElementIndex = (*tmpSetIter);
				string tmpAluElementId = tmpAluTotalInfo.return_id_withIndex(tmpAluElementIndex);
				readName2aluElementId_ofs << tmpAluElementId << ",";
			}
			readName2aluElementId_ofs << endl;
		}
	}
	SAM_ifs.close();
	readName2aluElementId_ofs.close();
	log_ofs.close();
	free(chrom);
	delete indexInfo;	
	return 0;
}