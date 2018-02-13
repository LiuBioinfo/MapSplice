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
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "../../general/index_info.h"

using namespace std;

int getReadNOfromReadName(string& tmpReadNameStr)
{
	string readNOstr;
	int dotLoc = tmpReadNameStr.find(".");
	int slashLoc = tmpReadNameStr.find("/");
	if(slashLoc == string::npos)
		readNOstr = tmpReadNameStr.substr(dotLoc + 1);
	else
		readNOstr = tmpReadNameStr.substr(dotLoc + 1, slashLoc - 1 - dotLoc - 1 + 1);
	return atoi(readNOstr.c_str());
}

bool end1orEnd2(int tmpFlag)
{
	bool end1_or_not = (tmpFlag & 0x40);
	bool end2_or_not = (tmpFlag & 0x80);
	if(end1_or_not && (!end2_or_not))
		return true;
	else if((!end1_or_not) && end2_or_not)
		return false;
	else
	{
		cout << "error in end1orEnd2" << endl;
		cout << "tmpFlag: " << tmpFlag << endl;
		cout << "end1_or_not: " << end1_or_not << endl;
		cout << "end2_or_not: " << end2_or_not << endl;
	}
}

int main(int argc, char**argv)
{
	if(argc != 4)
	{
		cout << "Executable readPairNum inputSam outputFolder" << endl;
		exit(1);
	}
	string readPairNumStr = argv[1];
	int readPairNum = atoi(readPairNumStr.c_str());
	int totalAlignmentNumber = readPairNum * 2;

	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	vector<int> alignmentNumVec_end1;
	vector<int> alignmentNumVec_end2;
	//vector<int> samVec_end1;
	//vector<int> samVec_end2;
	for(int tmp = 0; tmp < readPairNum; tmp++)
	{
		alignmentNumVec_end1.push_back(0);
		alignmentNumVec_end2.push_back(0);
	}
	cout << "start to count samNum for each read" << endl;
	log_ofs << "start to count samNum for each read" << endl;
	string inputSAMpath = argv[2];
	ifstream sam_ifs(inputSAMpath.c_str());
	int tmpLineNO = 0;
	while(!(sam_ifs.eof()))
	{
		string samStr;
		getline(sam_ifs, samStr);
		 if(sam_ifs.eof()||(samStr == ""))
		 	break;
		if(samStr.at(0) == '@')
			continue;
		//totalAlignNum ++;
		//cout << "samStr: " << samStr << endl;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 100000;
		if(tmpLineNO == tmpThousandIndex * 100000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec.push_back(samStr.substr(startLoc));
		string tmpReadNameStr = samFieldVec[0];
		string flagStr = samFieldVec[1];
		int tmpFlag = atoi(flagStr.c_str());
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(mappedOrNot_bool)
		{
			bool tmpRead_end1_or_end2_bool = end1orEnd2(tmpFlag);
			int tmpRead_NO = getReadNOfromReadName(tmpReadNameStr);
			if(tmpRead_NO > readPairNum)
			{
				cout << "invalid tmpRead_NO: " << tmpRead_NO << endl;
				exit(1);
			}
			int tmpRead_NO_index = tmpRead_NO - 1;
			if(tmpRead_end1_or_end2_bool)
				alignmentNumVec_end1[tmpRead_NO_index] ++;
			else
				alignmentNumVec_end2[tmpRead_NO_index] ++;
		}
	}
	sam_ifs.close();
	cout << "start to output uniqueMultiReadNum" << endl;
	log_ofs << "start to output uniqueMultiReadNum" << endl;
	int uniqueReadNum = 0;
	int multiReadNum = 0;
	int unmapReadNum = 0;
	for(int tmp = 0; tmp < readPairNum; tmp ++)
	{
		int tmpSamNum_end1 = alignmentNumVec_end1[tmp];
		if(tmpSamNum_end1 == 0)
			unmapReadNum ++;
		else if(tmpSamNum_end1 == 1)
			uniqueReadNum ++;
		else
			multiReadNum ++;
		int tmpSamNum_end2 = alignmentNumVec_end2[tmp];
		if(tmpSamNum_end2 == 0)
			unmapReadNum ++;
		else if(tmpSamNum_end2 == 1)
			uniqueReadNum ++;
		else
			multiReadNum ++;
	}
	double uniqueReadPerc = ((double)uniqueReadNum/(double)totalAlignmentNumber)*100;
	double multiReadPerc = ((double)multiReadNum/(double)totalAlignmentNumber)*100;
	double unmapReadPerc = ((double)unmapReadNum/(double)totalAlignmentNumber)*100;
	log_ofs << "Total read #:\t" << totalAlignmentNumber << endl;
	log_ofs << "Unique read #:\t" << uniqueReadNum << "\t" << uniqueReadPerc << "%" << endl;
	log_ofs << "Multi read #:\t" << multiReadNum << "\t" << multiReadPerc << "%" << endl;
	log_ofs << "Unmap read #:\t" << unmapReadNum << "\t" << unmapReadPerc << "%" << endl;
	log_ofs.close();

	string uniqueSamFile = outputFolderStr + "/unique.sam";
	string multiSamFile = outputFolderStr + "/multi.sam";
	string unmapSamFile = outputFolderStr + "/unmap.sam";
	ofstream uniqueSam_ofs(uniqueSamFile.c_str());
	ofstream multiSam_ofs(multiSamFile.c_str());
	ofstream unmapSam_ofs(unmapSamFile.c_str());
	ifstream sam_ifs_2(inputSAMpath.c_str());
	tmpLineNO = 0;
	while(!(sam_ifs_2.eof()))
	{
		string samStr;
		getline(sam_ifs_2, samStr);
		 if(sam_ifs_2.eof()||(samStr == ""))
		 	break;
		if(samStr.at(0) == '@')
			continue;
		//totalAlignNum ++;
		//cout << "samStr: " << samStr << endl;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 100000;
		if(tmpLineNO == tmpThousandIndex * 100000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec.push_back(samStr.substr(startLoc));
		string tmpReadNameStr = samFieldVec[0];
		string flagStr = samFieldVec[1];
		int tmpFlag = atoi(flagStr.c_str());
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(mappedOrNot_bool)
		{
			bool tmpRead_end1_or_end2_bool = end1orEnd2(tmpFlag);
			int tmpRead_NO = getReadNOfromReadName(tmpReadNameStr);
			if(tmpRead_NO > readPairNum)
			{
				cout << "invalid tmpRead_NO: " << tmpRead_NO << endl;
				exit(1);
			}
			int tmpRead_NO_index = tmpRead_NO - 1;
			int tmpAlignmentNum = 0;
			if(tmpRead_end1_or_end2_bool)
				tmpAlignmentNum = alignmentNumVec_end1[tmpRead_NO_index];
			else
				tmpAlignmentNum = alignmentNumVec_end2[tmpRead_NO_index];
			if(tmpAlignmentNum == 0)
			{
				cout << "error ! tmpAlignmentNum == 0, but mappedOrNot_bool = true" << endl;
				exit(1);
				unmapSam_ofs << samStr << endl;
			}
			else if(tmpAlignmentNum == 1)
				uniqueSam_ofs << samStr << endl;
			else
				multiSam_ofs << samStr << endl;
		}
		else
			unmapSam_ofs << samStr << endl;
	}
	sam_ifs_2.close();

	uniqueSam_ofs.close();
	multiSam_ofs.close();
	unmapSam_ofs.close();
	return  0;
}