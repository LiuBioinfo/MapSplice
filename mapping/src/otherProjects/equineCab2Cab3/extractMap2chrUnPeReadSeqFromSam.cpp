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

// bool primaryOrNot(int tmpFlag)
// {
// 	if(tmpFlag & 0x100)
// 		return false;
// 	else
// 		return true;
// }

string covertChar2ReverseComplement(const string& Ori_Char)
{
	if(Ori_Char == "A")
	{
		return "T";
	}
	else if(Ori_Char == "T")
	{
		return "A";
	}
	else if(Ori_Char == "G")
	{
		return "C";
	}
	else if(Ori_Char == "C")
	{
		return "G";
	}
	else if(Ori_Char == "N")
	{
		return "N";
	}
	else
	{
		cout << "incorrect Ori_Char in covertCharToReverseComplement" << endl;
		exit(1);
		return "X";
	}
}

string convertReadSeq2ReverseComplement(const string& originalString)
{
	int stringLength = originalString.size();
	string resultString = covertChar2ReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertChar2ReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

string convertQualSeq2Reverse(const string& originalQualityScoreString)
{
	int stringLength = originalQualityScoreString.size();
	string resultString = originalQualityScoreString.substr(stringLength-1, 1);//covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + originalQualityScoreString.substr(stringLength-1-tmp, 1);
			//covertCharToReverseComplement(originalString.substr(stringLength-1-tmp, 1));
	}
	return resultString;
}

bool for_or_rev(int tmpFlag)
{
	if(tmpFlag & 0x10)
		return false;
	else
		return true;
}

bool end1_or_end2(int tmpFlag)
{
	if(tmpFlag & 0x40)
		return true;
	else
		return false;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSam" << endl;
		cout << "#2 outputFilePrefix" << endl;
		exit(1);
	}

	string inputSam = argv[1];
	string outputFilePrefix = argv[2];
	string outputMap2chrUnSamFile = outputFilePrefix + ".map2chrUn.sam";
	string outputReadFile_end1 = outputFilePrefix + ".map2chrUn.1.fastq";
	string outputReadFile_end2 = outputFilePrefix + ".map2chrUn.2.fastq";
	ofstream map2chrUnSam_ofs(outputMap2chrUnSamFile.c_str());
	ofstream read_ofs_1(outputReadFile_end1.c_str());
	ofstream read_ofs_2(outputReadFile_end2.c_str());
	ifstream sam_ifs(inputSam.c_str());
	int tmpLineNO = 0;
	while(!sam_ifs.eof())
	{
		string tmpSamStr_nor;
		getline(sam_ifs, tmpSamStr_nor);
		if(tmpSamStr_nor == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 1000000;
		if(tmpLineNO == tmpThousandIndex * 1000000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		if(tmpSamStr_nor.at(0) == '@')
			continue;
		string tmpSamStr_rcm;
		getline(sam_ifs, tmpSamStr_rcm);
		
		vector<string> samFieldVec_nor;
		vector<string> samFieldVec_rcm;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpSamStr_nor.find("\t", startLoc);
			string tmpSamField = tmpSamStr_nor.substr(startLoc, tabLoc-startLoc);
			samFieldVec_nor.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpSamStr_rcm.find("\t", startLoc);
			string tmpSamField = tmpSamStr_rcm.substr(startLoc, tabLoc-startLoc);
			samFieldVec_rcm.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		if((samFieldVec_nor[2] == "chrUn")&&(samFieldVec_rcm[2] == "chrUn"))
		{
			string tmpFlagIntStr_nor = samFieldVec_nor[1];
			string tmpFlagIntStr_rcm = samFieldVec_rcm[1];
			int tmpFlagInt_nor = atoi(tmpFlagIntStr_nor.c_str());
			int tmpFlagInt_rcm = atoi(tmpFlagIntStr_rcm.c_str());
			bool primaryOrNot_bool_nor = primaryOrNot(tmpFlagInt_nor);
			bool primaryOrNot_bool_rcm = primaryOrNot(tmpFlagInt_rcm);	
			if(!(primaryOrNot_bool_nor && primaryOrNot_bool_rcm))
				continue;
			bool norRead_end1_or_end2_bool = end1_or_end2(tmpFlagInt_nor);
			bool rcmRead_end1_or_end2_bool = end1_or_end2(tmpFlagInt_rcm);
			string tmpReadId_nor = samFieldVec_nor[0];
			string tmpReadSeq_nor = samFieldVec_nor[9];
			string tmpReadQual_nor = samFieldVec_nor[10];
			string tmpReadId_rcm = samFieldVec_rcm[0];
			string tmpReadSeq_rcm = samFieldVec_rcm[9];
			string tmpReadQual_rcm = samFieldVec_rcm[10];

			string tmpId_end1, tmpId_end2, tmpSeq_end1, tmpSeq_end2, tmpQual_end1, tmpQual_end2;
			if(norRead_end1_or_end2_bool && (!rcmRead_end1_or_end2_bool)) // norRead--end1; rcmRead--end2
			{
				tmpId_end1 = "@" + tmpReadId_nor + "/1";
				tmpSeq_end1 = tmpReadSeq_nor;
				tmpQual_end1 = tmpReadQual_nor;
				tmpId_end2 = "@" + tmpReadId_rcm + "/2";
				tmpSeq_end2 = convertReadSeq2ReverseComplement(tmpReadSeq_rcm);
				tmpQual_end2 = convertQualSeq2Reverse(tmpReadQual_rcm);
			}
			else if((!norRead_end1_or_end2_bool) && rcmRead_end1_or_end2_bool) // norRead--end2; rcmRead--end1
			{
				tmpId_end1 = "@" + tmpReadId_rcm + "/1";
				tmpSeq_end1 = convertReadSeq2ReverseComplement(tmpReadSeq_rcm);
				tmpQual_end1 = convertQualSeq2Reverse(tmpReadQual_rcm);
				tmpId_end2 = "@" + tmpReadId_nor + "/2";
				tmpSeq_end2 = tmpReadSeq_nor;
				tmpQual_end2 = tmpReadQual_nor;
			}
			else
			{
				cout << "ERROR!" << endl;
				cout << "norRead_end1_or_end2_bool: " << norRead_end1_or_end2_bool << endl;
				cout << "rcmRead_end1_or_end2_bool: " << rcmRead_end1_or_end2_bool << endl;
				exit(1);
			}
			map2chrUnSam_ofs << tmpSamStr_nor << endl << tmpSamStr_rcm << endl;
			read_ofs_1 << tmpId_end1 << endl << tmpSeq_end1 << endl << "+" << endl << tmpQual_end1 << endl;
			read_ofs_2 << tmpId_end2 << endl << tmpSeq_end2 << endl << "+" << endl << tmpQual_end2 << endl;
		}
	}
	map2chrUnSam_ofs.close();
	sam_ifs.close();
	read_ofs_2.close();
	read_ofs_1.close();
	return 0;
}