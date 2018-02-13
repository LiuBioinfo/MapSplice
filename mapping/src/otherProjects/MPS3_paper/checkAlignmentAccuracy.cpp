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

bool paired(int tmpFlag)
{
	if(tmpFlag & 0x2)
		return true;
	else
		return false;
}

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

bool correctlyMapped(vector<int>& gt_chrNameInt_vec, vector<int>& gt_startPos_vec, vector<int>& gt_endPos_vec, 
	int tmpRead_NO, int mapChrNameInt, int mapChrPos, int mapChrEndPos)
{
	int tmpGT_chrName = gt_chrNameInt_vec[tmpRead_NO];
	int tmpGT_startPos = gt_startPos_vec[tmpRead_NO];
	int tmpGT_endPos = gt_endPos_vec[tmpRead_NO];
	//cout << "tmpGT_chrName: " << tmpGT_chrName << endl;
	//cout << "tmpGT_startPos: " << tmpGT_startPos << endl;
	//cout << "tmpGT_endPos: " << tmpGT_endPos << endl;
	#ifdef PERFECT
	if((mapChrNameInt == tmpGT_chrName)&&(tmpGT_startPos == mapChrPos)&&(tmpGT_endPos == mapChrEndPos))
		return true;
	else
		return false;
	#endif
	#ifdef WITHIN
	if((mapChrNameInt == tmpGT_chrName)&&(tmpGT_startPos <= mapChrPos)&&(tmpGT_endPos >= mapChrEndPos))
		return true;
	else
		return false;
	#endif

	if((mapChrNameInt != tmpGT_chrName)||(tmpGT_startPos > mapChrEndPos)||(tmpGT_endPos < mapChrPos))
		return false;
	else
		return true;
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

int getSJnumFromCigarString(string& tmpCigarStringStr)
{
	vector<Jump_Code> tmpJumpCodeVec;
	//cout << "cigarString " << cigarString << endl;
	cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);	
	int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
	int tmpSJnum = 0;
	for(int tmpJumpCodeIndex = 0; tmpJumpCodeIndex < tmpJumpCodeVecSize; tmpJumpCodeIndex++)
	{
		string tmpJumpCodeType = tmpJumpCodeVec[tmpJumpCodeIndex].type;
		if(tmpJumpCodeType == "N")
			tmpSJnum ++;
	}
	return tmpSJnum;
}

bool parseBeersSam(string& tmpStr, int& tmpRead_NO, bool& tmpRead_end1_or_end2_bool, 
	int& tmpRead_chrNameInt, int& tmpRead_startPos, int& tmpRead_endPos, Index_Info* indexInfo)
{
	vector<string> tmpFieldStrVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		// if(tabLoc == string::npos)
		// {	
		// 	tmpFieldStrVec.push_back(tmpAluStr.substr(startLoc));
		// 	break;
		// }
		tmpFieldStrVec.push_back(tmpStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}	
	string tmpReadNameStr = tmpFieldStrVec[0];
	int dotLocInReadName = tmpReadNameStr.find(".");
	int slashLocInReadName = tmpReadNameStr.find("/");
	string tmpReadNOstr = tmpReadNameStr.substr(dotLocInReadName + 1, slashLocInReadName - dotLocInReadName - 1);
	//cout << "tmpReadNOstr: " << tmpReadNOstr << endl;
	string tmpReadEnd1OrEnd2Str = tmpReadNameStr.substr(slashLocInReadName + 1, 1);
	//cout << "tmpReadEnd1OrEnd2Str: " << tmpReadEnd1OrEnd2Str << endl;
	tmpRead_NO = atoi(tmpReadNOstr.c_str());
	//cout << "tmpRead_NO: " << tmpRead_NO << endl;
	if(tmpReadEnd1OrEnd2Str == "1")
		tmpRead_end1_or_end2_bool = true;
	else if(tmpReadEnd1OrEnd2Str == "2")
		tmpRead_end1_or_end2_bool = false;
	else
	{
		cout << "incorrect tmpReadEnd1OrEnd2Str: " << tmpReadEnd1OrEnd2Str << endl;
		exit(1);
	}
	string tmpChrNameStr = tmpFieldStrVec[2];
	//cout << "tmpChrNameStr: " << tmpChrNameStr << endl;
	tmpRead_chrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
	//cout << "tmpRead_chrNameInt: " << tmpRead_chrNameInt << endl;
	if(tmpRead_chrNameInt < 0)
		return false;
	string tmpChrPosStr = tmpFieldStrVec[3];
	//cout << "tmpChrPosStr: " << tmpChrPosStr << endl;
	tmpRead_startPos = atoi(tmpChrPosStr.c_str());
	string tmpCigarString = tmpFieldStrVec[5];
	//cout << "tmpCigarString: " << tmpCigarString << endl;
	tmpRead_endPos = getEndMapPos(tmpRead_startPos, tmpCigarString);
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputGroundTruthSam input2cmpSam outputFolderPath readPairNum" << endl;
		exit(1);
	}
	string inputIndexFolderPath = argv[1];
	string input2cmpSam = argv[3];
	string readPairNumStr = argv[5];
	int readPairNum = atoi(readPairNumStr.c_str());
	int totalAlignmentNumber = readPairNum * 2;

	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	log_ofs << "start to initiate indexInfo" << endl;
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;	

	cout << "start to initiate gtInfoVec" << endl;
	log_ofs << "start to initiate gtInfoVec" << endl;		
	vector<int> gt_chrNameInt_end1;
	vector<int> gt_chrNameInt_end2;
	vector<int> gt_startPos_end1;
	vector<int> gt_startPos_end2;
	vector<int> gt_endPos_end1;
	vector<int> gt_endPos_end2;	
	for(int tmp = 0; tmp < readPairNum; tmp++)
	{
		gt_chrNameInt_end1.push_back(0);
		gt_chrNameInt_end2.push_back(0);
		gt_startPos_end1.push_back(0);
		gt_startPos_end2.push_back(0);
		gt_endPos_end1.push_back(0);
		gt_endPos_end2.push_back(0);
	}
	cout << "start to read ground truth SAM file and load info 2 gtInfoVec" << endl;
	log_ofs << "start to read ground truth SAM file and load info 2 gtInfoVec" << endl;
	int tmpLineNO = 0;
	string inputGroundTruthSam = argv[2];	
	ifstream gt_ifs(inputGroundTruthSam.c_str());
	while(!gt_ifs.eof())
	{
		string tmpStr;
		getline(gt_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 100000;
		if(tmpLineNO == tmpThousandIndex * 100000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		int tmpRead_NO;
		bool tmpRead_end1_or_end2_bool;
		int tmpRead_chrNameInt;
		int tmpRead_startPos;
		int tmpRead_endPos;
		bool tmpRead_parse_bool = parseBeersSam(tmpStr, tmpRead_NO, tmpRead_end1_or_end2_bool, 
			tmpRead_chrNameInt, tmpRead_startPos, tmpRead_endPos, indexInfo);
		//cout << "tmpRead_parse_bool: " << tmpRead_parse_bool << endl;
		if(!tmpRead_parse_bool)
			continue;
		if((tmpRead_NO < 1)||(tmpRead_NO > readPairNum))
		{
			cout << "tmpRead_NO error: " << tmpRead_NO << endl;
			cout << "tmpRead_sam: " << tmpStr << endl;
			exit(1);
		}
		int tmpRead_NO_index = tmpRead_NO - 1;
		//cout << "tmpRead_NO_index: " << tmpRead_NO_index << endl;
		if(tmpRead_end1_or_end2_bool)
		{
			gt_chrNameInt_end1[tmpRead_NO_index] = tmpRead_chrNameInt;
			gt_startPos_end1[tmpRead_NO_index] = tmpRead_startPos;
			gt_endPos_end1[tmpRead_NO_index] = tmpRead_endPos;
		}
		else
		{
			gt_chrNameInt_end2[tmpRead_NO_index] = tmpRead_chrNameInt;
			gt_startPos_end2[tmpRead_NO_index] = tmpRead_startPos;
			gt_endPos_end2[tmpRead_NO_index] = tmpRead_endPos;			
		}
	}
	gt_ifs.close();

	cout << "start to read intput2cmp SAM file and compare2gt" << endl;
	log_ofs << "start to read intput2cmp SAM file and compare2gt" << endl;
	int unmapped_alignmentNum_inFile = 0;
	int notPrimary_alignmentNum = 0;

	int primary_paired_total_alignmentNum = 0;
	int primary_paired_unspliced_alignmentNum = 0;
	int primary_paired_singleSpliced_alignmentNum = 0;
	int primary_paired_multiSpliced_alignmentNum = 0;
	int primary_unpaired_total_alignmentNum = 0;
	int primary_unpaired_unspliced_alignmentNum = 0;
	int primary_unpaired_singleSpliced_alignmentNum = 0;
	int primary_unpaired_multiSpliced_alignmentNum = 0;

	int primary_paired_total_alignmentNum_correct = 0;
	int primary_paired_unspliced_alignmentNum_correct = 0;
	int primary_paired_singleSpliced_alignmentNum_correct = 0;
	int primary_paired_multiSpliced_alignmentNum_correct = 0;
	int primary_unpaired_total_alignmentNum_correct = 0;
	int primary_unpaired_unspliced_alignmentNum_correct = 0;
	int primary_unpaired_singleSpliced_alignmentNum_correct = 0;
	int primary_unpaired_multiSpliced_alignmentNum_correct = 0;

	int primary_paired_total_alignmentNum_incorrect = 0;
	int primary_paired_unspliced_alignmentNum_incorrect = 0;
	int primary_paired_singleSpliced_alignmentNum_incorrect = 0;
	int primary_paired_multiSpliced_alignmentNum_incorrect = 0;
	int primary_unpaired_total_alignmentNum_incorrect = 0;
	int primary_unpaired_unspliced_alignmentNum_incorrect = 0;
	int primary_unpaired_singleSpliced_alignmentNum_incorrect = 0;
	int primary_unpaired_multiSpliced_alignmentNum_incorrect = 0;

	string inputSAMpath = argv[3];
	ifstream sam_ifs(inputSAMpath.c_str());
	tmpLineNO = 0;
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
		for(int tmp = 0; tmp < 6; tmp++)
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
		string mapChrNameStr = samFieldVec[2];
		int mapChrNameInt = indexInfo->convertStringToInt(mapChrNameStr);
		string mapChrPosStr = samFieldVec[3];
		int mapChrPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];			
	
		//cout << "tmpFlag" << tmpFlag << endl;
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		bool primaryOrNot_bool = primaryOrNot(tmpFlag);
		bool paired_bool = paired(tmpFlag);
		if(!mappedOrNot_bool)
			unmapped_alignmentNum_inFile ++;
		else
		{
			if(primaryOrNot_bool)
			{		
				int tmpRead_NO = getReadNOfromReadName(tmpReadNameStr);
				int tmpRead_NO_index = tmpRead_NO - 1;
				//cout << "tmpRead_NO: " << tmpRead_NO << endl;
				bool tmpRead_end1_or_end2_bool = end1orEnd2(tmpFlag);
				//cout << "tmpRead_end1_or_end2_bool: " << tmpRead_end1_or_end2_bool << endl;
				int tmpRead_endPos = getEndMapPos(mapChrPos, cigarString);
				//cout << "tmpRead_chrNameInt: " << mapChrNameInt << endl;
				//cout << "tmpRead_startPos: " << mapChrPos << endl;
				//cout << "tmpRead_endPos: " << tmpRead_endPos << endl;
				bool correctlyMappedOrNot_bool; 
				if(tmpRead_end1_or_end2_bool) 
					correctlyMappedOrNot_bool = correctlyMapped(gt_chrNameInt_end1, gt_startPos_end1, 
						gt_endPos_end1, tmpRead_NO_index, mapChrNameInt, mapChrPos, tmpRead_endPos);
				else
					correctlyMappedOrNot_bool = correctlyMapped(gt_chrNameInt_end2, gt_startPos_end2, 
						gt_endPos_end2, tmpRead_NO_index, mapChrNameInt, mapChrPos, tmpRead_endPos);
				//cout << "correctlyMappedOrNot_bool: " << correctlyMappedOrNot_bool << endl;
				int tmpAlignmentSJnum = getSJnumFromCigarString(cigarString);
				if(paired_bool)
				{
					primary_paired_total_alignmentNum ++;
					if(tmpAlignmentSJnum == 0)
						primary_paired_unspliced_alignmentNum ++;
					else if(tmpAlignmentSJnum == 1)
						primary_paired_singleSpliced_alignmentNum ++;
					else
						primary_paired_multiSpliced_alignmentNum ++;
					if(correctlyMappedOrNot_bool)
					{
						primary_paired_total_alignmentNum_correct ++;
						if(tmpAlignmentSJnum == 0)
							primary_paired_unspliced_alignmentNum_correct ++;
						else if(tmpAlignmentSJnum == 1)
							primary_paired_singleSpliced_alignmentNum_correct ++;
						else
							primary_paired_multiSpliced_alignmentNum_correct ++;						
					}
					else
					{
						primary_paired_total_alignmentNum_incorrect ++;
						if(tmpAlignmentSJnum == 0)
							primary_paired_unspliced_alignmentNum_incorrect ++;
						else if(tmpAlignmentSJnum == 1)
							primary_paired_singleSpliced_alignmentNum_incorrect ++;
						else
							primary_paired_multiSpliced_alignmentNum_incorrect ++;						
					}
				}
				else
				{
					primary_unpaired_total_alignmentNum ++;
					if(tmpAlignmentSJnum == 0)
						primary_unpaired_unspliced_alignmentNum ++;
					else if(tmpAlignmentSJnum == 1)
						primary_unpaired_singleSpliced_alignmentNum ++;
					else
						primary_unpaired_multiSpliced_alignmentNum ++;
					if(correctlyMappedOrNot_bool)
					{
						primary_unpaired_total_alignmentNum_correct ++;
						if(tmpAlignmentSJnum == 0)
							primary_unpaired_unspliced_alignmentNum_correct ++;
						else if(tmpAlignmentSJnum == 1)
							primary_unpaired_singleSpliced_alignmentNum_correct ++;
						else
							primary_unpaired_multiSpliced_alignmentNum_correct ++;						
					}
					else
					{
						primary_unpaired_total_alignmentNum_incorrect ++;
						if(tmpAlignmentSJnum == 0)
							primary_unpaired_unspliced_alignmentNum_incorrect ++;
						else if(tmpAlignmentSJnum == 1)
							primary_unpaired_singleSpliced_alignmentNum_incorrect ++;
						else
							primary_unpaired_multiSpliced_alignmentNum_incorrect ++;							
					}
				}			
			}
			else
				notPrimary_alignmentNum ++;
		}
	}	

	double unmapped_alignmentPerc_inFile = ((double)unmapped_alignmentNum_inFile/(double)totalAlignmentNumber) * 100;

	double primary_paired_total_alignmentPerc = ((double)primary_paired_total_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_paired_unspliced_alignmentPerc = ((double)primary_paired_unspliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_paired_singleSpliced_alignmentPerc = ((double)primary_paired_singleSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_paired_multiSpliced_alignmentPerc = ((double)primary_paired_multiSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_total_alignmentPerc = ((double)primary_unpaired_total_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_unspliced_alignmentPerc = ((double)primary_unpaired_unspliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_singleSpliced_alignmentPerc = ((double)primary_unpaired_singleSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_multiSpliced_alignmentPerc = ((double)primary_unpaired_multiSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;

	double primary_paired_total_alignmentPerc_correct = ((double)primary_paired_total_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_paired_unspliced_alignmentPerc_correct = ((double)primary_paired_unspliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_paired_singleSpliced_alignmentPerc_correct = ((double)primary_paired_singleSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_paired_multiSpliced_alignmentPerc_correct = ((double)primary_paired_multiSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_total_alignmentPerc_correct = ((double)primary_unpaired_total_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_unspliced_alignmentPerc_correct = ((double)primary_unpaired_unspliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_singleSpliced_alignmentPerc_correct = ((double)primary_unpaired_singleSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_multiSpliced_alignmentPerc_correct = ((double)primary_unpaired_multiSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;

	double primary_paired_total_alignmentPerc_incorrect = ((double)primary_paired_total_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_paired_unspliced_alignmentPerc_incorrect = ((double)primary_paired_unspliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_paired_singleSpliced_alignmentPerc_incorrect = ((double)primary_paired_singleSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_paired_multiSpliced_alignmentPerc_incorrect = ((double)primary_paired_multiSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_total_alignmentPerc_incorrect = ((double)primary_unpaired_total_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_unspliced_alignmentPerc_incorrect = ((double)primary_unpaired_unspliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_singleSpliced_alignmentPerc_incorrect = ((double)primary_unpaired_singleSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_unpaired_multiSpliced_alignmentPerc_incorrect = ((double)primary_unpaired_multiSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;	


	int primary_total_alignmentNum = primary_paired_total_alignmentNum + primary_unpaired_total_alignmentNum;
	int primary_unspliced_alignmentNum = primary_paired_unspliced_alignmentNum + primary_unpaired_unspliced_alignmentNum;
	int primary_singleSpliced_alignmentNum = primary_paired_singleSpliced_alignmentNum + primary_unpaired_singleSpliced_alignmentNum;
	int primary_multiSpliced_alignmentNum = primary_paired_multiSpliced_alignmentNum + primary_unpaired_multiSpliced_alignmentNum;	
	double primary_total_alignmentPerc = ((double)primary_total_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_unspliced_alignmentPerc = ((double)primary_unspliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_singleSpliced_alignmentPerc = ((double)primary_singleSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;
	double primary_multiSpliced_alignmentPerc = ((double)primary_multiSpliced_alignmentNum/(double)totalAlignmentNumber) * 100;

	int primary_total_alignmentNum_correct = primary_paired_total_alignmentNum_correct + primary_unpaired_total_alignmentNum_correct;
	int primary_unspliced_alignmentNum_correct = primary_paired_unspliced_alignmentNum_correct + primary_unpaired_unspliced_alignmentNum_correct;
	int primary_singleSpliced_alignmentNum_correct = primary_paired_singleSpliced_alignmentNum_correct + primary_unpaired_singleSpliced_alignmentNum_correct;
	int primary_multiSpliced_alignmentNum_correct = primary_paired_multiSpliced_alignmentNum_correct + primary_unpaired_multiSpliced_alignmentNum_correct;	
	double primary_total_alignmentPerc_correct = ((double)primary_total_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_unspliced_alignmentPerc_correct = ((double)primary_unspliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_singleSpliced_alignmentPerc_correct = ((double)primary_singleSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;
	double primary_multiSpliced_alignmentPerc_correct = ((double)primary_multiSpliced_alignmentNum_correct/(double)totalAlignmentNumber) * 100;

	int primary_total_alignmentNum_incorrect = primary_paired_total_alignmentNum_incorrect + primary_unpaired_total_alignmentNum_incorrect;
	int primary_unspliced_alignmentNum_incorrect = primary_paired_unspliced_alignmentNum_incorrect + primary_unpaired_unspliced_alignmentNum_incorrect;
	int primary_singleSpliced_alignmentNum_incorrect = primary_paired_singleSpliced_alignmentNum_incorrect + primary_unpaired_singleSpliced_alignmentNum_incorrect;
	int primary_multiSpliced_alignmentNum_incorrect = primary_paired_multiSpliced_alignmentNum_incorrect + primary_unpaired_multiSpliced_alignmentNum_incorrect;	
	double primary_total_alignmentPerc_incorrect = ((double)primary_total_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_unspliced_alignmentPerc_incorrect = ((double)primary_unspliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_singleSpliced_alignmentPerc_incorrect = ((double)primary_singleSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;
	double primary_multiSpliced_alignmentPerc_incorrect = ((double)primary_multiSpliced_alignmentNum_incorrect/(double)totalAlignmentNumber) * 100;	
	//double twoEndsNotProperlyMapped_alignmentPerc = ((double)twoEndsNotProperlyMapped_alignmentNum/(double)totalAlignmentNumber) * 100;

	double notPrimary_alignmentPerc = ((double)notPrimary_alignmentNum/(double)totalAlignmentNumber) * 100;
	
	int unmappedAlignmentNum_real = totalAlignmentNumber - primary_total_alignmentNum;
	double unmappedAlignmentPerc_real = ((double)unmappedAlignmentNum_real/(double)totalAlignmentNumber) * 100;

	string alignmentTypeNumPercFile = outputFolderStr + "alignmentTypeNumPerc.txt";
	ofstream alignmentTypeNumPerc_ofs(alignmentTypeNumPercFile.c_str());

	alignmentTypeNumPerc_ofs << "Total_readNum:\t" << totalAlignmentNumber << endl;
	alignmentTypeNumPerc_ofs << "primary_total_alignmentNumPerc:\t" << primary_total_alignmentNum << "\t" << primary_total_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_unspliced_alignmentNumPerc:\t" << primary_unspliced_alignmentNum << "\t" << primary_unspliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_singleSpliced_alignmentNumPerc:\t" << primary_singleSpliced_alignmentNum << "\t" << primary_singleSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_multiSpliced_alignmentNumPerc:\t" << primary_multiSpliced_alignmentNum << "\t" << primary_multiSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_total_alignmentNumPerc:\t" << primary_paired_total_alignmentNum << "\t" << primary_paired_total_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_unspliced_alignmentNumPerc:\t" << primary_paired_unspliced_alignmentNum << "\t" << primary_paired_unspliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_singleSpliced_alignmentNumPerc:\t" << primary_paired_singleSpliced_alignmentNum << "\t" << primary_paired_singleSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_multiSpliced_alignmentNumPerc:\t" << primary_paired_multiSpliced_alignmentNum << "\t" << primary_paired_multiSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_total_alignmentNumPerc:\t" << primary_unpaired_total_alignmentNum << "\t" << primary_unpaired_total_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_unspliced_alignmentNumPerc:\t" << primary_unpaired_unspliced_alignmentNum << "\t" << primary_unpaired_unspliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_singleSpliced_alignmentNumPerc:\t" << primary_unpaired_singleSpliced_alignmentNum << "\t" << primary_unpaired_singleSpliced_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_multiSpliced_alignmentNumPerc:\t" << primary_unpaired_multiSpliced_alignmentNum << "\t" << primary_unpaired_multiSpliced_alignmentPerc << endl;
	
	alignmentTypeNumPerc_ofs << endl << "#####################################" << endl << endl;
	alignmentTypeNumPerc_ofs << "primary_total_alignmentNumPerc_correct:\t" << primary_total_alignmentNum_correct << "\t" << primary_total_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_unspliced_alignmentNumPerc_correct:\t" << primary_unspliced_alignmentNum_correct << "\t" << primary_unspliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_singleSpliced_alignmentNumPerc_correct:\t" << primary_singleSpliced_alignmentNum_correct << "\t" << primary_singleSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_multiSpliced_alignmentNumPerc_correct:\t" << primary_multiSpliced_alignmentNum_correct << "\t" << primary_multiSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_total_alignmentNumPerc_correct:\t" << primary_paired_total_alignmentNum_correct << "\t" << primary_paired_total_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_unspliced_alignmentNumPerc_correct:\t" << primary_paired_unspliced_alignmentNum_correct << "\t" << primary_paired_unspliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_singleSpliced_alignmentNumPerc_correct:\t" << primary_paired_singleSpliced_alignmentNum_correct << "\t" << primary_paired_singleSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_multiSpliced_alignmentNumPerc_correct:\t" << primary_paired_multiSpliced_alignmentNum_correct << "\t" << primary_paired_multiSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_total_alignmentNumPerc_correct:\t" << primary_unpaired_total_alignmentNum_correct << "\t" << primary_unpaired_total_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_unspliced_alignmentNumPerc_correct:\t" << primary_unpaired_unspliced_alignmentNum_correct << "\t" << primary_unpaired_unspliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_singleSpliced_alignmentNumPerc_correct:\t" << primary_unpaired_singleSpliced_alignmentNum_correct << "\t" << primary_unpaired_singleSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_multiSpliced_alignmentNumPerc_correct:\t" << primary_unpaired_multiSpliced_alignmentNum_correct << "\t" << primary_unpaired_multiSpliced_alignmentPerc_correct << endl;
	alignmentTypeNumPerc_ofs << endl;

	alignmentTypeNumPerc_ofs << endl << "#####################################" << endl << endl;
	alignmentTypeNumPerc_ofs << "primary_total_alignmentNumPerc_incorrect:\t" << primary_total_alignmentNum_incorrect << "\t" << primary_total_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_unspliced_alignmentNumPerc_incorrect:\t" << primary_unspliced_alignmentNum_incorrect << "\t" << primary_unspliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_singleSpliced_alignmentNumPerc_incorrect:\t" << primary_singleSpliced_alignmentNum_incorrect << "\t" << primary_singleSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_multiSpliced_alignmentNumPerc_incorrect:\t" << primary_multiSpliced_alignmentNum_incorrect << "\t" << primary_multiSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_total_alignmentNumPerc_incorrect:\t" << primary_paired_total_alignmentNum_incorrect << "\t" << primary_paired_total_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_unspliced_alignmentNumPerc_incorrect:\t" << primary_paired_unspliced_alignmentNum_incorrect << "\t" << primary_paired_unspliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_singleSpliced_alignmentNumPerc_incorrect:\t" << primary_paired_singleSpliced_alignmentNum_incorrect << "\t" << primary_paired_singleSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_paired_multiSpliced_alignmentNumPerc_incorrect:\t" << primary_paired_multiSpliced_alignmentNum_incorrect << "\t" << primary_paired_multiSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_total_alignmentNumPerc_incorrect:\t" << primary_unpaired_total_alignmentNum_incorrect << "\t" << primary_unpaired_total_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_unspliced_alignmentNumPerc_incorrect:\t" << primary_unpaired_unspliced_alignmentNum_incorrect << "\t" << primary_unpaired_unspliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_singleSpliced_alignmentNumPerc_incorrect:\t" << primary_unpaired_singleSpliced_alignmentNum_incorrect << "\t" << primary_unpaired_singleSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << "primary_unpaired_multiSpliced_alignmentNumPerc_incorrect:\t" << primary_unpaired_multiSpliced_alignmentNum_incorrect << "\t" << primary_unpaired_multiSpliced_alignmentPerc_incorrect << endl;
	alignmentTypeNumPerc_ofs << endl;

	alignmentTypeNumPerc_ofs << endl << "#####################################" << endl << endl;
	alignmentTypeNumPerc_ofs << "notPrimary_alignmentPerc:\t" << notPrimary_alignmentNum << "\t" << notPrimary_alignmentPerc << endl;
	alignmentTypeNumPerc_ofs << endl;
	alignmentTypeNumPerc_ofs << "unmappedAlignmentPerc_real:\t" << unmappedAlignmentNum_real << "\t" << unmappedAlignmentPerc_real << endl;
	alignmentTypeNumPerc_ofs << "unmapped_alignmentPerc_inFile:\t" << unmapped_alignmentNum_inFile << "\t" << unmapped_alignmentPerc_inFile << endl;

	alignmentTypeNumPerc_ofs.close();
	log_ofs.close();
	free(chrom);
	delete indexInfo;
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}