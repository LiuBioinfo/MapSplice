// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNINFER_INFO_H
#define ALIGNINFER_INFO_H

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

using namespace std;

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

bool mapped_or_not(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
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

// bool bothEndsMapped(int tmpFlag)
// {
// 	if(tmpFlag & 0x2)
// 		return true;
// 	else
// 		return false;
// }

bool parseSamStr2readInfo(string& tmpSamStr, bool& end1orEnd2_bool, string& tmpReadName, string& tmpReadSeq, string& tmpQualSeq)
{
	vector<string> samFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 11; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpMafStr find ......" << endl;
			exit(1);
		}
		string tmpSamField = tmpSamStr.substr(startLoc, tabLoc-startLoc);
		samFieldVec.push_back(tmpSamField);
		startLoc = tabLoc + 1;
	}	
	tmpReadName = samFieldVec[0];
	string tmpFlagStr = samFieldVec[1];
	int tmpFlagInt = atoi(tmpFlagStr.c_str());
	end1orEnd2_bool = end1_or_end2(tmpFlagInt);
	bool tmpSam_for_or_rcm_bool = for_or_rev(tmpFlagInt);
	if(tmpSam_for_or_rcm_bool)
	{
		tmpReadSeq = samFieldVec[9];
		tmpQualSeq = samFieldVec[10];
	}
	else
	{
		string tmpReadSeq_inSam = samFieldVec[9];
		string tmpQualSeq_inSam = samFieldVec[10];
		tmpReadSeq = convertReadSeq2ReverseComplement(tmpReadSeq_inSam);
		tmpQualSeq = convertQualSeq2Reverse(tmpQualSeq_inSam);
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

bool parseSamStr2samInfo(string& tmpSamStr, bool& end1orEnd2_bool, bool& forOrRcm_bool, 
	int& tmpChrNameInt, int& tmpStartPos, vector<Jump_Code>& tmpJumpCodeVec, Index_Info* indexInfo)
{
	vector<string> samFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpMafStr find ......" << endl;
			exit(1);
		}
		string tmpSamField = tmpSamStr.substr(startLoc, tabLoc-startLoc);
		samFieldVec.push_back(tmpSamField);
		startLoc = tabLoc + 1;
	}	
	string tmpFlagStr = samFieldVec[1];
	int tmpFlagInt = atoi(tmpFlagStr.c_str());
	bool mappedOrNot_bool = mapped_or_not(tmpFlagInt);
	if(mappedOrNot_bool)
	{	
		end1orEnd2_bool = end1_or_end2(tmpFlagInt);
		forOrRcm_bool = for_or_rev(tmpFlagInt);
		string tmpChrNameStr = samFieldVec[2];
		string tmpChrPosStr = samFieldVec[3];
		string tmpCigarStringStr = samFieldVec[5];
		tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		if(tmpChrNameInt < 0)
			return false;
		tmpStartPos = atoi(tmpChrPosStr.c_str());
		cigarString2jumpCodeVec(tmpCigarStringStr, tmpJumpCodeVec);
		return true;
	}
	else
		return false;
}

class Xenograft_Read_Info
{
private:	
	bool SE_or_PE_bool;

	string readName_1;
	string readSeq_1;
	string readSeq_rcm_1;
	string qualSeq_1;
	string qualSeq_rev_1;

	string readName_2;
	string readSeq_2;
	string readSeq_rcm_2;
	string qualSeq_2;
	string qualSeq_rev_2;

public:
	Xenograft_Read_Info()
	{}

	int return_readLength_1()
	{
		return readSeq_1.length();
	}

	int return_readLength_2()
	{
		return readSeq_2.length();
	}

	bool initiate_readInfo_SE(string& samStr_SE)
	{
		SE_or_PE_bool = true;
		bool tmp_end1_or_end2_bool;
		parseSamStr2readInfo(samStr_SE, tmp_end1_or_end2_bool, readName_1, readSeq_1, qualSeq_1);
		readSeq_rcm_1 = convertReadSeq2ReverseComplement(readSeq_1);
		qualSeq_rev_1 = convertQualSeq2Reverse(qualSeq_1);
	}

	bool initiate_readInfo_PE(vector<string>& samStrVec)
	{
		SE_or_PE_bool = false;
		int index_sam_1stEnd1 = -1;
		int index_sam_1stEnd2 = -1;
		get1stSamIndex(samStrVec, index_sam_1stEnd1, index_sam_1stEnd2);
		string samStr_1stEnd1 = samStrVec[index_sam_1stEnd1];
		string samStr_1stEnd2 = samStrVec[index_sam_1stEnd2];
		if((samStr_1stEnd1 < 0)||(samStr_1stEnd2 < 0))
			return false;
		bool tmp_end1_or_end2_bool;
		parseSamStr2readInfo(samStr_1stEnd1, tmp_end1_or_end2_bool, readName_1, readSeq_1, qualSeq_1);
		parseSamStr2readInfo(samStr_1stEnd2, tmp_end1_or_end2_bool, readName_2, readSeq_2, qualSeq_2);		
		readSeq_rcm_1 = convertReadSeq2ReverseComplement(readSeq_1);
		qualSeq_rev_1 = convertQualSeq2Reverse(qualSeq_1);
		readSeq_rcm_2 = convertReadSeq2ReverseComplement(readSeq_2);
		qualSeq_rev_2 = convertQualSeq2Reverse(qualSeq_2);		
	}
};

class Xenograft_Sam_Info
{
private:
	bool SE_or_PE_bool;
	bool graft_or_host_bool;

	vector<int> chrNameIntVec_1;
	vector<int> startPosVec_1;
	vector< vector<Jump_Code> > jumpCodeVecVec_1;
	vector<bool> forOrRcmBoolVec_1; // true--for; false--rcm

	vector<int> chrNameIntVec_2;
	vector<int> startPosVec_2;
	vector< vector<Jump_Code> > jumpCodeVecVec_2;
	vector<bool> forOrRcmBoolVec_2; // true--for; false--rcm
public:
	Xenograft_Sam_Info()
	{}

	int return_alignmentNum_1()
	{
		return chrNameIntVec_1.size();
	}
	int return_alignmentNum_2()
	{
		return chrNameIntVec_2.size();
	}

	bool initiate_samInfo_SE(vector<string>& samStrVec, Index_Info* indexInfo)
	{
		int samStrVecSize = samStrVec.size();
		for(int tmp = 0; tmp < samStrVecSize; tmp++)
		{
			string tmpSamStr = samStrVec[tmp];
			bool tmp_end1OrEnd2_bool, tmp_forOrRcm_bool;
			int tmpChrNameInt, tmpStartPos;
			vector<Jump_Code> tmpJumpCodeVec;
			bool tmp_initiateSam_bool = parseSamStr2samInfo(tmpSamStr, tmp_end1OrEnd2_bool,
				tmp_forOrRcm_bool, tmpChrNameInt, tmpStartPos, tmpJumpCodeVec, indexInfo);
			if(!tmp_end1OrEnd2_bool)
			{
				cout << "error ! end2 sam found ! not valid for SE" << endl;
				exit(1); 
			}
			if(tmp_initiateSam_bool)
			{
				chrNameIntVec_1.push_back(tmpChrNameInt);
				startPosVec_1.push_back(tmpStartPos);
				jumpCodeVecVec_1.push_back(tmpJumpCodeVec);
				forOrRcmBoolVec_1.push_back(tmp_forOrRcm_bool);
			}
		}
	}

	bool initiate_samInfo_PE(vector<string>& samStrVec, Index_Info* indexInfo)
	{
		int samStrVecSize = samStrVec.size();
		for(int tmp = 0; tmp < samStrVecSize; tmp++)
		{
			string tmpSamStr = samStrVec[tmp];
			bool tmp_end1OrEnd2_bool, tmp_forOrRcm_bool;
			int tmpChrNameInt, tmpStartPos;
			vector<Jump_Code> tmpJumpCodeVec;
			bool tmp_initiateSam_bool = parseSamStr2samInfo(tmpSamStr, tmp_end1OrEnd2_bool,
				tmp_forOrRcm_bool, tmpChrNameInt, tmpStartPos, tmpJumpCodeVec, indexInfo);
			if(tmp_initiateSam_bool)
			{
				if(tmp_end1OrEnd2_bool)
				{	
					chrNameIntVec_1.push_back(tmpChrNameInt);
					startPosVec_1.push_back(tmpStartPos);
					jumpCodeVecVec_1.push_back(tmpJumpCodeVec);
					forOrRcmBoolVec_1.push_back(tmp_forOrRcm_bool);
				}
				else
				{
					chrNameIntVec_2.push_back(tmpChrNameInt);
					startPosVec_2.push_back(tmpStartPos);
					jumpCodeVecVec_2.push_back(tmpJumpCodeVec);
					forOrRcmBoolVec_2.push_back(tmp_forOrRcm_bool);
				}
			}
		}
	}
};

class Xenograft_Read_Sam_Info
{
	bool SE_or_PE_bool;
	Xenograft_Read_Info readInfo;
	Xenograft_Sam_Info samInfo_graft;	
	Xenograft_Sam_Info samInfo_host;
private:
	Xenograft_Read_Sam_Info()
	{

	}

	void initiate_SE(vector<string>& samStrVec_graft, Index_Info* indexInfo_graft,
		vector<string>& samStrVec_host, Index_Info* indexInfo_host)
	{
		SE_or_PE_bool = true;
		int samStrVec_graft_size = samStrVec_graft.size();
		int samStrVec_host_size = samStrVec_host.size();
		if(samStrVec_graft_size + samStrVec_host_size == 0)
		{
			cout << "error ! samStrVec_graft_size + samStrVec_host_size == 0 " << endl;
			exit(1);
		}
		// initaite readInfo
		if(samStrVec_graft_size > 0)
			readInfo.initiate_readInfo_SE(samStrVec_graft[0]);
		else if(samStrVec_host_size > 0)
			readInfo.initiate_readInfo_SE(samStrVec_host[0]);
		else
		{
			cout << "error in initiate_SE of Xenograft_Read_Sam_Info: "
				<< samStrVec_graft_size << "\t" << samStrVec_host_size <<  endl;
			exit(1);
		}
		// initiate samInfo
		// initiate samInfo_graft
		if(samStrVec_graft_size > 0)
			samInfo_graft.initiate_samInfo_SE(samStrVec_graft, indexInfo_graft);
		// initiate samInfo_host
		if(samStrVec_host_size > 0)
			samInfo_host.initiate_samInfo_SE(samStrVec_host, indexInfo_host);
	}

	void initiate_PE(vector<string>& samStrVec_host, Index_Info* indexInfo_host, 
		vector<string>& samStrVec_graft, Index_Info* indexInfo_graft)
	{
		SE_or_PE_bool = false;
		int samStrVec_graft_size = samStrVec_graft.size();
		int samStrVec_host_size = samStrVec_host.size();
		if(samStrVec_graft_size + samStrVec_host_size < 2)
		{
			cout << "error ! samStrVec_graft_size + samStrVec_host_size == " 
				<< samStrVec_graft_size + samStrVec_host_size << endl;
			exit(1);
		}			
		int samStrVec_graft_size_dividedBy2 = samStrVec_graft_size/2;
		int samStrVec_host_size_dividedBy2 = samStrVec_host_size/2;
		if((samStrVec_graft_size_dividedBy2 * 2 != samStrVec_graft_size)
			||(samStrVec_host_size_dividedBy2 * 2 != samStrVec_host_size))
		{
			cout << "(samStrVec_graft_size_dividedBy2 * 2 != samStrVec_graft_size)" << endl;
			cout << "or (samStrVec_host_size_dividedBy2 * 2 != samStrVec_host_size)" << endl;
			cout << "samStrVec_graft_size: " << samStrVec_graft_size << endl;
			cout << "samStrVec_host_size: " << samStrVec_host_size << endl;
			exit(1);
		}
		// initaite readInfo
		if(samStrVec_graft_size >= 2)
			readInfo.initiate_readInfo_PE(samStrVec_graft);
		else if(samStrVec_host_size >= 2)
			readInfo.initiate_readInfo_PE(samStrVec_host);
		else
		{
			cout << "error in initiate_SE of Xenograft_Read_Sam_Info: "
				<< samStrVec_graft_size << "\t" << samStrVec_host_size <<  endl;
			exit(1);
		}

		// initiate samInfo
		// initiate samInfo_graft
		if(samStrVec_graft_size >= 2)
			samInfo_graft.initiate_samInfo_PE(samStrVec_graft, indexInfo_graft);
		// initiate samInfo_host
		if(samStrVec_host_size >= 2)
			samInfo_host.initiate_samInfo_PE(samStrVec_host, indexInfo_host);
	}
};
#endif