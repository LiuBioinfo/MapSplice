#ifndef PESAM_INFO_H
#define PESAM_INFO_H

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

class PeSam_Info
{
private:
	int chrNameInt;

	string readName_upstreamEnd;
	int flag_upstreamEnd;
	int startPos_upstreamEnd;
	vector<Jump_Code> jumpCodeVec_upstreamEnd;
	string mapSeq_upstreamEnd;
	string mapQual_upstreamEnd;
	int NM_upstreamEnd;
	int IH_upstreamEnd;
	int HI_upstreamEnd;
	string samOtherField_upstreamEnd;
	string oriSam_upstreamEnd;

	string readName_downstreamEnd;
	int flag_downstreamEnd;
	int startPos_downstreamEnd;
	vector<Jump_Code> jumpCodeVec_downstreamEnd;
	string mapSeq_downstreamEnd;
	string mapQual_downstreamEnd;
	int NM_downstreamEnd;
	int IH_downstreamEnd;
	int HI_downstreamEnd;
	string samOtherField_downstreamEnd;
	string oriSam_downstreamEnd;

	bool end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool;
public:
	PeSam_Info()
	{}

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

	int returnEndPos_downstreamRead()
	{
		int lastMatchJumpCodeIndex = jumpCodeVec_downstreamEnd.size() - 2;
		if((jumpCodeVec_downstreamEnd[lastMatchJumpCodeIndex].type != "M")
			||(jumpCodeVec_downstreamEnd[lastMatchJumpCodeIndex + 1].type != "S"))
		{
			cout << "error! not M for the lastMatchJumpCodeIndex" << endl;
			cout << "error! not S for the lastJumpCodeIndex" << endl;
			exit(1);
		}
		int lastPos = this->getEndPosOfSpecificJumpCode(startPos_downstreamEnd, 
			jumpCodeVec_downstreamEnd, lastMatchJumpCodeIndex);
		return lastPos;
	}

	int returnStartPos_upstreamRead()
	{
		return startPos_upstreamEnd;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int return_flag_upstreamEnd()
	{
		return flag_upstreamEnd;
	}

	int return_flag_downstreamEnd()
	{
		return flag_downstreamEnd;
	}

	string return_oriSam_upstreamEnd()
	{
		return oriSam_upstreamEnd;
	}

	string return_oriSam_downstreamEnd()
	{
		return oriSam_downstreamEnd;
	}

	void update_mainBody_chrPos_jumpCodeVec_fixFusionAtUpstreamEndHead(int tmpToFixFusionSeq_acceptorAnchorLength,
		int& tmp_updatedMainBodySam_chrPos_upstreamEnd, vector<Jump_Code>& tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd)
	{
		int oriSoftClipHeadLength_upstreamEndSam = jumpCodeVec_upstreamEnd[0].len;
		int oriFirstMatchLength_upstreamEndSam = jumpCodeVec_upstreamEnd[1].len;
		if(tmpToFixFusionSeq_acceptorAnchorLength < 0)
		{
			cout << "tmpToFixFusionSeq_acceptorAnchorLength < 0, error !" << endl;
			exit(1);
		}
		else if(tmpToFixFusionSeq_acceptorAnchorLength > (oriSoftClipHeadLength_upstreamEndSam + oriFirstMatchLength_upstreamEndSam))
		{
			cout << "tmpToFixFusionSeq_acceptorAnchorLength > (oriSoftClipHeadLength_upstreamEndSam + oriFirstMatchLength_upstreamEndSam)" << endl;
			exit(1);
		}
		else
		{
			tmp_updatedMainBodySam_chrPos_upstreamEnd = startPos_upstreamEnd + oriFirstMatchLength_upstreamEndSam - tmpToFixFusionSeq_acceptorAnchorLength;			
			int updatedSoftClippedHeadLength_upstreamEndSam = oriSoftClipHeadLength_upstreamEndSam 
				+ oriFirstMatchLength_upstreamEndSam - tmpToFixFusionSeq_acceptorAnchorLength;
			if(updatedSoftClippedHeadLength_upstreamEndSam != 0)
			{
				Jump_Code updatedSoftClipJumpCode(updatedSoftClippedHeadLength_upstreamEndSam, "S");
				tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd.push_back(updatedSoftClipJumpCode);				
			}
			if(tmpToFixFusionSeq_acceptorAnchorLength != 0)
			{
				Jump_Code acceptorJumpCode(tmpToFixFusionSeq_acceptorAnchorLength, "M");
				tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd.push_back(acceptorJumpCode);
			}
			for(int tmp = 2; tmp < jumpCodeVec_upstreamEnd.size(); tmp++)
				tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd.push_back(jumpCodeVec_upstreamEnd[tmp]);
		}
	}

	void update_mainBody_jumpCodeVec_fixFusionAtDownstreamEndTail(int tmpToFixFusionSeq_donerAnchorLength,
		vector<Jump_Code>& tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd)
	{
		int oriSoftClipTailLength_downstreamEndSam = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].len;
		int oriLastMatchLength_downstreamEndSam = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 2].len;
		if(tmpToFixFusionSeq_donerAnchorLength < 0)
		{
			cout << "tmpToFixFusionSeq_donerAnchorLength < 0, error !" << endl;
			exit(1);
		}
		else if(tmpToFixFusionSeq_donerAnchorLength > (oriSoftClipTailLength_downstreamEndSam + oriLastMatchLength_downstreamEndSam))
		{
			cout << "tmpToFixFusionSeq_donerAnchorLength > (oriSoftClipTailLength_downstreamEndSam + oriLastMatchLength_downstreamEndSam)" << endl;
			exit(1);
		}
		else
		{
			for(int tmp = 0; tmp < jumpCodeVec_downstreamEnd.size() - 2; tmp++)
				tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd.push_back(jumpCodeVec_downstreamEnd[tmp]);
			int updatedSoftClippedTailLength_downstreamEndSam = oriSoftClipTailLength_downstreamEndSam
				+ oriLastMatchLength_downstreamEndSam - tmpToFixFusionSeq_donerAnchorLength;
			if(tmpToFixFusionSeq_donerAnchorLength != 0)
			{
				Jump_Code donerJumpCode(tmpToFixFusionSeq_donerAnchorLength, "M");
				tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd.push_back(donerJumpCode);
			}
			if(updatedSoftClippedTailLength_downstreamEndSam != 0)
			{
				Jump_Code updatedSoftClipJumpCode(updatedSoftClippedTailLength_downstreamEndSam, "S");
				tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd.push_back(updatedSoftClipJumpCode);
			}
		}
	}

	string return_mapSeq_upstreamEnd()
	{
		return mapSeq_upstreamEnd;
	}

	string return_mapSeq_downstreamEnd()
	{
		return mapSeq_downstreamEnd;
	}

	string return_mapQual_upstreamEnd()
	{
		return mapQual_upstreamEnd;
	}

	string return_mapQual_downstreamEnd()
	{
		return mapQual_downstreamEnd;
	}

	int return_toFixFusionSeqLength_atUpstreamEndHead()
	{
		return jumpCodeVec_upstreamEnd[0].len + jumpCodeVec_upstreamEnd[1].len;
	}

	int return_toFixFusionSeqLength_atDownstreamEndTail()
	{
		return jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 2].len
			+ jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].len;
	}

	string returnReadName_upstreamEnd()
	{
		return readName_upstreamEnd;
	}

	string returnReadName_downstreamEnd()
	{
		return readName_downstreamEnd;
	}	

	string returnUnfixedUpstreamReadHead_readSeq()
	{
		string firstJumpCodeType = jumpCodeVec_upstreamEnd[0].type;
		int firstJumpCodeLen = jumpCodeVec_upstreamEnd[0].len;
		return mapSeq_upstreamEnd.substr(0, firstJumpCodeLen);
	}

	string returnUnfixedDownstreamReadTail_readSeq()
	{
		string lastJumpCodeType = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].type;
		int lastJumpCodeLen = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].len;
		return mapSeq_downstreamEnd.substr(mapSeq_downstreamEnd.length() - lastJumpCodeLen);
	}

	string returnUnfixedUpstreamReadHead_qualSeq(bool fasta_or_fastq_bool)
	{
		if(fasta_or_fastq_bool)
			return "*";
		string firstJumpCodeType = jumpCodeVec_upstreamEnd[0].type;
		int firstJumpCodeLen = jumpCodeVec_upstreamEnd[0].len;		
		return mapQual_upstreamEnd.substr(0, firstJumpCodeLen);
	}

	string returnUnfixedDownstreamReadTail_qualSeq(bool fasta_or_fastq_bool)
	{
		if(fasta_or_fastq_bool)
			return "*";
		string lastJumpCodeType = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].type;
		int lastJumpCodeLen = jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size() - 1].len;
		return mapQual_downstreamEnd.substr(mapQual_downstreamEnd.length() - lastJumpCodeLen);
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
				tmpJumpCodeLength = atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	void pushJumpCodeVec2target_upstreamEnd(vector<Jump_Code>& targetJumpCodeVec)
	{
		for(int tmp = 0; tmp < jumpCodeVec_upstreamEnd.size(); tmp++)
			targetJumpCodeVec.push_back(jumpCodeVec_upstreamEnd[tmp]);
	}

	void pushJumpCodeVec2target_downstreamEnd(vector<Jump_Code>& targetJumpCodeVec)
	{
		for(int tmp = 0; tmp < jumpCodeVec_downstreamEnd.size(); tmp++)
			targetJumpCodeVec.push_back(jumpCodeVec_downstreamEnd[tmp]);		
	}

	int return_startPos_upstreamEnd()
	{
		return startPos_upstreamEnd;
	}

	int return_startPos_downstreamEnd()
	{
		return startPos_downstreamEnd;
	}

	bool return_end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool()
	{
		return end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool;
	}

	bool unique_or_not_bool()
	{
		if((HI_upstreamEnd == 1)&&(HI_downstreamEnd == 1))
			return true;
		else
			return false;
	}

    string returnRawReadName_1()
    {
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return readName_upstreamEnd;
    	else
    		return readName_downstreamEnd;
    }

    string returnRawReadName_2()
    {
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return readName_downstreamEnd;
    	else
    		return readName_upstreamEnd;    	
    }

    string returnRawReadSeq_1() 
	{
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return mapSeq_upstreamEnd;
    	else
    		return covertStringToReverseComplement(mapSeq_downstreamEnd);
	}

	string returnRawReadSeq_2()
	{
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return covertStringToReverseComplement(mapSeq_downstreamEnd);
    	else
    		return mapSeq_upstreamEnd;    				
	}

	string returnRawQualSeq_1()
	{
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return mapQual_upstreamEnd;
    	else
    		return convertQualityScoreString2Reverse(mapQual_downstreamEnd);
	}

	string returnRawQualSeq_2()
	{
    	if(end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool)
    		return mapQual_downstreamEnd;
    	else
    		return convertQualityScoreString2Reverse(mapQual_upstreamEnd);
	}

	bool end1forMapEnd2RevMap_or_end2forMapEnd1RevMap(int tmpFlag_upstreamEnd, int tmpFlag_downstreamEnd)
	{
		bool upstreamEnd_end1_or_end2_bool = (tmpFlag_upstreamEnd & 0x40);
		bool downstreamEnd_end1_or_end2_bool = (tmpFlag_downstreamEnd & 0x40);
		if((upstreamEnd_end1_or_end2_bool)&&(!downstreamEnd_end1_or_end2_bool))
			return true;
		else if((!upstreamEnd_end1_or_end2_bool)&&(downstreamEnd_end1_or_end2_bool))
			return false;
		else
		{
			cout << "error in end1forMapEnd2RevMap_or_end2forMapEnd1RevMap" << endl;
			cout << "tmpFlag_upstreamEnd: " << tmpFlag_upstreamEnd << endl;
			cout << "tmpFlag_downstreamEnd: " << tmpFlag_downstreamEnd << endl;
			exit(1);
		}
	}

	bool initiateWith2samStr(string& sam_upstreamEnd, string& sam_downstreamEnd, Index_Info* indexInfo)
	{
		oriSam_upstreamEnd = sam_upstreamEnd;
		oriSam_downstreamEnd = sam_downstreamEnd;

		vector<string> samFieldVec_upstreamEnd;
		vector<string> samFieldVec_downstreamEnd;
		int startLoc = 0;
		for(int tmp = 0; tmp < 14; tmp++)
		{
			int tabLoc = sam_upstreamEnd.find("\t", startLoc);
			if(tabLoc == string::npos)
				return false;
			string tmpSamField = sam_upstreamEnd.substr(startLoc, tabLoc-startLoc);
			samFieldVec_upstreamEnd.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samOtherField_upstreamEnd = sam_upstreamEnd.substr(startLoc);

		startLoc = 0;
		for(int tmp = 0; tmp < 14; tmp++)
		{
			int tabLoc = sam_downstreamEnd.find("\t", startLoc);
			if(tabLoc == string::npos)
				return false;
			string tmpSamField = sam_downstreamEnd.substr(startLoc, tabLoc-startLoc);
			samFieldVec_downstreamEnd.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samOtherField_downstreamEnd = sam_downstreamEnd.substr(startLoc);

		string chrNameStr_upstreamEnd = samFieldVec_upstreamEnd[2];
		string chrNameStr_downstreamEnd = samFieldVec_downstreamEnd[2];
		if(chrNameStr_upstreamEnd != chrNameStr_downstreamEnd)
			return false;
		chrNameInt = indexInfo->convertStringToInt(chrNameStr_upstreamEnd);
		if(chrNameInt < 0)
			return false;

		readName_upstreamEnd = samFieldVec_upstreamEnd[0];
		readName_downstreamEnd = samFieldVec_downstreamEnd[0];

		string flagStr_upstreamEnd = samFieldVec_upstreamEnd[1];
		string flagStr_downstreamEnd = samFieldVec_downstreamEnd[1];
		flag_upstreamEnd = atoi(flagStr_upstreamEnd.c_str());
		flag_downstreamEnd = atoi(flagStr_downstreamEnd.c_str());		

		string startPosStr_upstreamEnd = samFieldVec_upstreamEnd[3];
		string startPosStr_downstreamEnd = samFieldVec_downstreamEnd[3];
		startPos_upstreamEnd = atoi(startPosStr_upstreamEnd.c_str());
		startPos_downstreamEnd = atoi(startPosStr_downstreamEnd.c_str());

		string cigarString_upstreamEnd = samFieldVec_upstreamEnd[5];
		string cigarString_downstreamEnd = samFieldVec_downstreamEnd[5];
		this->cigarString2jumpCodeVec(cigarString_upstreamEnd, jumpCodeVec_upstreamEnd);
		this->cigarString2jumpCodeVec(cigarString_downstreamEnd, jumpCodeVec_downstreamEnd);

		mapSeq_upstreamEnd = samFieldVec_upstreamEnd[9];
		mapSeq_downstreamEnd = samFieldVec_downstreamEnd[9];

		mapQual_upstreamEnd = samFieldVec_upstreamEnd[10];
		mapQual_downstreamEnd = samFieldVec_downstreamEnd[10];

		string NMstr_upstreamEnd = (samFieldVec_upstreamEnd[11]).substr(5);
		string NMstr_downstreamEnd = (samFieldVec_downstreamEnd[11]).substr(5);
		NM_upstreamEnd = atoi(NMstr_upstreamEnd.c_str());
		NM_downstreamEnd = atoi(NMstr_downstreamEnd.c_str());

		string IHstr_upstreamEnd = (samFieldVec_upstreamEnd[12]).substr(5);
		string IHstr_downstreamEnd = (samFieldVec_downstreamEnd[12]).substr(5);
		IH_upstreamEnd = atoi(IHstr_upstreamEnd.c_str());
		IH_downstreamEnd = atoi(IHstr_downstreamEnd.c_str());

		string HIstr_upstreamEnd = (samFieldVec_upstreamEnd[13]).substr(5);
		string HIstr_downstreamEnd = (samFieldVec_downstreamEnd[13]).substr(5);
		HI_upstreamEnd = atoi(HIstr_upstreamEnd.c_str());
		HI_downstreamEnd = atoi(HIstr_downstreamEnd.c_str());

		end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool = end1forMapEnd2RevMap_or_end2forMapEnd1RevMap(
			flag_upstreamEnd, flag_downstreamEnd);
		return true;
	}

	bool return_unfixedHeadAtUpstreamEnd_exists_bool()
	{
		return (jumpCodeVec_upstreamEnd[0].type == "S");
	}

	bool return_unfixedTailAtDownstreamEnd_exists_bool()
	{
		return (jumpCodeVec_downstreamEnd[jumpCodeVec_downstreamEnd.size()-1].type == "S");
	}
};
#endif