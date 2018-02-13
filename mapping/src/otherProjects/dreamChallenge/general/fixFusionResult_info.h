#ifndef FIXFUSIONRESULT_INFO_H
#define FIXFUSIONRESULT_INFO_H

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

class FixFusionResult_Info
{
private:
	bool SE_or_PE_bool;

	// PE
	bool fusionFixed_unfixedHeadAtUpstreamEnd_success_bool;
	int fusionFixed_unfixedLeftReadHead_donerAnchorLength;
	int fusionFixed_unfixedLeftReadHead_foreignInsertionLength;
	int fusionFixed_unfixedLeftReadHead_acceptorAnchorLength;
	int fusionFixed_unfixedLeftReadHead_mismatchNum;
	int fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt;
	int fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos;
	bool fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool;
	vector<Jump_Code> fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec;
	string fusionFixed_unfixedLeftReadHead_clippedSeg_readSeq;
	string fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq;

	bool fusionFixed_unfixedTailAtDownstreamEnd_success_bool;
	int fusionFixed_unfixedRightReadTail_donerAnchorLength;
	int fusionFixed_unfixedRightReadTail_foreignInsertionLength;
	int fusionFixed_unfixedRightReadTail_acceptorAnchorLength;
	int fusionFixed_unfixedRightReadTail_mismatchNum;
	int fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt;
	int fusionFixed_unfixedRightReadTail_clippedSeg_chrPos;
	bool fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool;
	vector<Jump_Code> fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec;
	string fusionFixed_unfixedRightReadTail_clippedSeg_readSeq;
	string fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq;

	// SE
	bool fusionFixed_unfixedHeadAtSEread_success_bool;
	bool fusionFixed_unfixedTailAtSEread_success_bool;
public:
	FixFusionResult_Info()
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

	int getEndPosOfJumpCodeVec(int tmpStartPos, vector<Jump_Code>& tmpJumpCodeVec)
	{
		int lastJumpCodeIndex = tmpJumpCodeVec.size() - 1;
		return getEndPosOfSpecificJumpCode(tmpStartPos, tmpJumpCodeVec, lastJumpCodeIndex);
	}

	void update_unfixedEndSam_mainBodySam_breakPointPairStr(PeSam_Info& tmpPeSamInfo, PE_Read_Info& tmpPeReadInfo,
		string& outputClippedSam_unfixedHeadAtUpstreamEnd, string& outputClippedSam_unfixedTailAtDownstreamEnd,
		string& outputMainBodySam_upstreamEnd, string& outputMainBodySam_downstreamEnd,
		string& outputBreakPointPair_upstreamEnd, string& outputBreakPointPair_downstreamEnd, Index_Info* indexInfo, 
		GeneAnnEntry_Hash_Info& tmpGeneAnnEntryHashInfo)
	{
		if(fusionFixed_unfixedHeadAtUpstreamEnd_success_bool)
		{
			int tmpToFixFusionSeq_foreighInsertionLength = fusionFixed_unfixedLeftReadHead_foreignInsertionLength;
			int tmpToFixFusionSeq_donerAnchorLength = fusionFixed_unfixedLeftReadHead_donerAnchorLength;
			int tmpToFixFusionSeq_acceptorAnchorLength = fusionFixed_unfixedLeftReadHead_acceptorAnchorLength;
			int tmpToFixFusionSeq_mismatchNum = fusionFixed_unfixedLeftReadHead_mismatchNum;
			//cout << "fusionFixed_unfixedHeadAtUpstreamEnd_success_bool: " << fusionFixed_unfixedHeadAtUpstreamEnd_success_bool << endl;
			//cout << "tmpToFixFusionSeq_donerAnchorLength: " << tmpToFixFusionSeq_donerAnchorLength << endl;
			//cout << "tmpToFixFusionSeq_acceptorAnchorLength: " << tmpToFixFusionSeq_acceptorAnchorLength << endl;
			//cout << "tmpToFixFusionSeq_foreighInsertionLength: " << tmpToFixFusionSeq_foreighInsertionLength << endl;
			//cout << "tmpToFixFusionSeq_mismatchNum: " << tmpToFixFusionSeq_mismatchNum << endl;
			string tmpReadName = tmpPeSamInfo.returnReadName_upstreamEnd();
			int tmpFlag_unfixedLeftReadHead_clippedSeg;
			if(fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool)
				tmpFlag_unfixedLeftReadHead_clippedSeg = 0;
			else
				tmpFlag_unfixedLeftReadHead_clippedSeg = 16;
			
			outputClippedSam_unfixedHeadAtUpstreamEnd = tmpReadName
				+ "\t" + int_to_str(tmpFlag_unfixedLeftReadHead_clippedSeg)
				+ "\t" + indexInfo->returnChrNameStr(fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt) 
				+ "\t" + int_to_str(fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos) 
				+ "\t255\t" + this->jumpCodeVec2Str(fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec) 
				+ "\t*\t0\t0\t" + fusionFixed_unfixedLeftReadHead_clippedSeg_readSeq 
				+ "\t" + fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq + "\tIH:i:1\tHI:i:1";
			
			int tmp_updatedMainBodySam_chrPos_upstreamEnd;
			vector<Jump_Code> tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd;
			tmpPeSamInfo.update_mainBody_chrPos_jumpCodeVec_fixFusionAtUpstreamEndHead(tmpToFixFusionSeq_acceptorAnchorLength,
				tmp_updatedMainBodySam_chrPos_upstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd);
			
			outputMainBodySam_upstreamEnd = tmpReadName
				+ "\t" + int_to_str(tmpPeSamInfo.return_flag_upstreamEnd())
				+ "\t" + indexInfo->returnChrNameStr(tmpPeSamInfo.returnChrNameInt()) 
				+ "\t" + int_to_str(tmp_updatedMainBodySam_chrPos_upstreamEnd) 
				+ "\t255\t" + this->jumpCodeVec2Str(tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd) 
				+ "\t*\t0\t0\t" + tmpPeSamInfo.return_mapSeq_upstreamEnd() 
				+ "\t" + tmpPeSamInfo.return_mapQual_upstreamEnd() + "\tIH:i:1\tHI:i:1";
			
			outputMainBodySam_downstreamEnd = tmpPeSamInfo.return_oriSam_downstreamEnd();
			
			int tmpFusionBreakPointAtClippedSeg, tmpFusionMaxRangeOtherEndPos_mainBody, tmpFusionMaxRangeOtherEndPos_clippedSeg;
			string tmpFusionFlankString, tmpFusionStrand_left, tmpFusionStrand_right, tmpFusionCaseStr;
			if(fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool) // case 2,5
			{	
				tmpFusionBreakPointAtClippedSeg = this->getEndPosOfJumpCodeVec(fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos, fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec);
				tmpFusionMaxRangeOtherEndPos_mainBody = this->getEndPosOfJumpCodeVec(tmp_updatedMainBodySam_chrPos_upstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd);
				tmpFusionMaxRangeOtherEndPos_clippedSeg = fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos;
				tmpFusionFlankString = indexInfo->returnFormattedFusionJuncFlankString(fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt,
					tmpPeSamInfo.returnChrNameInt(), tmpFusionBreakPointAtClippedSeg, tmp_updatedMainBodySam_chrPos_upstreamEnd, true, true);

				string tmpFusion_leftGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt, tmpFusionBreakPointAtClippedSeg, indexInfo);
				string tmpFusion_rightGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					tmpPeSamInfo.returnChrNameInt(), tmp_updatedMainBodySam_chrPos_upstreamEnd, indexInfo);
				string tmpFusion_strand_inGeneAnnEntryHash = tmpFusion_leftGene_strand_inGeneAnnEntryHash + tmpFusion_rightGene_strand_inGeneAnnEntryHash;
				
				//breakPointDetection_valid_orNot_bool = true;
				if((tmpFusion_leftGene_strand_inGeneAnnEntryHash == "O")||(tmpFusion_rightGene_strand_inGeneAnnEntryHash == "O"))
				{
					tmpFusionStrand_left = "O";
					tmpFusionStrand_right = "O";
					tmpFusionCaseStr = "2,5,";					
				}
				// else if((tmpFusion_strand_inGeneAnnEntryHash == "+-")||(tmpFusion_strand_inGeneAnnEntryHash == "-+"))
				// {
				// 	tmpFusionStrand_left = "X";
				// 	tmpFusionStrand_right = "X";
				// 	tmpFusionCaseStr = "2,5,";					
				// }
				else if((tmpFusion_strand_inGeneAnnEntryHash == "++")||(tmpFusion_strand_inGeneAnnEntryHash == "X+")||(tmpFusion_strand_inGeneAnnEntryHash == "+X"))
				{
					tmpFusionStrand_left = "+";
					tmpFusionStrand_right = "+";
					tmpFusionCaseStr = "2,";					
				}
				else if((tmpFusion_strand_inGeneAnnEntryHash == "--")||(tmpFusion_strand_inGeneAnnEntryHash == "--")||(tmpFusion_strand_inGeneAnnEntryHash == "--"))
				{
					tmpFusionStrand_left = "-";
					tmpFusionStrand_right = "-";
					tmpFusionCaseStr = "5,";					
				}
				else // tmpFusion_strand_inGeneAnnEntryHash == XX
				{
					if((tmpFusionFlankString == "GTAG")||(tmpFusionFlankString == "GCAG")||(tmpFusionFlankString == "ATAC")) // case 2
					{
						tmpFusionStrand_left = "+";
						tmpFusionStrand_right = "+";
						tmpFusionCaseStr = "2,";
					}
					else if((tmpFusionFlankString == "CTAC")||(tmpFusionFlankString == "CTGC")||(tmpFusionFlankString == "GTAT")) // case 5
					{
						tmpFusionStrand_left = "-";
						tmpFusionStrand_right = "-";
						tmpFusionCaseStr = "5,";
					}
					else // case 2,5
					{
						tmpFusionStrand_left = "N";
						tmpFusionStrand_right = "N";
						tmpFusionCaseStr = "2,5,";
					}
				}
			}
			else // case 10,11
			{
				tmpFusionBreakPointAtClippedSeg = fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos;
				tmpFusionMaxRangeOtherEndPos_mainBody = this->getEndPosOfJumpCodeVec(tmp_updatedMainBodySam_chrPos_upstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_upstreamEnd);
				tmpFusionMaxRangeOtherEndPos_clippedSeg = this->getEndPosOfJumpCodeVec(fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos, fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec);
				tmpFusionFlankString = indexInfo->returnFormattedFusionJuncFlankString(fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt,
					tmpPeSamInfo.returnChrNameInt(), tmpFusionBreakPointAtClippedSeg, tmp_updatedMainBodySam_chrPos_upstreamEnd, true, false);
				
				string tmpFusion_leftGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt, tmpFusionBreakPointAtClippedSeg, indexInfo);
				string tmpFusion_rightGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					tmpPeSamInfo.returnChrNameInt(), tmp_updatedMainBodySam_chrPos_upstreamEnd, indexInfo);
				string tmpFusion_strand_inGeneAnnEntryHash = tmpFusion_leftGene_strand_inGeneAnnEntryHash + tmpFusion_rightGene_strand_inGeneAnnEntryHash;
				
				if((tmpFusion_leftGene_strand_inGeneAnnEntryHash == "O")||(tmpFusion_rightGene_strand_inGeneAnnEntryHash == "O")) // not in any gene
				{
					tmpFusionStrand_left = "O";
					tmpFusionStrand_right = "O";
					tmpFusionCaseStr = "10,11,";					
				}			
				// else if((tmpFusion_strand_inGeneAnnEntryHash == "++")||(tmpFusion_strand_inGeneAnnEntryHash == "--"))// conflict
				// {
				// 	tmpFusionStrand_left = "X";
				// 	tmpFusionStrand_right = "X";
				// 	tmpFusionCaseStr = "10,11,";					
				// }
				else if((tmpFusion_strand_inGeneAnnEntryHash == "+-")||(tmpFusion_strand_inGeneAnnEntryHash == "+X")||(tmpFusion_strand_inGeneAnnEntryHash == "X-")) // case 11
				{
					tmpFusionStrand_left = "+";
					tmpFusionStrand_right = "-";
					tmpFusionCaseStr = "11,";					
				}
				else if((tmpFusion_strand_inGeneAnnEntryHash == "-+")||(tmpFusion_strand_inGeneAnnEntryHash == "-X")||(tmpFusion_strand_inGeneAnnEntryHash == "X+"))
				{
					tmpFusionStrand_left = "-";
					tmpFusionStrand_right = "+";
					tmpFusionCaseStr = "10,";
				}
				else // tmpFusion_strand_inGeneAnnEntryHash == XX
				{	
					if((tmpFusionFlankString == "GTAG")||(tmpFusionFlankString == "GCAG")||(tmpFusionFlankString == "ATAC")) // case 10
					{
						tmpFusionStrand_left = "-";
						tmpFusionStrand_right = "+";
						tmpFusionCaseStr = "10,";
					}
					else if((tmpFusionFlankString == "CTAC")||(tmpFusionFlankString == "CTGC")||(tmpFusionFlankString == "GTAT")) // case 11
					{
						tmpFusionStrand_left = "+";
						tmpFusionStrand_right = "-";
						tmpFusionCaseStr = "11,";
					}
					else // case 10,11
					{
						tmpFusionStrand_left = "N";
						tmpFusionStrand_right = "N";
						tmpFusionCaseStr = "10,11,";
					}
				}
			}
			
			outputBreakPointPair_upstreamEnd = tmpReadName
				+ "\t" + indexInfo->returnChrNameStr(fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt)
				+ "\t" + indexInfo->returnChrNameStr(tmpPeSamInfo.returnChrNameInt()) 
				+ "\t" + int_to_str(tmpFusionBreakPointAtClippedSeg)
				+ "\t" + int_to_str(tmp_updatedMainBodySam_chrPos_upstreamEnd)
				+ "\t" + tmpFusionStrand_left + "\t" + tmpFusionStrand_right + "\t" + tmpFusionFlankString 
				+ "\t" + int_to_str(tmpToFixFusionSeq_donerAnchorLength) + "\t" + int_to_str(tmpToFixFusionSeq_acceptorAnchorLength)
				+ "\t" + int_to_str(tmpToFixFusionSeq_foreighInsertionLength) 
				+ "\t" + tmpFusionCaseStr + "\t" + int_to_str(tmpFusionMaxRangeOtherEndPos_clippedSeg) + "\t" + int_to_str(tmpFusionMaxRangeOtherEndPos_mainBody)
				+ "\t" + int_to_str(tmpToFixFusionSeq_mismatchNum);
		}

		if(fusionFixed_unfixedTailAtDownstreamEnd_success_bool)
		{
			int tmpToFixFusionSeq_foreighInsertionLength = fusionFixed_unfixedRightReadTail_foreignInsertionLength;
			int tmpToFixFusionSeq_donerAnchorLength = fusionFixed_unfixedRightReadTail_donerAnchorLength;
			int tmpToFixFusionSeq_acceptorAnchorLength = fusionFixed_unfixedRightReadTail_acceptorAnchorLength;			
			int tmpToFixFusionSeq_mismatchNum = fusionFixed_unfixedRightReadTail_mismatchNum;
			//cout << "fusionFixed_unfixedTailAtDownstreamEnd_success_bool: " << fusionFixed_unfixedTailAtDownstreamEnd_success_bool << endl;
			//cout << "tmpToFixFusionSeq_donerAnchorLength: " << tmpToFixFusionSeq_donerAnchorLength << endl;
			//cout << "tmpToFixFusionSeq_acceptorAnchorLength: " << tmpToFixFusionSeq_acceptorAnchorLength << endl;
			//cout << "tmpToFixFusionSeq_foreighInsertionLength: " << tmpToFixFusionSeq_foreighInsertionLength << endl;
			//cout << "tmpToFixFusionSeq_mismatchNum: " << tmpToFixFusionSeq_mismatchNum << endl;
			string tmpReadName = tmpPeSamInfo.returnReadName_downstreamEnd();
			int tmpFlag_unfixedRightReadTail_clippedSeg;
			if(fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool)
				tmpFlag_unfixedRightReadTail_clippedSeg = 0;
			else
				tmpFlag_unfixedRightReadTail_clippedSeg = 16;
			
			outputClippedSam_unfixedTailAtDownstreamEnd = tmpReadName
				+ "\t" + int_to_str(tmpFlag_unfixedRightReadTail_clippedSeg)
				+ "\t" + indexInfo->returnChrNameStr(fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt)
				+ "\t" + int_to_str(fusionFixed_unfixedRightReadTail_clippedSeg_chrPos)
				+ "\t255\t" + this->jumpCodeVec2Str(fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec)
				+ "\t*\t0\t0\t" + fusionFixed_unfixedRightReadTail_clippedSeg_readSeq
				+ "\t" + fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq + "\tIH:i:1\tHI:i:1";
			
			int tmp_chrPos_downstreamEnd = tmpPeSamInfo.return_startPos_downstreamEnd();
			vector<Jump_Code> tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd;
			tmpPeSamInfo.update_mainBody_jumpCodeVec_fixFusionAtDownstreamEndTail(tmpToFixFusionSeq_donerAnchorLength,
				tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd);
			
			outputMainBodySam_upstreamEnd = tmpPeSamInfo.return_oriSam_upstreamEnd();
			
			outputMainBodySam_downstreamEnd = tmpReadName
				+ "\t" + int_to_str(tmpPeSamInfo.return_flag_downstreamEnd())
				+ "\t" + indexInfo->returnChrNameStr(tmpPeSamInfo.returnChrNameInt())
				+ "\t" + int_to_str(tmp_chrPos_downstreamEnd)
				+ "\t255\t" + this->jumpCodeVec2Str(tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd)
				+ "\t*\t0\t0\t" + tmpPeSamInfo.return_mapSeq_downstreamEnd() 
				+ "\t" + tmpPeSamInfo.return_mapQual_downstreamEnd() + "\tIH:i:1\tHI:i:1";

			int tmpFusionBreakPointAtClippedSeg, tmpFusionMaxRangeOtherEndPos_mainBody, tmpFusionMaxRangeOtherEndPos_clippedSeg;
			string tmpFusionFlankString, tmpFusionStrand_left, tmpFusionStrand_right, tmpFusionCaseStr;
			if(fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool) // case 1,4
			{
				tmpFusionBreakPointAtClippedSeg = fusionFixed_unfixedRightReadTail_clippedSeg_chrPos;
				tmpFusionMaxRangeOtherEndPos_mainBody = tmp_chrPos_downstreamEnd;
				tmpFusionMaxRangeOtherEndPos_clippedSeg = this->getEndPosOfJumpCodeVec(fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec);
				tmpFusionFlankString = indexInfo->returnFormattedFusionJuncFlankString(
					tmpPeSamInfo.returnChrNameInt(), fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt, 
					this->getEndPosOfJumpCodeVec(tmp_chrPos_downstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd),
					fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, false, true);

				string tmpFusion_leftGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					tmpPeSamInfo.returnChrNameInt(), this->getEndPosOfJumpCodeVec(tmp_chrPos_downstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd), indexInfo);
				string tmpFusion_rightGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt, fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, indexInfo);
				string tmpFusion_strand_inGeneAnnEntryHash = tmpFusion_leftGene_strand_inGeneAnnEntryHash + tmpFusion_rightGene_strand_inGeneAnnEntryHash;

				if((tmpFusion_leftGene_strand_inGeneAnnEntryHash == "O")||(tmpFusion_rightGene_strand_inGeneAnnEntryHash == "O"))
				{
					tmpFusionStrand_left = "O";
					tmpFusionStrand_right = "O";
					tmpFusionCaseStr = "1,4,";			
				}
				// else if((tmpFusion_strand_inGeneAnnEntryHash == "+-")||(tmpFusion_strand_inGeneAnnEntryHash == "-+"))
				// {
				// 	tmpFusionStrand_left = "X";
				// 	tmpFusionStrand_right = "X";
				// 	tmpFusionCaseStr = "1,4,";							
				// }
				else if((tmpFusion_strand_inGeneAnnEntryHash == "++")||(tmpFusion_strand_inGeneAnnEntryHash == "+X")||(tmpFusion_strand_inGeneAnnEntryHash == "X+"))
				{
					tmpFusionStrand_left = "+";
					tmpFusionStrand_right = "+";
					tmpFusionCaseStr = "1,";
				}
				else if((tmpFusion_strand_inGeneAnnEntryHash == "--")||(tmpFusion_strand_inGeneAnnEntryHash == "-X")||(tmpFusion_strand_inGeneAnnEntryHash == "X-"))
				{
					tmpFusionStrand_left = "-";
					tmpFusionStrand_right = "-";
					tmpFusionCaseStr = "4,";					
				}
				else // tmpFusion_strand_inGeneAnnEntryHash == XX
				{	
					if((tmpFusionFlankString == "GTAG")||(tmpFusionFlankString == "GCAG")||(tmpFusionFlankString == "ATAC")) // case 1
					{
						tmpFusionStrand_left = "+";
						tmpFusionStrand_right = "+";
						tmpFusionCaseStr = "1,";
					}
					else if((tmpFusionFlankString == "CTAC")||(tmpFusionFlankString == "CTGC")||(tmpFusionFlankString == "GTAT")) // case 4
					{
						tmpFusionStrand_left = "-";
						tmpFusionStrand_right = "-";
						tmpFusionCaseStr = "4,";
					}
					else // case 1,4
					{
						tmpFusionStrand_left = "N";
						tmpFusionStrand_right = "N";
						tmpFusionCaseStr = "1,4,";
					}
				}
			}
			else // case 7,8
			{
				tmpFusionBreakPointAtClippedSeg = this->getEndPosOfJumpCodeVec(fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec);
				tmpFusionMaxRangeOtherEndPos_mainBody = tmp_chrPos_downstreamEnd;
				tmpFusionMaxRangeOtherEndPos_clippedSeg = fusionFixed_unfixedRightReadTail_clippedSeg_chrPos;
				tmpFusionFlankString = indexInfo->returnFormattedFusionJuncFlankString(
					tmpPeSamInfo.returnChrNameInt(), fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt, 
					this->getEndPosOfJumpCodeVec(tmp_chrPos_downstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd),
					fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, false, false);

				string tmpFusion_leftGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					tmpPeSamInfo.returnChrNameInt(), this->getEndPosOfJumpCodeVec(tmp_chrPos_downstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd), indexInfo);
				string tmpFusion_rightGene_strand_inGeneAnnEntryHash = tmpGeneAnnEntryHashInfo.searchAndReturnGeneAnnEntryStrand_within(
					fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt, fusionFixed_unfixedRightReadTail_clippedSeg_chrPos, indexInfo);
				string tmpFusion_strand_inGeneAnnEntryHash = tmpFusion_leftGene_strand_inGeneAnnEntryHash + tmpFusion_rightGene_strand_inGeneAnnEntryHash;

				if((tmpFusion_leftGene_strand_inGeneAnnEntryHash == "O")||(tmpFusion_rightGene_strand_inGeneAnnEntryHash == "O"))
				{
					tmpFusionStrand_left = "O";
					tmpFusionStrand_right = "O";
					tmpFusionCaseStr = "7,8,";	
				}				
				// else if((tmpFusion_strand_inGeneAnnEntryHash == "++")||(tmpFusion_strand_inGeneAnnEntryHash == "--"))
				// {
				// 	tmpFusionStrand_left = "X";
				// 	tmpFusionStrand_right = "X";
				// 	tmpFusionCaseStr = "7,8,";			
				// }
				else if((tmpFusion_strand_inGeneAnnEntryHash == "+-")||(tmpFusion_strand_inGeneAnnEntryHash == "+X")||(tmpFusion_strand_inGeneAnnEntryHash == "X-"))
				{
					tmpFusionStrand_left = "+";
					tmpFusionStrand_right = "-";
					tmpFusionCaseStr = "8,";					
				}
				else if((tmpFusion_strand_inGeneAnnEntryHash == "-+")||(tmpFusion_strand_inGeneAnnEntryHash == "-X")||(tmpFusion_strand_inGeneAnnEntryHash == "X+"))
				{
					tmpFusionStrand_left = "-";
					tmpFusionStrand_right = "+";
					tmpFusionCaseStr = "7,";					
				}
				else // tmpFusion_strand_inGeneAnnEntryHash == XX
				{	
					if((tmpFusionFlankString == "GTAG")||(tmpFusionFlankString == "GCAG")||(tmpFusionFlankString == "ATAC")) // case 8
					{
						tmpFusionStrand_left = "+";
						tmpFusionStrand_right = "-";
						tmpFusionCaseStr = "8,";
					}
					else if((tmpFusionFlankString == "CTAC")||(tmpFusionFlankString == "CTGC")||(tmpFusionFlankString == "GTAT")) // case 7
					{
						tmpFusionStrand_left = "-";
						tmpFusionStrand_right = "+";
						tmpFusionCaseStr = "7,";
					}
					else // case 1,4
					{
						tmpFusionStrand_left = "N";
						tmpFusionStrand_right = "N";
						tmpFusionCaseStr = "7,8,";
					}
				}
			}
			
			outputBreakPointPair_downstreamEnd = tmpReadName
				+ "\t" + indexInfo->returnChrNameStr(tmpPeSamInfo.returnChrNameInt())
				+ "\t" + indexInfo->returnChrNameStr(fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt)
				+ "\t" + int_to_str(this->getEndPosOfJumpCodeVec(tmp_chrPos_downstreamEnd, tmp_updatedMainBodySam_jumpCodeVec_downstreamEnd))				
				+ "\t" + int_to_str(fusionFixed_unfixedRightReadTail_clippedSeg_chrPos)
				+ "\t" + tmpFusionStrand_left + "\t" + tmpFusionStrand_right + "\t" + tmpFusionFlankString 
				+ "\t" + int_to_str(tmpToFixFusionSeq_donerAnchorLength) + "\t" + int_to_str(tmpToFixFusionSeq_acceptorAnchorLength)
				+ "\t" + int_to_str(tmpToFixFusionSeq_foreighInsertionLength) 
				+ "\t" + tmpFusionCaseStr + "\t" + int_to_str(tmpFusionMaxRangeOtherEndPos_mainBody) + "\t" + int_to_str(tmpFusionMaxRangeOtherEndPos_clippedSeg)
				+ "\t" + int_to_str(tmpToFixFusionSeq_mismatchNum);
			//cout << "outputClippedSam_unfixedTailAtDownstreamEnd: " << outputClippedSam_unfixedTailAtDownstreamEnd << endl;
			//cout << "outputMainBodySam_downstreamEnd: " << outputMainBodySam_downstreamEnd << endl;
			//cout << "outputBreakPointPair_downstreamEnd: " << outputBreakPointPair_downstreamEnd << endl;
		}
	}

	string jumpCodeVec2Str(vector<Jump_Code>& tmpJumpCodeVec)
	{
		string str;
		for(int tmp = 0; tmp < tmpJumpCodeVec.size(); tmp++)
		{	
			str += tmpJumpCodeVec[tmp].toString();
		}
		return str;
	}

	bool return_fusionFixed_unfixedHeadAtUpstreamEnd_success_bool()
	{
		return fusionFixed_unfixedHeadAtUpstreamEnd_success_bool;
	}

	bool return_fusionFixed_unfixedTailAtDownstreamEnd_success_bool()
	{
		return fusionFixed_unfixedTailAtDownstreamEnd_success_bool;
	}

	void initiate_SE()
	{
		fusionFixed_unfixedHeadAtSEread_success_bool = false;
		fusionFixed_unfixedTailAtSEread_success_bool = false;		
	}

	void initiate_PE()
	{
		fusionFixed_unfixedHeadAtUpstreamEnd_success_bool = false;
		fusionFixed_unfixedTailAtDownstreamEnd_success_bool = false;		
	}

	void updateFixedFusionResult_unfixedHeadAtUpstreamEnd(PE_Read_Info& tmpPeReadInfo_unfixedLeftReadHead, 
		PE_Read_Alignment_Info& tmpPeAlignInfo_unfixedLeftReadHead,
		FixFusion_PeRead& tmpFixFusionPeReadInfo_unfixedLeftReadHead, bool fasta_or_fastq_bool, Index_Info* indexInfo)
	{
		//cout << "start to do updateFixedFusionResult_unfixedHeadAtUpstreamEnd( " << endl;
		fusionFixed_unfixedHeadAtUpstreamEnd_success_bool = true;
		fusionFixed_unfixedLeftReadHead_donerAnchorLength = tmpFixFusionPeReadInfo_unfixedLeftReadHead.return_fusionFixed_donerAnchorLength();
		fusionFixed_unfixedLeftReadHead_foreignInsertionLength = tmpFixFusionPeReadInfo_unfixedLeftReadHead.return_fusionFixed_foreignInsertionLength();
		fusionFixed_unfixedLeftReadHead_acceptorAnchorLength = tmpFixFusionPeReadInfo_unfixedLeftReadHead.return_fusionFixed_acceptorAnchorLength();
		fusionFixed_unfixedLeftReadHead_mismatchNum = tmpFixFusionPeReadInfo_unfixedLeftReadHead.return_fusionFixed_mismatchNum();
		fusionFixed_unfixedLeftReadHead_clippedSeg_chrNameInt = tmpPeAlignInfo_unfixedLeftReadHead.returnChrNameInt_SE_uniqAlign(indexInfo);
		fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool = tmpPeAlignInfo_unfixedLeftReadHead.returnForMapOrRcmMap_SE_uniqAlign();
		// cout << "fusionFixed_unfixedLeftReadHead_donerAnchorLength: " << fusionFixed_unfixedLeftReadHead_donerAnchorLength << endl;
		// cout << "fusionFixed_unfixedLeftReadHead_foreignInsertionLength: " << fusionFixed_unfixedLeftReadHead_foreignInsertionLength << endl;
		// cout << "fusionFixed_unfixedLeftReadHead_acceptorAnchorLength: " << fusionFixed_unfixedLeftReadHead_acceptorAnchorLength << endl;
		// cout << "fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool: " << fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool << endl;

		int rawChrPos_SE_uniqAlign_clippedSeg = tmpPeAlignInfo_unfixedLeftReadHead.returnChrPos_SE_uniqAlign();
		vector<Jump_Code> rawJumpCodeVec_SE_uniqAlign_clippedSeg;
		tmpPeAlignInfo_unfixedLeftReadHead.pushBackJumpCodeVec2target_SE_uniqAlign(rawJumpCodeVec_SE_uniqAlign_clippedSeg);
		//cout << "rawChrPos_SE_uniqAlign_clippedSeg: " << rawChrPos_SE_uniqAlign_clippedSeg << endl;
		//cout << "rawJumpCodeVec_SE_uniqAlign_clippedSeg: " << this->jumpCodeVec2Str(rawJumpCodeVec_SE_uniqAlign_clippedSeg) << endl;		
		this->update_clippedSeg_jumpCodeVec_chrPos(fusionFixed_unfixedLeftReadHead_clippedSeg_chrPos,
			fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec, rawChrPos_SE_uniqAlign_clippedSeg, rawJumpCodeVec_SE_uniqAlign_clippedSeg,
			fusionFixed_unfixedLeftReadHead_donerAnchorLength, true, fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool);

		int fusionFixed_unfixedLeftReadHead_clippedSeg_seqLength = returnJumpCodeVecSeqLength(fusionFixed_unfixedLeftReadHead_clippedSeg_jumpCodeVec);
		if(fusionFixed_unfixedLeftReadHead_forMapOrRcmMap_bool)
		{	
			fusionFixed_unfixedLeftReadHead_clippedSeg_readSeq = tmpPeReadInfo_unfixedLeftReadHead.returnReadSeq_SE();
			if(fasta_or_fastq_bool)
				fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq = "*";
			else
				fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq = tmpPeReadInfo_unfixedLeftReadHead.returnQualitySeq_SE();
		}
		else
		{
			fusionFixed_unfixedLeftReadHead_clippedSeg_readSeq = tmpPeReadInfo_unfixedLeftReadHead.returnRcmReadSeq_SE();
			if(fasta_or_fastq_bool)
				fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq = "*";
			else
				fusionFixed_unfixedLeftReadHead_clippedSeg_qualSeq = tmpPeReadInfo_unfixedLeftReadHead.returnRcmQualitySeq_SE();
		}
	}

	void updateFixedFusionResult_unfixedTailAtDownstreamEnd(PE_Read_Info& tmpPeReadInfo_unfixedRightReadTail, 
		PE_Read_Alignment_Info& tmpPeAlignInfo_unfixedRightReadTail, 
		FixFusion_PeRead& tmpFixFusionPeReadInfo_unfixedRightReadTail, bool fasta_or_fastq_bool, Index_Info* indexInfo)
	{
		//cout << endl << "start to do updateFixedFusionResult_unfixedTailAtDownstreamEnd" << endl;
		fusionFixed_unfixedTailAtDownstreamEnd_success_bool = true;
		fusionFixed_unfixedRightReadTail_donerAnchorLength = tmpFixFusionPeReadInfo_unfixedRightReadTail.return_fusionFixed_donerAnchorLength();
		fusionFixed_unfixedRightReadTail_foreignInsertionLength = tmpFixFusionPeReadInfo_unfixedRightReadTail.return_fusionFixed_foreignInsertionLength();
		fusionFixed_unfixedRightReadTail_acceptorAnchorLength = tmpFixFusionPeReadInfo_unfixedRightReadTail.return_fusionFixed_acceptorAnchorLength();
		fusionFixed_unfixedRightReadTail_mismatchNum = tmpFixFusionPeReadInfo_unfixedRightReadTail.return_fusionFixed_mismatchNum();
		fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool = tmpPeAlignInfo_unfixedRightReadTail.returnForMapOrRcmMap_SE_uniqAlign();
		fusionFixed_unfixedRightReadTail_clippedSeg_chrNameInt = tmpPeAlignInfo_unfixedRightReadTail.returnChrNameInt_SE_uniqAlign(indexInfo);
		//cout << "fusionFixed_unfixedRightReadTail_donerAnchorLength: " << fusionFixed_unfixedRightReadTail_donerAnchorLength << endl;
		//cout << "fusionFixed_unfixedRightReadTail_foreignInsertionLength: " << fusionFixed_unfixedRightReadTail_foreignInsertionLength << endl;
		//cout << "fusionFixed_unfixedRightReadTail_acceptorAnchorLength: " << fusionFixed_unfixedRightReadTail_acceptorAnchorLength << endl;

		int rawChrPos_SE_uniqAlign_clippedSeg = tmpPeAlignInfo_unfixedRightReadTail.returnChrPos_SE_uniqAlign();
		vector<Jump_Code> rawJumpCodeVec_SE_uniqAlign_clippedSeg;
		tmpPeAlignInfo_unfixedRightReadTail.pushBackJumpCodeVec2target_SE_uniqAlign(rawJumpCodeVec_SE_uniqAlign_clippedSeg);
		this->update_clippedSeg_jumpCodeVec_chrPos(fusionFixed_unfixedRightReadTail_clippedSeg_chrPos,
			fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec, rawChrPos_SE_uniqAlign_clippedSeg, rawJumpCodeVec_SE_uniqAlign_clippedSeg,
			fusionFixed_unfixedRightReadTail_acceptorAnchorLength, false, fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool);

		int fusionFixed_unfixedRightReadTail_clippedSeg_seqLength = returnJumpCodeVecSeqLength(fusionFixed_unfixedRightReadTail_clippedSeg_jumpCodeVec);
		if(fusionFixed_unfixedRightReadTail_forMapOrRcmMap_bool)
		{
			fusionFixed_unfixedRightReadTail_clippedSeg_readSeq = tmpPeReadInfo_unfixedRightReadTail.returnReadSeq_SE();
			if(fasta_or_fastq_bool)
				fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq = "*";
			else
				fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq = tmpPeReadInfo_unfixedRightReadTail.returnQualitySeq_SE();
		}
		else
		{
			fusionFixed_unfixedRightReadTail_clippedSeg_readSeq = tmpPeReadInfo_unfixedRightReadTail.returnRcmReadSeq_SE();
			if(fasta_or_fastq_bool)
				fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq = "*";
			else
				fusionFixed_unfixedRightReadTail_clippedSeg_qualSeq = tmpPeReadInfo_unfixedRightReadTail.returnRcmQualitySeq_SE();
		}
	}

	void update_clippedSeg_jumpCodeVec_chrPos(int& updatedChrPos, vector<Jump_Code>& updatedJumpCodeVec, 
		int rawChrPos, vector<Jump_Code>& rawJumpCodeVec, int correspondingClippedSegAnchorSeqLength,
		bool fusionAtUpstreamEndOrDownstreamEnd_bool, bool forMapOrRcmMap_cmp2mainBody_bool)
	{
		// cout << endl << "start to do update_clippedSeg_jumpCodeVec_chrPos(" << endl;
		// cout << "fusionAtUpstreamEndOrDownstreamEnd_bool: " << fusionAtUpstreamEndOrDownstreamEnd_bool << endl;
		// cout << "forMapOrRcmMap_cmp2mainBody_bool: " << forMapOrRcmMap_cmp2mainBody_bool << endl;
		//int rawJumpCodeVecSize = rawJumpCodeVec.size();
		//for(int tmp = 0; tmp < rawJumpCodeVecSize; tmp++)
		//	updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);

		if(fusionAtUpstreamEndOrDownstreamEnd_bool) // case 2,5,10,11 
		{
			if(forMapOrRcmMap_cmp2mainBody_bool) // case 2,5 
			{	
				if(rawJumpCodeVec[rawJumpCodeVec.size()-1].type == "S")
				{
					for(int tmp = 0; tmp < rawJumpCodeVec.size()-1; tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else if(rawJumpCodeVec[rawJumpCodeVec.size()-1].type == "M")
				{
					for(int tmp = 0; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else
				{
					cout << "error rawJumpCodeVec[rawJumpCodeVec.size()-1].type != S or N" << endl;
					exit(1);
				}
				updatedChrPos = rawChrPos;
				updatedJumpCodeVec[updatedJumpCodeVec.size()-1].len = correspondingClippedSegAnchorSeqLength;
			}
			else // case 10,11
			{
				int oriEndMatchLength;
				if(rawJumpCodeVec[0].type == "S")
				{
					oriEndMatchLength = rawJumpCodeVec[1].len;
					for(int tmp = 1; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else if(rawJumpCodeVec[0].type == "M")
				{
					oriEndMatchLength = rawJumpCodeVec[0].len;
					for(int tmp = 0; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else
				{
					cout << "error rawJumpCodeVec[0].type != S or N" << endl;
					exit(1);
				}
				updatedChrPos = rawChrPos + oriEndMatchLength - correspondingClippedSegAnchorSeqLength;
				updatedJumpCodeVec[0].len = correspondingClippedSegAnchorSeqLength;
			}
		}
		else // case 1,4,7,8
		{
			if(forMapOrRcmMap_cmp2mainBody_bool) // case 1,4
			{
				int oriEndMatchLength;
				if(rawJumpCodeVec[0].type == "S")
				{
					oriEndMatchLength = rawJumpCodeVec[1].len;
					for(int tmp = 1; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else if(rawJumpCodeVec[0].type == "M")
				{
					oriEndMatchLength = rawJumpCodeVec[0].len;
					for(int tmp = 0; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else
				{
					cout << "error rawJumpCodeVec[0].type != S or N" << endl;
					exit(1);
				}
				updatedChrPos = rawChrPos + oriEndMatchLength - correspondingClippedSegAnchorSeqLength;
				updatedJumpCodeVec[0].len = correspondingClippedSegAnchorSeqLength;
			}
			else // case 7,8
			{
				if(rawJumpCodeVec[rawJumpCodeVec.size()-1].type == "S")
				{
					for(int tmp = 0; tmp < rawJumpCodeVec.size()-1; tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else if(rawJumpCodeVec[rawJumpCodeVec.size()-1].type == "M")
				{
					for(int tmp = 0; tmp < rawJumpCodeVec.size(); tmp++)
						updatedJumpCodeVec.push_back(rawJumpCodeVec[tmp]);
				}
				else
				{
					cout << "error rawJumpCodeVec[rawJumpCodeVec.size()-1].type != S or N" << endl;
					exit(1);
				}
				updatedChrPos = rawChrPos;
				updatedJumpCodeVec[updatedJumpCodeVec.size()-1].len = correspondingClippedSegAnchorSeqLength;
			}
		}
	}

	int returnJumpCodeVecSeqLength(vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpSeqLength = 0;
		for(int tmpIndex = 0; tmpIndex < cigarStringJumpCodeVec.size(); tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;		
			if((tmpJumpCodeType == "S")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "I")||(tmpJumpCodeType == "H"))
				tmpSeqLength += tmpJumpCodeLength;
			else if((tmpJumpCodeType == "D")||(tmpJumpCodeType == "N"))
			{}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}			
		}
		return tmpSeqLength;
	}

};
#endif