// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXFUSION_PEREAD_H
#define FIXFUSION_PEREAD_H

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
#include "../../../general/fixDoubleAnchorFusion_info.h"

using namespace std;

class FixFusion_PeRead
{
private:
	int chrNameInt_mainBody;
	int chrPos_mainBody;
	vector<Jump_Code> jumpCodeVec_mainBody;
	bool mainBody_forOrRevMap_bool;

	bool clippedSegAt_downstreamOrUpstream_end_bool;
	bool clippedSeg_forOrRevMap_bool;

	int chrNameInt_clippedSeg;
	int chrPos_clippedSeg;
	vector<Jump_Code> jumpCodeVec_clippedSeg;

	string inProcessTotalReadSeq;

	//results:
	bool fusionFixed_bool;
	int fusionFixed_donerAnchorLength;
	int fusionFixed_acceptorAnchorLength;
	int fusionFixed_foreignInsertionLength;
	int fusionFixed_mismatchNum;
	int fusionFixed_flankStringCase;
	bool fusionFixed_siteInAnnBoundary_bool;
public:
	FixFusion_PeRead()
	{}

	bool return_fusionFixed_bool()
	{
		return fusionFixed_bool;
	}

	int return_fusionFixed_acceptorAnchorLength()
	{
		return fusionFixed_acceptorAnchorLength;
	}

	int return_fusionFixed_donerAnchorLength()
	{
		return fusionFixed_donerAnchorLength;
	}

	int return_fusionFixed_foreignInsertionLength()
	{
		return fusionFixed_foreignInsertionLength;
	}

	int return_fusionFixed_mismatchNum()
	{
		return fusionFixed_mismatchNum;
	}

	int return_fusionFixed_flankStringCase()
	{
		return fusionFixed_flankStringCase;
	}

	bool return_fusionFixed_siteInAnnBoundary_bool()
	{
		return fusionFixed_siteInAnnBoundary_bool;
	}

	void fixFusion(PE_Read_Info& readInfo, GeneAnnEntry_Hash_Info& geneAnnEntryHashInfo, Index_Info* indexInfo, 
		bool geneAnnEntryBoundaryOnly_bool, bool geneAnnIncorporated_bool, bool insertionAllowed_bool,
		int maximum_allowed_mismatchNum, int minimum_allowed_insertionLength, int maximum_allowed_insertionLength)
	{
		string matchReadSeq_left, matchReadSeq_right, toFixReadSeq_mid, toFixChrSeq_left, toFixChrSeq_right;
		int toFixChrName_left, toFixChrName_right;
		int toFixChrSeq_left_startPos, toFixChrSeq_left_endPos;
		int toFixChrSeq_right_startPos, toFixChrSeq_right_endPos;
		bool mapStrandSameAsMainBody_left_bool, mapStrandSameAsMainBody_right_bool;
		if(clippedSegAt_downstreamOrUpstream_end_bool) // case 1,4,7,8
		{
			if(clippedSeg_forOrRevMap_bool) // case 1,4
			{
				//cout << "case 1, 4" << endl;
				int jumpCodeIndex_lastMatch_mainBody = jumpCodeVec_mainBody.size() - 2;
				if(jumpCodeVec_mainBody[jumpCodeIndex_lastMatch_mainBody].type != "M")
				{
					fusionFixed_bool = false;
					return;
				}
				int seqLength_lastMatch_mainBody = jumpCodeVec_mainBody[jumpCodeIndex_lastMatch_mainBody].len;
				int startLocInRead_lastMatch_mainBody = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec_mainBody, jumpCodeIndex_lastMatch_mainBody) - seqLength_lastMatch_mainBody + 1;
				int startLocInRead_toFixReadSeq_mid = inProcessTotalReadSeq.length() - this->returnJumpCodeVecSeqLength(jumpCodeVec_clippedSeg) + 1;
				int seqLength_toFixReadSeq_mid;
				int startLocInRead_firstMatch_clippedSeg; 
				int seqLength_firstMatch_clippedSeg;
				int jumpCodeIndex_firstMatch_clippedSeg;
				if(jumpCodeVec_clippedSeg[0].type == "S") // toFixReadSeq_mid exists
				{
					seqLength_toFixReadSeq_mid = jumpCodeVec_clippedSeg[0].len;
					startLocInRead_firstMatch_clippedSeg = startLocInRead_toFixReadSeq_mid + seqLength_toFixReadSeq_mid;
					seqLength_firstMatch_clippedSeg = jumpCodeVec_clippedSeg[1].len;
					jumpCodeIndex_firstMatch_clippedSeg = 1;
				}
				else if(jumpCodeVec_clippedSeg[0].type == "M") // toFixReadSeq_mid does not exist
				{
					seqLength_toFixReadSeq_mid = 0;
					startLocInRead_firstMatch_clippedSeg = startLocInRead_toFixReadSeq_mid;
					seqLength_firstMatch_clippedSeg = jumpCodeVec_clippedSeg[0].len;
					jumpCodeIndex_firstMatch_clippedSeg = 0;
				}
				else
				{
					cout << "invalid jumpCodeVec_clippedSeg[0].type: " << jumpCodeVec_clippedSeg[0].type << endl;
					exit(1);
				}

				matchReadSeq_left = inProcessTotalReadSeq.substr(startLocInRead_lastMatch_mainBody - 1, seqLength_lastMatch_mainBody); 
				matchReadSeq_right = inProcessTotalReadSeq.substr(startLocInRead_firstMatch_clippedSeg - 1, seqLength_firstMatch_clippedSeg);
				//string toFixReadSeq_mid;
				if(seqLength_toFixReadSeq_mid == 0)
					toFixReadSeq_mid = "";
				else
					toFixReadSeq_mid = inProcessTotalReadSeq.substr(startLocInRead_toFixReadSeq_mid - 1, seqLength_toFixReadSeq_mid);
				
				int toFixDoubleAnchorFusionReadSubSeqLength = seqLength_lastMatch_mainBody + seqLength_toFixReadSeq_mid + seqLength_firstMatch_clippedSeg;
				int startChrPos_lastMatch_mainBody = this->getEndPosOfSpecificJumpCode(chrPos_mainBody, jumpCodeVec_mainBody, jumpCodeIndex_lastMatch_mainBody) - seqLength_lastMatch_mainBody + 1;
				toFixChrSeq_left = indexInfo->returnChromStrSubstr(chrNameInt_mainBody, startChrPos_lastMatch_mainBody, toFixDoubleAnchorFusionReadSubSeqLength);
				int endChrPos_firstMatch_clippedSeg = this->getEndPosOfSpecificJumpCode(chrPos_clippedSeg, jumpCodeVec_clippedSeg, jumpCodeIndex_firstMatch_clippedSeg);
				toFixChrSeq_right = indexInfo->returnChromStrSubstr(chrNameInt_clippedSeg, endChrPos_firstMatch_clippedSeg - toFixDoubleAnchorFusionReadSubSeqLength + 1, toFixDoubleAnchorFusionReadSubSeqLength);
				toFixChrSeq_left_startPos = startChrPos_lastMatch_mainBody;
				toFixChrSeq_left_endPos = startChrPos_lastMatch_mainBody + toFixDoubleAnchorFusionReadSubSeqLength - 1;
				toFixChrSeq_right_startPos = endChrPos_firstMatch_clippedSeg - toFixDoubleAnchorFusionReadSubSeqLength + 1;
				toFixChrSeq_right_endPos = endChrPos_firstMatch_clippedSeg;
				mapStrandSameAsMainBody_left_bool = true;
				mapStrandSameAsMainBody_right_bool = true;
				toFixChrName_left = chrNameInt_mainBody; 
				toFixChrName_right = chrNameInt_clippedSeg;
			}
			else // case 7,8
			{
				//cout << "case 7, 8" << endl;
				//cout << "jumpCodeVec_mainBody: ";
				// for(int tmp = 0; tmp < jumpCodeVec_mainBody.size(); tmp++)
				// 	cout << jumpCodeVec_mainBody[tmp].len << jumpCodeVec_mainBody[tmp].type;
				// cout << endl;
				// cout << "jumpCodeVec_clippedSeg: ";
				// for(int tmp = 0; tmp < jumpCodeVec_clippedSeg.size(); tmp++)
				// 	cout << jumpCodeVec_clippedSeg[tmp].len << jumpCodeVec_clippedSeg[tmp].type;
				// cout << endl;				

				int seqLength_lastMatch_mainBody = jumpCodeVec_mainBody[jumpCodeVec_mainBody.size() - 2].len;
				int startLocInRead_lastMatch_mainBody = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec_mainBody, jumpCodeVec_mainBody.size() - 2) - seqLength_lastMatch_mainBody + 1;
				int seqLength_lastMatch_clippedSeg;
				int seqLength_toFixReadSeq_mid;
				int endLocInRead_toFixReadSeq_mid_rev;
				int endLocInRead_lastMatch_clippedSeg_rev;
				int startChrPos_lastMatch_clippedSeg_rev;
				//cout << "jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type: " << jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type << endl;
				if(jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type == "S")
				{
 					seqLength_lastMatch_clippedSeg = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 2].len;
					seqLength_toFixReadSeq_mid = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].len;
					endLocInRead_toFixReadSeq_mid_rev = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec_mainBody, jumpCodeVec_mainBody.size() - 2) + 1;
					endLocInRead_lastMatch_clippedSeg_rev = endLocInRead_toFixReadSeq_mid_rev + seqLength_toFixReadSeq_mid;			
					startChrPos_lastMatch_clippedSeg_rev = this->getEndPosOfSpecificJumpCode(chrPos_clippedSeg, jumpCodeVec_clippedSeg, jumpCodeVec_clippedSeg.size() - 2) - seqLength_lastMatch_clippedSeg + 1;
				}
				else if(jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type == "M")
				{
 					seqLength_lastMatch_clippedSeg = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].len;
					seqLength_toFixReadSeq_mid = 0;
					//endLocInRead_toFixReadSeq_mid_rev = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec_mainBody, jumpCodeVec_mainBody.size() - 2) + 1;
					endLocInRead_lastMatch_clippedSeg_rev = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec_mainBody, jumpCodeVec_mainBody.size() - 2) + 1;
					startChrPos_lastMatch_clippedSeg_rev = this->getEndPosOfSpecificJumpCode(chrPos_clippedSeg, jumpCodeVec_clippedSeg, jumpCodeVec_clippedSeg.size() - 1) - seqLength_lastMatch_clippedSeg + 1; 
				}
				else
				{
					cout << "invalid jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type: " << jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type << endl;
					exit(1);
				}
				//cout << "endLocInRead_lastMatch_clippedSeg_rev: " << endLocInRead_lastMatch_clippedSeg_rev << endl;
				//cout << "seqLength_lastMatch_clippedSeg: " << seqLength_lastMatch_clippedSeg << endl;
				matchReadSeq_left = inProcessTotalReadSeq.substr(startLocInRead_lastMatch_mainBody - 1, seqLength_lastMatch_mainBody); 
				matchReadSeq_right = inProcessTotalReadSeq.substr(endLocInRead_lastMatch_clippedSeg_rev - 1, seqLength_lastMatch_clippedSeg); 
				//string toFixReadSeq_mid;
				if(seqLength_toFixReadSeq_mid == 0)
					toFixReadSeq_mid = "";
				else
					toFixReadSeq_mid = inProcessTotalReadSeq.substr(startLocInRead_lastMatch_mainBody + seqLength_lastMatch_mainBody, seqLength_toFixReadSeq_mid);
				int toFixDoubleAnchorFusionReadSubSeqLength = seqLength_lastMatch_mainBody + seqLength_toFixReadSeq_mid + seqLength_lastMatch_clippedSeg;
				int startChrPos_lastMatch_mainBody = this->getEndPosOfSpecificJumpCode(chrPos_mainBody, jumpCodeVec_mainBody, jumpCodeVec_mainBody.size() - 2) - seqLength_lastMatch_mainBody + 1;
				toFixChrSeq_left = indexInfo->returnChromStrSubstr(chrNameInt_mainBody, startChrPos_lastMatch_mainBody, toFixDoubleAnchorFusionReadSubSeqLength);
				toFixChrSeq_right = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_clippedSeg, startChrPos_lastMatch_clippedSeg_rev, toFixDoubleAnchorFusionReadSubSeqLength));
				toFixChrSeq_left_startPos = startChrPos_lastMatch_mainBody; 
				toFixChrSeq_left_endPos = startChrPos_lastMatch_mainBody + toFixDoubleAnchorFusionReadSubSeqLength - 1;
				toFixChrSeq_right_startPos = startChrPos_lastMatch_clippedSeg_rev;
				toFixChrSeq_right_endPos = startChrPos_lastMatch_clippedSeg_rev + toFixDoubleAnchorFusionReadSubSeqLength - 1;
				mapStrandSameAsMainBody_left_bool = true;
				mapStrandSameAsMainBody_right_bool = false;
				toFixChrName_left = chrNameInt_mainBody; 
				toFixChrName_right = chrNameInt_clippedSeg;
			}
		}
		else // case 2,5,10,11
		{
			if(clippedSeg_forOrRevMap_bool) // case 2,5
			{
				//cout << "case 2, 5" << endl;
				int jumpCodeIndex_firstMatch_mainBody = 1;
				if(jumpCodeVec_mainBody[1].type != "M")
				{	
					fusionFixed_bool = false;
					return;
				}
				int seqLength_firstMatch_mainBody = jumpCodeVec_mainBody[1].len;
				int startLocInRead_firstMatch_mainBody = jumpCodeVec_mainBody[0].len + 1;
				//cout << "seqLength_firstMatch_mainBody: " << seqLength_firstMatch_mainBody << endl;
				//cout << "startLocInRead_firstMatch_mainBody: " << startLocInRead_firstMatch_mainBody << endl;
				int startLocInRead_toFixReadSeq_mid;
				int seqLength_toFixReadSeq_mid;
				int startLocInRead_lastMatch_clippedSeg;
				int seqLength_lastMatch_clippedSeg;
				int startChrPos_lastMatch_clippedSeg;
				if(jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type == "S")
				{
					seqLength_toFixReadSeq_mid = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].len;
					startLocInRead_toFixReadSeq_mid = startLocInRead_firstMatch_mainBody - seqLength_toFixReadSeq_mid;
					seqLength_lastMatch_clippedSeg = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 2].len;
					startLocInRead_lastMatch_clippedSeg = startLocInRead_toFixReadSeq_mid - seqLength_lastMatch_clippedSeg;
					startChrPos_lastMatch_clippedSeg = this->getEndPosOfSpecificJumpCode(chrPos_clippedSeg, jumpCodeVec_clippedSeg, jumpCodeVec_clippedSeg.size() - 2) - seqLength_lastMatch_clippedSeg + 1;
				}
				else if(jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type == "M")
				{
					seqLength_toFixReadSeq_mid = 0;
					startLocInRead_toFixReadSeq_mid = - 1;
					seqLength_lastMatch_clippedSeg = jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].len;
					startLocInRead_lastMatch_clippedSeg = startLocInRead_firstMatch_mainBody - seqLength_lastMatch_clippedSeg;
					startChrPos_lastMatch_clippedSeg = this->getEndPosOfSpecificJumpCode(chrPos_clippedSeg, jumpCodeVec_clippedSeg, jumpCodeVec_clippedSeg.size() - 1) - seqLength_lastMatch_clippedSeg + 1;
				}
				else
				{
					cout << "invalid jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type: " << jumpCodeVec_clippedSeg[jumpCodeVec_clippedSeg.size() - 1].type << endl;
					exit(1);
				}
				// cout << "startLocInRead_toFixReadSeq_mid: " << startLocInRead_toFixReadSeq_mid << endl;
				// cout << "seqLength_toFixReadSeq_mid: " << seqLength_toFixReadSeq_mid << endl;
				// cout << "startLocInRead_lastMatch_clippedSeg: " << startLocInRead_lastMatch_clippedSeg << endl;
				// cout << "seqLength_lastMatch_clippedSeg: " << seqLength_lastMatch_clippedSeg << endl;
				// cout << "startChrPos_lastMatch_clippedSeg: " << startChrPos_lastMatch_clippedSeg << endl;
				matchReadSeq_left = inProcessTotalReadSeq.substr(startLocInRead_lastMatch_clippedSeg - 1, seqLength_lastMatch_clippedSeg); 
				matchReadSeq_right = inProcessTotalReadSeq.substr(startLocInRead_firstMatch_mainBody - 1, seqLength_firstMatch_mainBody);
				//cout << "matchReadSeq_left: " << matchReadSeq_left << endl;
				//cout << "matchReadSeq_right: " << matchReadSeq_right << endl;
				//toFixReadSeq_mid;
				if(seqLength_toFixReadSeq_mid == 0)
					toFixReadSeq_mid = "";
				else
					toFixReadSeq_mid = inProcessTotalReadSeq.substr(startLocInRead_toFixReadSeq_mid - 1, seqLength_toFixReadSeq_mid);
				int toFixDoubleAnchorFusionReadSubSeqLength = seqLength_firstMatch_mainBody + seqLength_toFixReadSeq_mid + seqLength_lastMatch_clippedSeg;
				//cout << "toFixDoubleAnchorFusionReadSubSeqLength: " << toFixDoubleAnchorFusionReadSubSeqLength << endl;
				toFixChrSeq_left = indexInfo->returnChromStrSubstr(chrNameInt_clippedSeg, startChrPos_lastMatch_clippedSeg, toFixDoubleAnchorFusionReadSubSeqLength);
				//cout << "toFixChrSeq_left: " << toFixChrSeq_left << endl;
				int endChrPos_firstMatch_mainBody = this->getEndPosOfSpecificJumpCode(chrPos_mainBody, jumpCodeVec_mainBody, 1);
				//cout << "endChrPos_firstMatch_mainBody: " << endChrPos_firstMatch_mainBody << endl;
				toFixChrSeq_right = indexInfo->returnChromStrSubstr(chrNameInt_mainBody, endChrPos_firstMatch_mainBody - toFixDoubleAnchorFusionReadSubSeqLength + 1, toFixDoubleAnchorFusionReadSubSeqLength);
				//cout << "toFixChrSeq_right: " << toFixChrSeq_right << endl;
				toFixChrSeq_left_startPos = startChrPos_lastMatch_clippedSeg; 
				toFixChrSeq_left_endPos = startChrPos_lastMatch_clippedSeg + toFixDoubleAnchorFusionReadSubSeqLength - 1;
				toFixChrSeq_right_startPos = endChrPos_firstMatch_mainBody - toFixDoubleAnchorFusionReadSubSeqLength + 1;
				toFixChrSeq_right_endPos = endChrPos_firstMatch_mainBody;
				mapStrandSameAsMainBody_left_bool = true;
				mapStrandSameAsMainBody_right_bool = true;
				toFixChrName_left = chrNameInt_clippedSeg; 
				toFixChrName_right = chrNameInt_mainBody;
			}
			else // case 10,11
			{
				//cout << "case 10, 11" << endl;
				int startLocInRead_firstMatch_mainBody = jumpCodeVec_mainBody[0].len + 1;
				int seqLength_firstMatch_mainBody = jumpCodeVec_mainBody[1].len; 
				int seqLength_firstMatch_clippedSeg;
				int seqLength_toFixReadSeq_mid;
				int endLocInRead_toFixReadSeq_mid_rev;
				int endLocInRead_firstMatch_clippedSeg_rev;
				if(jumpCodeVec_clippedSeg[0].type == "S")
				{
					seqLength_firstMatch_clippedSeg = jumpCodeVec_clippedSeg[1].len;
					seqLength_toFixReadSeq_mid = jumpCodeVec_clippedSeg[0].len;
					endLocInRead_toFixReadSeq_mid_rev = startLocInRead_firstMatch_mainBody - seqLength_toFixReadSeq_mid;
					endLocInRead_firstMatch_clippedSeg_rev = endLocInRead_toFixReadSeq_mid_rev - seqLength_firstMatch_clippedSeg;
				}
				else if(jumpCodeVec_clippedSeg[0].type == "M")
				{
					seqLength_firstMatch_clippedSeg = jumpCodeVec_clippedSeg[0].len;
					seqLength_toFixReadSeq_mid = 0;
					//endLocInRead_toFixReadSeq_mid_rev = startLocInRead_firstMatch_mainBody - 1;
					endLocInRead_firstMatch_clippedSeg_rev = startLocInRead_firstMatch_mainBody - seqLength_firstMatch_clippedSeg;
				}
				else
				{
					cout << "invalid jumpCodeVec_clippedSeg[0].type: " << jumpCodeVec_clippedSeg[0].type << endl;
					exit(1);
				}

				matchReadSeq_left = inProcessTotalReadSeq.substr(endLocInRead_firstMatch_clippedSeg_rev - 1, seqLength_firstMatch_clippedSeg);
				matchReadSeq_right = inProcessTotalReadSeq.substr(startLocInRead_firstMatch_mainBody - 1, seqLength_firstMatch_mainBody);
				//string toFixReadSeq_mid;
				if(seqLength_toFixReadSeq_mid == 0)
					toFixReadSeq_mid = "";
				else
					toFixReadSeq_mid = inProcessTotalReadSeq.substr(endLocInRead_toFixReadSeq_mid_rev - 1, seqLength_toFixReadSeq_mid);
				int toFixDoubleAnchorFusionReadSubSeqLength = seqLength_firstMatch_mainBody + seqLength_toFixReadSeq_mid + seqLength_firstMatch_clippedSeg;
				int endChrPos_firstMatch_clippedSeg = chrPos_clippedSeg + seqLength_firstMatch_clippedSeg - 1;
				toFixChrSeq_left = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_clippedSeg, 
					endChrPos_firstMatch_clippedSeg - toFixDoubleAnchorFusionReadSubSeqLength + 1, toFixDoubleAnchorFusionReadSubSeqLength)); 
				toFixChrSeq_right = indexInfo->returnChromStrSubstr(chrNameInt_mainBody, chrPos_mainBody + seqLength_firstMatch_mainBody - 1 
					- toFixDoubleAnchorFusionReadSubSeqLength + 1, toFixDoubleAnchorFusionReadSubSeqLength);
				toFixChrSeq_left_startPos = endChrPos_firstMatch_clippedSeg - toFixDoubleAnchorFusionReadSubSeqLength + 1;
				toFixChrSeq_left_endPos = endChrPos_firstMatch_clippedSeg;
				toFixChrSeq_right_startPos = chrPos_mainBody + seqLength_firstMatch_mainBody - 1 - toFixDoubleAnchorFusionReadSubSeqLength + 1;
				toFixChrSeq_right_endPos = chrPos_mainBody + seqLength_firstMatch_mainBody - 1;
				mapStrandSameAsMainBody_left_bool = false;
				mapStrandSameAsMainBody_right_bool = true;
				toFixChrName_left = chrNameInt_clippedSeg;
				toFixChrName_right = chrNameInt_mainBody;
			}
		}
		// cout << "start to do fixDoubleAnchorFusion_info ...." << endl;
		// cout << "matchReadSeq_left.length(): " << matchReadSeq_left.length() << endl;
		// cout << "matchReadSeq_right.length(): " << matchReadSeq_right.length() << endl;
		// cout << "toFixReadSeq_mid.length(): " << toFixReadSeq_mid.length() << endl;
		// cout << "toFixChrName_left: " << toFixChrName_left << endl;
		// cout << "toFixChrName_right: " << toFixChrName_right << endl;
		// cout << "toFixChrSeq_left: " << toFixChrSeq_left_startPos << " ~ " << toFixChrSeq_left_endPos << endl;
		// cout << "toFixChrSeq_right: " << toFixChrSeq_right_startPos << " ~ " << toFixChrSeq_right_endPos << endl;
		// cout << " mapStrandSameAsMainBody_left_bool: " << mapStrandSameAsMainBody_left_bool << endl;
		// cout << " mapStrandSameAsMainBody_right_bool: " << mapStrandSameAsMainBody_right_bool << endl;
		FixDoubleAnchor_Fusion_Info tmpFixFusionInfo;
		tmpFixFusionInfo.fixFusion(matchReadSeq_left, matchReadSeq_right, toFixReadSeq_mid, toFixChrSeq_left, toFixChrSeq_right,
			toFixChrName_left, toFixChrSeq_left_startPos, toFixChrSeq_left_endPos, mapStrandSameAsMainBody_left_bool,
			toFixChrName_right, toFixChrSeq_right_startPos, toFixChrSeq_right_endPos, mapStrandSameAsMainBody_right_bool,
			geneAnnEntryHashInfo, indexInfo, geneAnnEntryBoundaryOnly_bool, geneAnnIncorporated_bool, insertionAllowed_bool,
			maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength);
		//cout << "end of doing fixDoubleAnchorFusion_info ..." << endl;
		fusionFixed_bool = tmpFixFusionInfo.return_fusionFixed_bool();
		fusionFixed_donerAnchorLength = tmpFixFusionInfo.return_theBestFusionSite_donerLength();
		fusionFixed_acceptorAnchorLength = tmpFixFusionInfo.return_theBestFusionSite_acceptorLength();
		fusionFixed_foreignInsertionLength = tmpFixFusionInfo.return_theBestFusionSite_insLength();
		fusionFixed_mismatchNum = tmpFixFusionInfo.return_theBestFusionSite_mismatchNum();
		fusionFixed_flankStringCase = tmpFixFusionInfo.return_theBestFusionSite_flankStringCase();
		fusionFixed_siteInAnnBoundary_bool = tmpFixFusionInfo.return_theBestFusionSite_annotated_bool();
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
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	void initiate_clippedSegGlobalMap(PeSam_Info& peSamInfo, PE_Read_Alignment_Info& peAlignInfo,
		bool clippedSeg_unfixedHeadAtUpstreamEnd_or_unfixedTailAtDownstreamEnd_bool, PE_Read_Info& tmpPeReadInfo, Index_Info* indexInfo)
	{
		bool tmp_mainBody_forOrRevMap_bool = peSamInfo.return_end1forMapEnd2RevMap_end2forMapEnd1RevMap_bool();
		bool tmp_clippedSegAt_downstreamOrUpstream_end_bool = (!clippedSeg_unfixedHeadAtUpstreamEnd_or_unfixedTailAtDownstreamEnd_bool);
		bool tmp_clippedSeg_forOrRevMap_bool;
		
		int tmp_chrNameInt_mainBody = peSamInfo.returnChrNameInt();

		int tmp_chrPos_mainBody;
		if(tmp_clippedSegAt_downstreamOrUpstream_end_bool)
			tmp_chrPos_mainBody = peSamInfo.return_startPos_downstreamEnd();
		else
			tmp_chrPos_mainBody = peSamInfo.return_startPos_upstreamEnd();

		vector<Jump_Code> tmp_jumpCodeVec_mainBody;
		if(tmp_clippedSegAt_downstreamOrUpstream_end_bool)
			peSamInfo.pushJumpCodeVec2target_downstreamEnd(tmp_jumpCodeVec_mainBody);
		else
			peSamInfo.pushJumpCodeVec2target_upstreamEnd(tmp_jumpCodeVec_mainBody);

		int tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg;
		vector<Jump_Code> tmp_jumpCodeVec_clippedSeg;
		peAlignInfo.initiate_clippedSegSamForFusionPostAnalysis_SE_uniqAlign(tmp_clippedSeg_forOrRevMap_bool, 
			tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg, tmp_jumpCodeVec_clippedSeg, indexInfo);

		//PE_Read_Info tmpPeReadInfo;

		this->initiate(tmp_chrNameInt_mainBody, tmp_chrPos_mainBody, tmp_jumpCodeVec_mainBody, 
				tmp_mainBody_forOrRevMap_bool, tmp_clippedSegAt_downstreamOrUpstream_end_bool, tmp_clippedSeg_forOrRevMap_bool, 
				tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg, tmp_jumpCodeVec_clippedSeg, tmpPeReadInfo);
	}

	void initiate(int tmp_chrNameInt_mainBody, int tmp_chrPos_mainBody, vector<Jump_Code>& tmp_jumpCodeVec_mainBody,
		bool tmp_mainBody_forOrRevMap_bool, bool tmp_clippedSegAt_downstreamOrUpstream_end_bool, bool tmp_clippedSeg_forOrRevMap_bool,
		int tmp_chrNameInt_clippedSeg, int tmp_chrPos_clippedSeg, vector<Jump_Code>& tmp_jumpCodeVec_clippedSeg, PE_Read_Info& peReadInfo)
	{
		chrNameInt_mainBody = tmp_chrNameInt_mainBody;
		chrPos_mainBody = tmp_chrPos_mainBody;
		//cout << "chrNameInt_mainBody: " << chrNameInt_mainBody << endl;
		//cout << "chrPos_mainBody: " << chrPos_mainBody << endl;
		for(int tmp = 0; tmp < tmp_jumpCodeVec_mainBody.size(); tmp++)
			jumpCodeVec_mainBody.push_back(tmp_jumpCodeVec_mainBody[tmp]);

		mainBody_forOrRevMap_bool = tmp_mainBody_forOrRevMap_bool;
		clippedSegAt_downstreamOrUpstream_end_bool = tmp_clippedSegAt_downstreamOrUpstream_end_bool;
		clippedSeg_forOrRevMap_bool = tmp_clippedSeg_forOrRevMap_bool;
		//cout << "mainBody_forOrRevMap_bool: " << mainBody_forOrRevMap_bool << endl;
		//cout << "clippedSegAt_downstreamOrUpstream_end_bool: " << clippedSegAt_downstreamOrUpstream_end_bool << endl;
		//cout << "clippedSeg_forOrRevMap_bool: " << clippedSeg_forOrRevMap_bool << endl;

		chrNameInt_clippedSeg = tmp_chrNameInt_clippedSeg;
		chrPos_clippedSeg = tmp_chrPos_clippedSeg;
		//cout << "chrNameInt_clippedSeg: " << chrNameInt_clippedSeg << endl;
		//cout << "chrPos_clippedSeg: " << chrPos_clippedSeg << endl;

		for(int tmp = 0; tmp < tmp_jumpCodeVec_clippedSeg.size(); tmp++)
			jumpCodeVec_clippedSeg.push_back(tmp_jumpCodeVec_clippedSeg[tmp]);

		if(mainBody_forOrRevMap_bool) // upstreamEnd -- end1-for, downstreamEnd -- end2-rev 
		{
			if(clippedSegAt_downstreamOrUpstream_end_bool) // fusion in downstream
				inProcessTotalReadSeq = peReadInfo.returnRcmReadSeq_2();
			else  // fusion in uptream
				inProcessTotalReadSeq = peReadInfo.returnReadSeq_1();
		}
		else // upstreamEnd -- end2-for, downstreamEnd -- end1 -- rev
		{
			if(clippedSegAt_downstreamOrUpstream_end_bool) // fusion in downstream
				inProcessTotalReadSeq = peReadInfo.returnRcmReadSeq_1();
			else  // fusion in uptream
				inProcessTotalReadSeq = peReadInfo.returnReadSeq_2();
		}
		//cout << "inProcessTotalReadSeq: " << endl << inProcessTotalReadSeq << endl;
	}
};

#endif