// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORFUSION_INFO_H
#define FIXDOUBLEANCHORFUSION_INFO_H

//#include "fixDoubleAnchorSplice_complicate_info.h"
#include "nw_DP.h"
#include "../otherProjects/dreamChallenge/general/geneAnnEntryHash.h"

#define CANONICAL_FUSION_PENALTY 0
#define SEMICANONICAL_FUSION_PENALTY 2
#define NONCANONICAL_FUSION_PENALTY 3
#define PER_MISMATCH_PENALTY 1
#define PER_INSERTION_PENALTY 1
#define ANNOTATED_SCORE 5

#define ALLOWED_DELETION_IN_READ_DUETO_DREAMCHALLENGE_SIMULATIONERROR 0 // original valueis 1

using namespace std;

class FixDoubleAnchor_Fusion_Info
{
private:
	//int chrNameInt_doner, chrNameInt_acceptor;
	vector<int> donerMapCumulativeMismatchNumVec; 
	// donerMapCumulativeMismatchNumVec[tmp] = Y;  Y is mismatch # if readSeq.substr(toFixSeqLocInRead_start-1, tmp+1) 
	// or toFixSeq.substr(0, tmp+1)  mapped to chromSeq.substr(toFixSeqMapPos_doner - 1, tmp+1); 
	vector<int> donerMapMismatchPosVec;
	//vector<char> donerMapMismatchCharVec;

	vector<int> acceptorMapCumulativeMismatchNumVec;
	// acceptorMapCumulativeMismatchNumVec[tmp] = Y; Y is mismatch # if readSeq.substr(toFixSeqLocInRead_end - 1 - tmp, tmp+1) 
	// or toFixSeq.substr(toFixSeqLength - 1- tmp, tmp+1) mapped to chromSeq.substr(toFixSeqMapPos_end - tmp - 1, tmp+1);
	vector<int> acceptorMapMismatchPosVec;
	//vector<char> acceptorMapMismatchCharVec;

	//int max_extension_forward_doner;
	//int max_extension_backward_acceptor;

	// fusion site, anchorSeqLengthVec, insLengthVec, mismatchNumVec, signalScoreVec
	vector< pair<int,int> > fusionSiteAnchorSeqLengthPairVec;
	vector<int> fusionSiteInsLengthVec;
	vector<int> fusionSiteMismatchNumVec;
	vector<int> fusionSiteFlankStringCaseVec;

	vector<bool> fusionSiteAnnotatedBoolVec;

	vector<int> annotatedIndexVecInFusionSiteVec;

	// RESULTS:
	bool fusionFixed_bool;
	int theBestFusionSite_donerLength;
	int theBestFusionSite_acceptorLength;
	int theBestFusionSite_insLength;
	int theBestFusionSite_mismatchNum;
	int theBestFusionSite_flankStringCase;
	bool theBestFusionSite_annotated_bool;
public:
	FixDoubleAnchor_Fusion_Info()
	{
		fusionFixed_bool = false;
	}

	bool return_fusionFixed_bool()
	{
		return fusionFixed_bool;
	}

	int return_theBestFusionSite_donerLength()
	{
		return theBestFusionSite_donerLength;
	}

	int return_theBestFusionSite_acceptorLength()
	{
		return theBestFusionSite_acceptorLength;
	}

	int return_theBestFusionSite_insLength()
	{
		return theBestFusionSite_insLength;
	}

	int return_theBestFusionSite_mismatchNum()
	{
		return theBestFusionSite_mismatchNum;
	}

	int return_theBestFusionSite_flankStringCase()
	{
		return theBestFusionSite_flankStringCase;
	}

	int return_theBestFusionSite_annotated_bool()
	{
		return theBestFusionSite_annotated_bool;
	}

	void fixFusion(string& matchReadSeq_left, string& matchReadSeq_right, string& toFixReadSeq_mid, 
		string& toFixChrSeq_left, string& toFixChrSeq_right,
		int chrNameInt_1, int toFixChrSeq_1_startPos, int toFixChrSeq_1_endPos, bool mapStrandSameAsMainBody_1_bool,
		int chrNameInt_2, int toFixChrSeq_2_startPos, int toFixChrSeq_2_endPos, bool mapStrandSameAsMainBody_2_bool,		
		//int toFixReadSeq_1_match_startLocInRead, int toFixReadSeq_1_match_endLocInRead,
		//int toFixReadSeq_2_match_startLocInRead, int toFixReadSeq_2_match_endLocInRead,
		GeneAnnEntry_Hash_Info& geneAnnEntryHashInfo, Index_Info* indexInfo,//int toFixSeqLength,
		bool geneAnnEntryBoundaryOnly_bool, bool geneAnnIncorporated_bool, bool insertionAllowed_bool,
		int maximum_allowed_mismatchNum, int minimum_allowed_insertionLength, int maximum_allowed_insertionLength)
	{
		//cout << "fix fusion starts ..." << endl;
		int matchReadSeq_left_length = matchReadSeq_left.length();
		int matchReadSeq_right_length = matchReadSeq_right.length();
		int toFixReadSeq_mid_length = toFixReadSeq_mid.length();
		// always, toFixChrSeq_left_length == toFixChrSeq_right_length == toFixSeqLength
		int toFixChrSeq_left_length = toFixChrSeq_left.length();
		int toFixChrSeq_right_length = toFixChrSeq_right.length();
		int toFixSeqLength = toFixChrSeq_left_length;
		// cout << "matchReadSeq_left_length: " << matchReadSeq_left_length << endl;
		// cout << "matchReadSeq_right_length: " << matchReadSeq_right_length << endl;
		// cout << "toFixReadSeq_mid_length: " << toFixReadSeq_mid_length << endl;
		// cout << "toFixChrSeq_left_length: " << toFixChrSeq_left_length << endl;
		// cout << "toFixChrSeq_right_length: " << toFixChrSeq_right_length << endl;
		// cout << "toFixSeqLength: " << toFixSeqLength << endl;

		// examine parameter settings
		// always, toFixChrSeq_left_length == toFixChrSeq_right_length == toFixSeqLength
		if((toFixChrSeq_left_length != toFixChrSeq_right_length)||(matchReadSeq_left_length 
			+ matchReadSeq_right_length + toFixReadSeq_mid_length != toFixChrSeq_left_length))
		{
			cout << "(toFixChrSeq_left_length != toFixChrSeq_right_length)||(matchReadSeq_left_length + matchReadSeq_right_length + toFixReadSeq_mid_length != toFixChrSeq_left_length)" << endl;
			cout << "matchReadSeq_left_length: " << matchReadSeq_left_length << endl;
			cout << "matchReadSeq_right_length: " << matchReadSeq_right_length << endl;
			cout << "toFixReadSeq_mid_length: " << toFixReadSeq_mid_length << endl;
			cout << "toFixChrSeq_left_length: " << toFixChrSeq_left_length << endl;
			cout << "toFixChrSeq_right_length: " << toFixChrSeq_right_length << endl;
			exit(1);
		}

		string toFixReadSeq = matchReadSeq_left + toFixReadSeq_mid + matchReadSeq_right;
		//cout << "toFixReadSeq: " << toFixReadSeq << endl;
		this->scanGenomeAndReadSeq(toFixReadSeq, toFixChrSeq_left, toFixChrSeq_right);
		//cout << "ends of scanning genome and read seq" << endl;
		if(!geneAnnEntryBoundaryOnly_bool)
			this->generateFusionSiteVec_denovo(toFixReadSeq, toFixChrSeq_left, toFixChrSeq_right, indexInfo, 
				insertionAllowed_bool, maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength);
		//cout << "ends of generating fusion site vec - denovo" << endl;
		if(geneAnnIncorporated_bool)
			this->addNewFusionSite_generateAnnotatedIndexVec_withGeneAnn(
				toFixChrSeq_left, toFixChrSeq_right, geneAnnEntryHashInfo, toFixSeqLength,
				chrNameInt_1, toFixChrSeq_1_startPos, toFixChrSeq_1_endPos, mapStrandSameAsMainBody_1_bool,
				chrNameInt_2, toFixChrSeq_2_startPos, toFixChrSeq_2_endPos, mapStrandSameAsMainBody_2_bool,
				maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength);
		//cout << "ends of adding new fusion sites..." << endl;
		this->generateFusionSiteAnnotatedBoolVec(); // even no annotated is provided, pushBack false for each fusion site
		//cout << "ends of generating fusion site annotated bool vec " << endl;
		int theBestFusionSiteIndex = selectBestFusionSiteIndex(maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength); // best splice site, return index in fusionSiteAnchorSeqLengthPairVec
		//cout << "ends of selecting best fusion site index" << endl;
		if(theBestFusionSiteIndex < 0) // no fusion site is detected, return false;
		{
			fusionFixed_bool = false;
			return;
		}
		fusionFixed_bool = true;
		theBestFusionSite_donerLength = fusionSiteAnchorSeqLengthPairVec[theBestFusionSiteIndex].first;
		theBestFusionSite_acceptorLength = fusionSiteAnchorSeqLengthPairVec[theBestFusionSiteIndex].second;
		theBestFusionSite_insLength = fusionSiteInsLengthVec[theBestFusionSiteIndex];
		theBestFusionSite_mismatchNum = fusionSiteMismatchNumVec[theBestFusionSiteIndex];
		theBestFusionSite_flankStringCase = fusionSiteFlankStringCaseVec[theBestFusionSiteIndex];
		theBestFusionSite_annotated_bool = fusionSiteAnnotatedBoolVec[theBestFusionSiteIndex];
		// cout << "theBestFusionSite_donerLength: " << theBestFusionSite_donerLength << endl;
		// cout << "theBestFusionSite_acceptorLength: " << theBestFusionSite_acceptorLength << endl;
		// cout << "theBestFusionSite_insLength: " << theBestFusionSite_insLength << endl;
		// cout << "theBestFusionSite_mismatchNum: " << theBestFusionSite_mismatchNum << endl;
		// cout << "theBestFusionSite_flankStringCase: " << theBestFusionSite_flankStringCase << endl;
		// cout << "theBestFusionSite_annotated_bool: " << theBestFusionSite_annotated_bool << endl;
		return;
	}

	void scanGenomeAndReadSeq(string& toFixReadSeq, string& toFixChrSeq_left, string& toFixChrSeq_right)
	{
		int toFixSeqLength = toFixReadSeq.length();
		//max_extension_forward_doner = toFixSeqLength;
		//max_extension_backward_acceptor = toFixSeqLength;

		int tmpDonerMapCumulativeMismatchNum = 0;
		donerMapCumulativeMismatchNumVec.push_back(0);
		for(int tmp = 0; tmp < toFixSeqLength; tmp++)
		{
			char tmpCharInRead = toFixReadSeq.at(tmp);
			char tmpCharInChr = toFixChrSeq_left.at(tmp);
			if(tmpCharInRead != tmpCharInChr)
			{
				tmpDonerMapCumulativeMismatchNum ++;
				donerMapMismatchPosVec.push_back(tmp + 1);
				//donerMapMismatchCharVec.push_back(tmpCharInRead);
			}
			donerMapCumulativeMismatchNumVec.push_back(tmpDonerMapCumulativeMismatchNum);
		}

		int tmpAcceptorMapCumulativeMismatchNum = 0;
		acceptorMapCumulativeMismatchNumVec.push_back(0);
		for(int tmp = 0; tmp < toFixSeqLength; tmp++)
		{
			char tmpCharInRead = toFixReadSeq.at(toFixSeqLength - 1 - tmp);
			char tmpCharInChr = toFixChrSeq_right.at(toFixSeqLength - 1 - tmp);
			if(tmpCharInRead != tmpCharInChr)
			{
				tmpAcceptorMapCumulativeMismatchNum ++;
				acceptorMapMismatchPosVec.push_back(toFixSeqLength - tmp);
			}
			acceptorMapCumulativeMismatchNumVec.push_back(tmpAcceptorMapCumulativeMismatchNum);			
		}
	}

	void generateFusionSiteVec_denovo(string& toFixReadSeq, string& toFixChrSeq_left, string& toFixChrSeq_right, 
		Index_Info* indexInfo, bool insertionAllowed_bool, int max_allowed_mismatchNum, 
		int min_allowed_insertionLength, int max_allowed_insertionLength)
	{
		this->generateFusionSiteVec_noIns_denovo(toFixReadSeq, toFixChrSeq_left, toFixChrSeq_right, 
			indexInfo, max_allowed_mismatchNum);
		if(insertionAllowed_bool)
			this->generateFusionSiteVec_insOnly_denovo(toFixReadSeq, toFixChrSeq_left, toFixChrSeq_right, 
				indexInfo, max_allowed_mismatchNum, min_allowed_insertionLength, max_allowed_insertionLength);
	}

	void generateFusionSiteVec_noIns_denovo(string& toFixReadSeq, string& toFixChrSeq_left, string& toFixChrSeq_right, 
		Index_Info* indexInfo, int max_allowed_mismatchNum)
	{
		int toFixReadSeqLength = toFixReadSeq.length();
		for(int tmp_donerLength = 2; tmp_donerLength <= toFixReadSeqLength; tmp_donerLength ++)
		{
			int tmp_acceptorLength = toFixReadSeqLength - tmp_donerLength;
			if(tmp_acceptorLength < 2)
				continue;
			int index_donerMapCumulativeMismatchNumVec = tmp_donerLength - 1;
			int index_acceptorMapCumulativeMismatchNumVec = tmp_acceptorLength - 1;
			int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec];
			int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec];
			int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;
			if(tmp_mismatch_sum <= max_allowed_mismatchNum)
			{
				string tmp_flank_string_left = toFixChrSeq_left.substr(tmp_donerLength, 2);
				string tmp_flank_string_right = toFixChrSeq_right.substr(tmp_donerLength - 2, 2);
				string tmp_flank_string = tmp_flank_string_left + tmp_flank_string_right;
				int tmp_flank_string_case = this->returnFlankStringCase_forward(tmp_flank_string); 
				fusionSiteAnchorSeqLengthPairVec.push_back(pair<int,int>(tmp_donerLength, tmp_acceptorLength));
				fusionSiteInsLengthVec.push_back(0);
				fusionSiteMismatchNumVec.push_back(tmp_mismatch_sum);
				fusionSiteFlankStringCaseVec.push_back(tmp_flank_string_case);
			}
		}
	}
	
	void generateFusionSiteVec_insOnly_denovo(string& toFixReadSeq, string& toFixChrSeq_left, string& toFixChrSeq_right, 
		Index_Info* indexInfo, int max_allowed_mismatchNum, int min_allowed_insertionLength, int max_allowed_insertionLength)
	{
		int toFixReadSeqLength = toFixReadSeq.length();
		for(int tmp_donerLength = 2; tmp_donerLength <= toFixReadSeqLength; tmp_donerLength ++)
		{
			for(int tmp_acceptorLength = 2; tmp_acceptorLength <= toFixReadSeqLength; tmp_acceptorLength ++)
			{
				if(tmp_donerLength + tmp_acceptorLength >= toFixReadSeqLength)
					continue;
				int index_donerMapCumulativeMismatchNumVec = tmp_donerLength - 1;
				int index_acceptorMapCumulativeMismatchNumVec = tmp_acceptorLength - 1;
				int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec];
				int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec];
				int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;		
				if(tmp_mismatch_sum <= max_allowed_mismatchNum)
				{
					int tmp_insertion_length = toFixReadSeqLength - tmp_donerLength - tmp_acceptorLength;
					if((tmp_insertion_length <= max_allowed_insertionLength)&&(tmp_insertion_length >= min_allowed_insertionLength))
					{
						string tmp_flank_string_left = toFixChrSeq_left.substr(tmp_donerLength, 2);
						string tmp_flank_string_right = toFixChrSeq_right.substr(toFixReadSeqLength - tmp_acceptorLength - 2, 2);
						string tmp_flank_string = tmp_flank_string_left + tmp_flank_string_right;
						int tmp_flank_string_case = this->returnFlankStringCase_forward(tmp_flank_string); 
						fusionSiteAnchorSeqLengthPairVec.push_back(pair<int,int>(tmp_donerLength, tmp_acceptorLength));
						fusionSiteInsLengthVec.push_back(tmp_insertion_length);
						fusionSiteMismatchNumVec.push_back(tmp_mismatch_sum);
						fusionSiteFlankStringCaseVec.push_back(tmp_flank_string_case);
					}
				}
			}
		}
	}

	void addNewFusionSite_generateAnnotatedIndexVec_withGeneAnn(string& toFixChrSeq_left, string& toFixChrSeq_right,
		GeneAnnEntry_Hash_Info& geneAnnEntryHashInfo, int toFixSeqLength,
		int chrNameInt_1, int toFixChrSeq_1_startPos, int toFixChrSeq_1_endPos, bool mapStrandSameAsMainBody_left_bool,
		int chrNameInt_2, int toFixChrSeq_2_startPos, int toFixChrSeq_2_endPos, bool mapStrandSameAsMainBody_right_bool,
		int maximum_allowed_mismatchNum, int minimum_allowed_insertionLength, int maximum_allowed_insertionLength)
	{
		//cout << "addNewFusionSite_generateAnnotatedIndexVec_withGeneAnn( starts ..." << endl;
		//cout << "mapStrandSameAsMainBody_left_bool: " << mapStrandSameAsMainBody_left_bool << endl;
		//cout << "mapStrandSameAsMainBody_right_bool: " << mapStrandSameAsMainBody_right_bool << endl;
		vector<int> inBoundaryPosVec_1;
		vector<int> inBoundaryPosVec_2;
		//cout << "chrNameInt_1: " << chrNameInt_1 << " " << toFixChrSeq_1_startPos << " ~ " << toFixChrSeq_1_endPos << endl;
		//cout << "chrNameInt_2: " << chrNameInt_2 << " " << toFixChrSeq_2_startPos << " ~ " << toFixChrSeq_2_endPos << endl;
		geneAnnEntryHashInfo.returnRightBoundaryPosVec(inBoundaryPosVec_1, chrNameInt_1, toFixChrSeq_1_startPos, toFixChrSeq_1_endPos);
		geneAnnEntryHashInfo.returnLeftBoundaryPosVec(inBoundaryPosVec_2, chrNameInt_2, toFixChrSeq_2_startPos, toFixChrSeq_2_endPos);	
	
		vector<int> inBoundaryAnchorLengthVec_1;
		vector<int> inBoundaryAnchorLengthVec_2;
		// cout << "inBoundaryPosVec_1.size(): " << inBoundaryPosVec_1.size() << endl;
		// cout << "inBoundaryPosVec_2.size(): " << inBoundaryPosVec_2.size() << endl;
		for(int tmp = 0; tmp < inBoundaryPosVec_1.size(); tmp++)
		{
			int tmp_inBoundaryPos_1 = inBoundaryPosVec_1[tmp];
			//cout << "tmp_inBoundaryPos_1: " << tmp_inBoundaryPos_1 << endl;
			if(mapStrandSameAsMainBody_left_bool)
				inBoundaryAnchorLengthVec_1.push_back(tmp_inBoundaryPos_1 - toFixChrSeq_1_startPos + 1);
			else
				inBoundaryAnchorLengthVec_1.push_back(toFixChrSeq_1_endPos - tmp_inBoundaryPos_1 + 1);
		}
		//cout << "ends of building inBoundaryAnchorLengthVec_1 ..." << endl;
		for(int tmp = 0; tmp < inBoundaryPosVec_2.size(); tmp++)
		{
			int tmp_inBoundaryPos_2 = inBoundaryPosVec_2[tmp];
			//cout << "tmp_inBoundaryPos_2: " << tmp_inBoundaryPos_2 << endl;
			if(mapStrandSameAsMainBody_right_bool)
				inBoundaryAnchorLengthVec_2.push_back(toFixChrSeq_2_endPos - tmp_inBoundaryPos_2 + 1);
			else
				inBoundaryAnchorLengthVec_2.push_back(tmp_inBoundaryPos_2 - toFixChrSeq_2_startPos + 1);
		}
		// cout << "ends of building inBoundaryAnchorLengthVec_2 ..." << endl;
		// cout << "inBoundaryAnchorLengthVec_1.size(): " << inBoundaryAnchorLengthVec_1.size() << endl;
		// cout << "inBoundaryAnchorLengthVec_2.size(): " << inBoundaryAnchorLengthVec_2.size() << endl;
		int denovoDetectedFusionSiteNum = fusionSiteAnchorSeqLengthPairVec.size();
		//cout << "denovoDetectedFusionSiteNum: " << denovoDetectedFusionSiteNum << endl;
		for(int tmpDonerIndex = 0; tmpDonerIndex < inBoundaryAnchorLengthVec_1.size(); tmpDonerIndex++)
		{
			int tmpDonerLength = inBoundaryAnchorLengthVec_1[tmpDonerIndex];
			//cout << "tmpDonerIndex: " << tmpDonerIndex << "\ttmpDonerLength: " << tmpDonerLength << endl;
			if(tmpDonerLength < 2)
				continue;
			for(int tmpAcceptorIndex = 0; tmpAcceptorIndex < inBoundaryAnchorLengthVec_2.size(); tmpAcceptorIndex ++)
			{
				int tmpAcceptorLength = inBoundaryAnchorLengthVec_2[tmpAcceptorIndex];
				// cout << "\ttmpAcceptorIndex: " << tmpAcceptorIndex 
				// 	<< "\ttmpAcceptorLength: " << tmpAcceptorLength << endl;
				if(tmpAcceptorLength < 2)
					continue;
				if(tmpDonerLength + tmpAcceptorLength != toFixSeqLength)// + ALLOWED_DELETION_IN_READ_DUETO_DREAMCHALLENGE_SIMULATIONERROR)
					continue;
				bool alreadyExistInDenovoDetectedFusionSiteVec_bool = false;
				for(int tmpFusionSiteIndex = 0; tmpFusionSiteIndex < denovoDetectedFusionSiteNum; tmpFusionSiteIndex ++)
				{
					int tmpAnchorLength_1 = fusionSiteAnchorSeqLengthPairVec[tmpFusionSiteIndex].first; // doner anchor length in existing fusionSiteVec
					int tmpAnchorLength_2 = fusionSiteAnchorSeqLengthPairVec[tmpFusionSiteIndex].second; // acceptor anchor length in existing fusionSiteVec
					if((tmpDonerLength == tmpAnchorLength_1)&&(tmpAcceptorLength == tmpAnchorLength_2))// already denovoly detected
					{
						annotatedIndexVecInFusionSiteVec.push_back(tmpFusionSiteIndex);
						alreadyExistInDenovoDetectedFusionSiteVec_bool = true;
						break;
					}
				}
				if(!alreadyExistInDenovoDetectedFusionSiteVec_bool) // can not be denovoly detected
				{
					// add new fusion site to fuisonSiteVec
					int index_donerMapCumulativeMismatchNumVec_new = tmpDonerLength - 1;
					int index_acceptorMapCumulativeMismatchNumVec_new = tmpAcceptorLength - 1;
					int tmpDonerMismatchNum_new = donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec_new];
					int tmpAcceptorMismatchNum_new = acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec_new];
					int tmpTotalMismatchNum_new = tmpDonerMismatchNum_new + tmpAcceptorMismatchNum_new;
					int tmpInsertionLength_new = toFixSeqLength - tmpDonerLength - tmpAcceptorLength;
					if((tmpTotalMismatchNum_new <= maximum_allowed_mismatchNum)&&(tmpInsertionLength_new >= minimum_allowed_insertionLength)
						&&(tmpInsertionLength_new <= maximum_allowed_insertionLength))	
					{
						string tmpDonerFlankString_new = toFixChrSeq_left.substr(tmpDonerLength, 2);
						string tmpAcceptorFlankString_new = toFixChrSeq_right.substr(toFixSeqLength - tmpAcceptorLength - 2, 2);
						string tmpTotalFlankString_new = tmpDonerFlankString_new + tmpAcceptorFlankString_new;
						int tmpTotalFlankStringCase_new = this->returnFlankStringCase_forward(tmpTotalFlankString_new);

						fusionSiteAnchorSeqLengthPairVec.push_back(pair<int,int>(tmpDonerLength, tmpAcceptorLength));
						fusionSiteInsLengthVec.push_back(tmpInsertionLength_new);
						fusionSiteMismatchNumVec.push_back(tmpTotalMismatchNum_new);
						fusionSiteFlankStringCaseVec.push_back(tmpTotalFlankStringCase_new);

						annotatedIndexVecInFusionSiteVec.push_back(fusionSiteAnchorSeqLengthPairVec.size() - 1);
					}
				}
			}
		}
	}

	void generateFusionSiteAnnotatedBoolVec()
	{
		int fusionSiteVecSize = fusionSiteAnchorSeqLengthPairVec.size();
		for(int tmp = 0; tmp < fusionSiteVecSize; tmp++)
			fusionSiteAnnotatedBoolVec.push_back(false);
		int annotatedFusionSiteNum = annotatedIndexVecInFusionSiteVec.size();
		for(int tmp = 0; tmp < annotatedFusionSiteNum; tmp++)
		{
			int tmpAnnotatedFusionSiteIndex = annotatedIndexVecInFusionSiteVec[tmp];
			fusionSiteAnnotatedBoolVec[tmpAnnotatedFusionSiteIndex] = true;
		}
	}

	int selectBestFusionSiteIndex(int maximum_allowed_mismatchNum, int minimum_allowed_insertionLength, 
		int maximum_allowed_insertionLength)
	{
		double currentBestSJ_penalty = 100000.0;
		int currentBestSJ_index = -1;
		int fusionSiteVecSize = fusionSiteAnchorSeqLengthPairVec.size();
		for(int tmp = 0; tmp < fusionSiteVecSize; tmp++)
		{
			int tmpFusionSiteInsLength = fusionSiteInsLengthVec[tmp];
			int tmpFusionSiteMismatchNum = fusionSiteMismatchNumVec[tmp];
			if((tmpFusionSiteMismatchNum > maximum_allowed_mismatchNum)
				||(tmpFusionSiteInsLength < minimum_allowed_insertionLength)
				||(tmpFusionSiteInsLength > maximum_allowed_insertionLength))
				continue;
			int tmpFusionSiteFlankStringCase = fusionSiteFlankStringCaseVec[tmp];
			double tmpPenalty = abs(tmpFusionSiteInsLength) * PER_INSERTION_PENALTY + tmpFusionSiteMismatchNum * PER_MISMATCH_PENALTY 
				+ this->returnFusionSitePenalty(tmpFusionSiteFlankStringCase);
			bool tmpFusionSiteAnnotatedBool = fusionSiteAnnotatedBoolVec[tmp];
			if(tmpFusionSiteAnnotatedBool)
				tmpPenalty = tmpPenalty - ANNOTATED_SCORE;
			if(tmpPenalty < currentBestSJ_penalty)
			{
				currentBestSJ_penalty = tmpPenalty;
				currentBestSJ_index = tmp;
			}
		}
		return currentBestSJ_index;
	}	

	int returnFlankStringCase_forward(const string& flank_string)   
	{
		if(flank_string == "ATAC")
			return 1;
		else if(flank_string == "GCAG")
			return 2;
		else if(flank_string == "GTAG")
			return 3;
		else
			return 0;
	}

	double returnFusionSitePenalty(int flank_string_case)
	{
		return this->returnFusionSiteFlankStringCasePenalty(flank_string_case);
	}

	double returnFusionSiteFlankStringCasePenalty(int flank_string_case)
	{
		if(flank_string_case == 3)
			return CANONICAL_FUSION_PENALTY;
		else if((flank_string_case == 2)||(flank_string_case == 1))
			return SEMICANONICAL_FUSION_PENALTY;
		else if(flank_string_case == 0)
			return NONCANONICAL_FUSION_PENALTY;
		else
		{
			cout << "incorrect flank_string_case in FixDoubleAnchorSplice_Info.h" << endl;
			exit(1);
		}			
	}
};
#endif