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
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

int getLeastPenalty(vector<int>& tmpPenaltyVec)
{
	int tmpPenalty_least = 999;
	for(int tmp = 0; tmp < tmpPenaltyVec.size(); tmp ++)
	{
		int tmpPenalty = tmpPenaltyVec[tmp];
		if(tmpPenalty < tmpPenalty_least)
			tmpPenalty_least = tmpPenalty;
	}
	return tmpPenalty_least;
}

double getLeastPenaltyRatio(vector<double>& tmpPenaltyRatioVec)
{
	double tmpPenaltyRatio_least = 1.0;
	int tmpPenaltyRatioVecSize = tmpPenaltyRatioVec.size();
	for(int tmp = 0; tmp < tmpPenaltyRatioVecSize; tmp++)
	{
		double tmpPenaltyRatio = tmpPenaltyRatioVec[tmp];
		if(tmpPenaltyRatio < tmpPenaltyRatio_least)
			tmpPenaltyRatio_least = tmpPenaltyRatio;
	}
	return tmpPenaltyRatio_least;
}

double extensionSeqSimi_returnLeastPenaltyRatio(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo,
	vector<int>& tmpAnchorLengthVec, vector<int>& tmpPenaltyVec)
{
	int defaultAnchorSeqLen_1 = 5;
	int defaultAnchorSeqLen_2 = 10;
		
	int tmpJunc_chrNameInt = tmpChrNameInt;
	int tmpJunc_startPos = tmpStartPos;
	int tmpJunc_endPos = tmpEndPos;

	string tmpJunc_anchorSeq_doner_1 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos - defaultAnchorSeqLen_1 + 1, defaultAnchorSeqLen_1);
	string tmpJunc_anchorSeq_acceptor_1 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos, defaultAnchorSeqLen_1);
	string tmpJunc_extension_atDoner_1 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos + 1, defaultAnchorSeqLen_1);
	string tmpJunc_extension_atAcceptor_1 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos - defaultAnchorSeqLen_1, defaultAnchorSeqLen_1);
	string tmpJunc_anchorSeq_doner_2 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos - defaultAnchorSeqLen_2 + 1, defaultAnchorSeqLen_2);
	string tmpJunc_anchorSeq_acceptor_2 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos, defaultAnchorSeqLen_2);
	string tmpJunc_extension_atDoner_2 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos + 1, defaultAnchorSeqLen_2);
	string tmpJunc_extension_atAcceptor_2 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos - defaultAnchorSeqLen_2, defaultAnchorSeqLen_2);

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtDoner_1;
	tmpFixNWDPinfo_matchThroughAtDoner_1.doNWDP_withMismatchJumpCode(tmpJunc_extension_atDoner_1, tmpJunc_anchorSeq_acceptor_1);
	int penalty_matchThrough_doner_1 = tmpFixNWDPinfo_matchThroughAtDoner_1.getPenalty();
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtAcceptor_1;
	tmpFixNWDPinfo_matchThroughAtAcceptor_1.doNWDP_withMismatchJumpCode(tmpJunc_extension_atAcceptor_1, tmpJunc_anchorSeq_doner_1);
	int penalty_matchThrough_acceptor_1 = tmpFixNWDPinfo_matchThroughAtAcceptor_1.getPenalty();
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtDoner_2;
	tmpFixNWDPinfo_matchThroughAtDoner_2.doNWDP_withMismatchJumpCode(tmpJunc_extension_atDoner_2, tmpJunc_anchorSeq_acceptor_2);
	int penalty_matchThrough_doner_2 = tmpFixNWDPinfo_matchThroughAtDoner_2.getPenalty();
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtAcceptor_2;
	tmpFixNWDPinfo_matchThroughAtAcceptor_2.doNWDP_withMismatchJumpCode(tmpJunc_extension_atAcceptor_2, tmpJunc_anchorSeq_doner_2);
	int penalty_matchThrough_acceptor_2 = tmpFixNWDPinfo_matchThroughAtAcceptor_2.getPenalty();

	double penalty_matchThrough_doner_1_ratio = (double)penalty_matchThrough_doner_1 / (double)defaultAnchorSeqLen_1;
	double penalty_matchThrough_doner_2_ratio = (double)penalty_matchThrough_doner_2 / (double)defaultAnchorSeqLen_2;
	double penalty_matchThrough_acceptor_1_ratio = (double)penalty_matchThrough_acceptor_1 / (double)defaultAnchorSeqLen_1;
	double penalty_matchThrough_acceptor_2_ratio = (double)penalty_matchThrough_acceptor_2 / (double)defaultAnchorSeqLen_2;

	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_1);
	//tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_1);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_2);
	//tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_2);

	int penalty_matchThrough_1 = penalty_matchThrough_doner_1;
	if(penalty_matchThrough_acceptor_1 < penalty_matchThrough_1)
		penalty_matchThrough_1 = penalty_matchThrough_acceptor_1;
	int penalty_matchThrough_2 = penalty_matchThrough_doner_2;
	if(penalty_matchThrough_acceptor_2 < penalty_matchThrough_2)
		penalty_matchThrough_2 = penalty_matchThrough_acceptor_2;	
	tmpPenaltyVec.push_back(penalty_matchThrough_1);
	tmpPenaltyVec.push_back(penalty_matchThrough_2);

	vector<double> tmpPenaltyRatioVec;
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_doner_1_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_doner_2_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_acceptor_1_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_acceptor_2_ratio);

	double tmpPenaltyRatio_least = getLeastPenaltyRatio(tmpPenaltyRatioVec);
	return tmpPenaltyRatio_least;
}

double alterSpliceSeqSimi_returnLeastPenaltyRatio(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo,
	vector<int>& tmpAnchorLengthVec, vector<int>& tmpPenaltyVec, SJhash_Info* SJ)
{
	//cout << "alterSpliceSeqSimi_returnLeastPenaltyRatio( starts ..." << endl;
	int defaultAnchorSeqLen_1 = 5;
	int defaultAnchorSeqLen_2 = 10;

	int tmpSJ_chrNameInt = tmpChrNameInt;
	int tmpSJ_donerEndPos = tmpStartPos;
	int tmpSJ_acceptorStartPos = tmpEndPos;

	// generate anchor sequences
	string tmpSJ_anchorSeq_doner_1 = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_donerEndPos - defaultAnchorSeqLen_1 + 1, defaultAnchorSeqLen_1);
	string tmpSJ_anchorSeq_doner_2 = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_donerEndPos - defaultAnchorSeqLen_2 + 1, defaultAnchorSeqLen_2);
	string tmpSJ_anchorSeq_acceptor_1 = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_acceptorStartPos, defaultAnchorSeqLen_1);
	string tmpSJ_anchorSeq_acceptor_2 = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_acceptorStartPos, defaultAnchorSeqLen_2);
	// generate tmpAlterOtherAcceptorSiteVec & tmpAlterOtherDonerSiteVec
	vector<int> tmpAlterOtherAcceptorSiteVec;
	vector<int> tmpAlterOtherDonerSiteVec;
	vector<int> tmpAlterAcceptorSpliceSitePosVec;
	SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpSJ_chrNameInt,
		tmpSJ_donerEndPos, tmpAlterAcceptorSpliceSitePosVec);
	for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
		if(tmpAlterAcceptorSpliceSitePos != tmpSJ_acceptorStartPos)
			tmpAlterOtherAcceptorSiteVec.push_back(tmpAlterAcceptorSpliceSitePos);
	}
	vector<int> tmpAlterDonerSpliceSitePosVec;
	SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpSJ_chrNameInt,
		tmpSJ_acceptorStartPos, tmpAlterDonerSpliceSitePosVec);	
	for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
		if(tmpAlterDonerSpliceSitePos != tmpSJ_donerEndPos)
			tmpAlterOtherDonerSiteVec.push_back(tmpAlterDonerSpliceSitePos);
	}
	// generate tmpAlterOtherAcceptorSite_anchorSeqVec & tmpAlterOtherDonerSite_anchorSeqVec
	int tmpAlterOtherAcceptorSiteVecSize = tmpAlterOtherAcceptorSiteVec.size();
	vector<string> tmpAlterOtherAcceptorSite_anchorSeqVec_1;
	vector<string> tmpAlterOtherAcceptorSite_anchorSeqVec_2;
	for(int tmp = 0; tmp < tmpAlterOtherAcceptorSiteVecSize; tmp ++)
	{
		int tmpAlterOtherAcceptorSite = tmpAlterOtherAcceptorSiteVec[tmp];
		//cout << "tmpAlterOtherAcceptorSite: " << tmpAlterOtherAcceptorSite << endl;
		string tmpAlterOtherAcceptorSite_anchorSeq_1 
			= indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt, tmpAlterOtherAcceptorSite, defaultAnchorSeqLen_1);
		string tmpAlterOtherAcceptorSite_anchorSeq_2 
			= indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt, tmpAlterOtherAcceptorSite, defaultAnchorSeqLen_2);
		tmpAlterOtherAcceptorSite_anchorSeqVec_1.push_back(tmpAlterOtherAcceptorSite_anchorSeq_1);
		tmpAlterOtherAcceptorSite_anchorSeqVec_2.push_back(tmpAlterOtherAcceptorSite_anchorSeq_2);
	}

	int tmpAlterOtherDonerSiteVecSize = tmpAlterOtherDonerSiteVec.size();
	vector<string> tmpAlterOtherDonerSite_anchorSeqVec_1;
	vector<string> tmpAlterOtherDonerSite_anchorSeqVec_2;
	for(int tmp = 0; tmp < tmpAlterOtherDonerSiteVecSize; tmp ++)
	{
		int tmpAlterOtherDonerSite = tmpAlterOtherDonerSiteVec[tmp];
		//cout << "tmpAlterOtherDonerSite: " << tmpAlterOtherDonerSite << endl;
		string tmpAlterOtherDonerSite_anchorSeq_1 = indexInfo->returnChromStrSubstr(
			tmpSJ_chrNameInt, tmpAlterOtherDonerSite - defaultAnchorSeqLen_1 + 1, defaultAnchorSeqLen_1);
		string tmpAlterOtherDonerSite_anchorSeq_2 = indexInfo->returnChromStrSubstr(
			tmpSJ_chrNameInt, tmpAlterOtherDonerSite - defaultAnchorSeqLen_2 + 1, defaultAnchorSeqLen_2);
		tmpAlterOtherDonerSite_anchorSeqVec_1.push_back(tmpAlterOtherDonerSite_anchorSeq_1);
		tmpAlterOtherDonerSite_anchorSeqVec_2.push_back(tmpAlterOtherDonerSite_anchorSeq_2);
	}

	// generate tmpAlterOtherAcceptorSite_penaltyVec & tmpAlterOtherDonerSite_penaltyVec
	//vector<int> tmpAlterOtherAcceptorSite_penaltyVec_1;
	//vector<int> tmpAlterOtherAcceptorSite_penaltyVec_2;
	vector<int> tmpPenaltyVec_1;
	vector<int> tmpPenaltyVec_2;
	//cout << "start to check penalty in tmpAlterOtherAcceptorSiteVec" << endl;
	for(int tmp = 0; tmp < tmpAlterOtherAcceptorSiteVecSize; tmp ++)
	{
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1, tmpFixNWDPinfo_2;
		tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(tmpSJ_anchorSeq_acceptor_1, tmpAlterOtherAcceptorSite_anchorSeqVec_1[tmp]);
		int tmpPenalty_1 = tmpFixNWDPinfo_1.getPenalty();
		//cout << "penalty_1: " << tmpPenalty_1 << endl;
		//tmpAlterOtherAcceptorSite_penaltyVec_1.push_back(tmpPenalty_1);
		tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(tmpSJ_anchorSeq_acceptor_2, tmpAlterOtherAcceptorSite_anchorSeqVec_2[tmp]);
		int tmpPenalty_2 = tmpFixNWDPinfo_2.getPenalty();
		//cout << "penalty_2: " << tmpPenalty_2 << endl;
		//tmpAlterOtherAcceptorSite_penaltyVec_2.push_back(tmpPenalty_2);
		tmpPenaltyVec_1.push_back(tmpPenalty_1);
		tmpPenaltyVec_2.push_back(tmpPenalty_2);	
	}
	//vector<int> tmpAlterOtherDonerSite_penaltyVec_1;
	//vector<int> tmpAlterOtherDonerSite_penaltyVec_2;
	//cout << "start to check penalty in tmpAlterOtherDonerSiteVec" << endl;
	for(int tmp = 0; tmp < tmpAlterOtherDonerSiteVecSize; tmp ++)
	{
		FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1, tmpFixNWDPinfo_2;
		tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(tmpSJ_anchorSeq_doner_1, tmpAlterOtherDonerSite_anchorSeqVec_1[tmp]);
		int tmpPenalty_1 = tmpFixNWDPinfo_1.getPenalty();
		//cout << "penalty_1: " << tmpPenalty_1 << endl;
		//tmpAlterOtherDonerSite_penaltyVec_1.push_back(tmpPenalty_1);
		tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(tmpSJ_anchorSeq_doner_2, tmpAlterOtherDonerSite_anchorSeqVec_2[tmp]);
		int tmpPenalty_2 = tmpFixNWDPinfo_2.getPenalty();
		//cout << "penalty_2: " << tmpPenalty_2 << endl;
		//tmpAlterOtherDonerSite_penaltyVec_2.push_back(tmpPenalty_2);
		tmpPenaltyVec_1.push_back(tmpPenalty_1);
		tmpPenaltyVec_2.push_back(tmpPenalty_2);	
	}

	int tmpAlterSpliceNum = tmpAlterOtherAcceptorSiteVecSize + tmpAlterOtherDonerSiteVecSize;
	if(tmpAlterSpliceNum == 0)
		return 1.0;
	else
	{
		int tmpPenalty_least_1 = getLeastPenalty(tmpPenaltyVec_1);
		int tmpPenalty_least_2 = getLeastPenalty(tmpPenaltyVec_2);
		tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_1);
		tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_2);
		tmpPenaltyVec.push_back(tmpPenalty_least_1);
		tmpPenaltyVec.push_back(tmpPenalty_least_2);

		vector<double> tmpPenaltyRatioVec;
		double tmpPenaltyRatio_least_1 = (double)tmpPenalty_least_1 / (double)defaultAnchorSeqLen_1;
		double tmpPenaltyRatio_least_2 = (double)tmpPenalty_least_2 / (double)defaultAnchorSeqLen_2;
		tmpPenaltyRatioVec.push_back(tmpPenaltyRatio_least_1);
		tmpPenaltyRatioVec.push_back(tmpPenaltyRatio_least_2);		
		double tmpPenaltyRatio_least = getLeastPenaltyRatio(tmpPenaltyRatioVec);
		return tmpPenaltyRatio_least;	
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputIntropolisJuncFile outputFileWithExtensionSeqSimi" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;
	cout << endl << "[" << asctime(local) << "... start to initiate alignInferJunctionHashInfo" << endl;

	string inputIntropolisJuncFile = argv[2];
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePosOnly(inputIntropolisJuncFile, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate SJhashInfo" << endl;
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to process each junction" << endl;	
	ifstream intropolis_ifs(inputIntropolisJuncFile.c_str());
	string outputFile_extensionSeqSimi_alterSpliceAnchorSeqSimi = argv[3];
	ofstream extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs(outputFile_extensionSeqSimi_alterSpliceAnchorSeqSimi.c_str());
	int tmpJuncNum = 0;
	while(!intropolis_ifs.eof())
	{
		string tmpStr;
		getline(intropolis_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpJuncNum ++;
		int tmpThousandIndex = tmpJuncNum / 100000;
		if(tmpJuncNum == tmpThousandIndex * 100000)
			cout << "Processed Junc #: " << tmpJuncNum << endl;
		//cout << "tmpStr: " << endl << tmpStr << endl;
		vector<string> tmpJuncFieldVec;
		int tmpStartLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
			tmpJuncFieldVec.push_back(tmpFieldStr);
			tmpStartLoc = tmpTabLoc + 1;
		}
		string tmpJunc_chrName = tmpJuncFieldVec[0];
		string tmpJunc_startPosStr = tmpJuncFieldVec[1];
		string tmpJunc_endPosStr = tmpJuncFieldVec[2];
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		
		// extension penalty
		vector<int> tmpAnchorLengthVec_extension; 
		vector<int> tmpPenaltyVec_extension;
		double tmpJunc_extensionSeq_penaltyRatio_least = extensionSeqSimi_returnLeastPenaltyRatio(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, tmpAnchorLengthVec_extension, tmpPenaltyVec_extension);
		extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << tmpStr << "\t" << tmpJunc_extensionSeq_penaltyRatio_least; 
		int tmpVecSize_extension = tmpAnchorLengthVec_extension.size();
		for(int tmp = 0; tmp < tmpVecSize_extension; tmp ++)
			extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << "\t" << tmpAnchorLengthVec_extension[tmp] << ":" << tmpPenaltyVec_extension[tmp];
		
		// alter splice anchor penalty
		vector<int> tmpAnchorLengthVec_alterSplice; 
		vector<int> tmpPenaltyVec_alterSplice;
		double tmpJunc_alterSpliceSeq_penaltyRatio_least = alterSpliceSeqSimi_returnLeastPenaltyRatio(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, tmpAnchorLengthVec_alterSplice, tmpPenaltyVec_alterSplice, SJ);
		int tmpVecSize_alterSplice = tmpAnchorLengthVec_alterSplice.size();
		if(tmpVecSize_alterSplice > 0)
		{	
			extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << "\t" << tmpJunc_alterSpliceSeq_penaltyRatio_least;		
			for(int tmp = 0; tmp < tmpVecSize_alterSplice; tmp ++)
				extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << "\t" << tmpAnchorLengthVec_alterSplice[tmp] << ":" << tmpPenaltyVec_alterSplice[tmp];
		}
		else
			extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << "\t1.0\tNULL\tNULL";

		extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs << endl;
	}	
	extensionSeqSimi_alterSpliceAnchorSeqSimi_ofs.close();
	delete indexInfo;
	delete alignInferJunctionHashInfo;
	free(chrom);
	intropolis_ifs.close();
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}	