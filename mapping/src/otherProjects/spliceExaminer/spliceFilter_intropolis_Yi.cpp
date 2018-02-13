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
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

bool parseIntropolisJuncStr(string& tmpStr, int& tmpJunc_chrNameInt, int& tmpJunc_startPos, 
	int& tmpJunc_endPos, int& tmpJunc_sampleSupNum, int& tmpJunc_readSupNum, Index_Info* indexInfo)
{
	vector<string> tmpJuncFieldVec;
	int tmpStartLoc = 0;
	//cout << "tmpStr " << tmpStr << endl;
	for(int tmp = 0; ; tmp++)
	{
		int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
		//cout << "tmpTabLoc: " << tmpTabLoc << endl;
		if(tmpTabLoc == string::npos)
		{
			tmpJuncFieldVec.push_back(tmpStr.substr(tmpStartLoc));
			break;
		}
		string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
		tmpJuncFieldVec.push_back(tmpFieldStr);
		tmpStartLoc = tmpTabLoc + 1;
	}
	string tmpJunc_chrName = tmpJuncFieldVec[0];
	//cout << "tmpJunc_chrName: " << tmpJunc_chrName << endl;
	tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
	if(tmpJunc_chrNameInt < 0)
		return false;
	string tmpJunc_startPosStr = tmpJuncFieldVec[1];
	tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
	string tmpJunc_endPosStr = tmpJuncFieldVec[2];
	tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
	string tmpJunc_sampleSupNumStr = tmpJuncFieldVec[5];
	tmpJunc_sampleSupNum = atoi(tmpJunc_sampleSupNumStr.c_str());
	string tmpJunc_readSupNumStr = tmpJuncFieldVec[6];
	tmpJunc_readSupNum = atoi(tmpJunc_readSupNumStr.c_str());
	return true;
	exit(1);
}

void extensionSeqSimi(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo, 
	int tmpAchorLength_doner, int tmpAnchorLength_acceptor,
	int& tmpPenalty_matchThroughAtDoner, int& tmpPenalty_matchThroughAtAcceptor)
{
	string tmpAnchorSeq_doner = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
		tmpStartPos - tmpAchorLength_doner + 1, tmpAchorLength_doner);
	string tmpAnchorSeq_acceptor = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
		tmpEndPos, tmpAnchorLength_acceptor);

	string tmpExtensionSeq_atAcceptorSite = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
		tmpEndPos - tmpAchorLength_doner, tmpAchorLength_doner);
	string tmpExtensionSeq_atDonerSite = indexInfo->returnChromStrSubstr(tmpChrNameInt,
		tmpStartPos + 1, tmpAnchorLength_acceptor);
	
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtDoner;
	tmpFixNWDPinfo_matchThroughAtDoner.doNWDP_withMismatchJumpCode(tmpExtensionSeq_atDonerSite, tmpAnchorSeq_acceptor);
	tmpPenalty_matchThroughAtDoner = tmpFixNWDPinfo_matchThroughAtDoner.getPenalty();	

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtAcceptor;
	tmpFixNWDPinfo_matchThroughAtAcceptor.doNWDP_withMismatchJumpCode(tmpExtensionSeq_atAcceptorSite, tmpAnchorSeq_doner);
	tmpPenalty_matchThroughAtAcceptor = tmpFixNWDPinfo_matchThroughAtAcceptor.getPenalty();	
}

void alterSpliceAnchorSeqSimi_cmp2ref(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo,
	int tmpAnchorLength_doner, int tmpAnchorLength_acceptor, SJhash_Info* SJ_ref,
	vector<int>& tmpAlterSiteVec_sharedDonor, vector<int>& tmpAlterSiteVec_sharedAcceptor,
	vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor, vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor)
{
	string tmpAnchorSeq_doner = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
		tmpStartPos - tmpAnchorLength_doner + 1, tmpAnchorLength_doner);
	string tmpAnchorSeq_acceptor = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
		tmpEndPos, tmpAnchorLength_acceptor);
	// generate tmpAlterOtherAcceptorSiteVec & tmpAlterOtherDonerSiteVec
	vector<int> tmpAlterAcceptorSpliceSitePosVec;
	SJ_ref->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChrNameInt, tmpStartPos, tmpAlterAcceptorSpliceSitePosVec);
	vector<int> tmpAlterDonerSpliceSitePosVec;
	SJ_ref->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChrNameInt, tmpEndPos, tmpAlterDonerSpliceSitePosVec);

	for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
		if(tmpAlterAcceptorSpliceSitePos != tmpEndPos)
		{
			tmpAlterSiteVec_sharedDonor.push_back(tmpAlterAcceptorSpliceSitePos);
			string tmpAlterAnchorSeq_sharedDonor = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
				tmpAlterAcceptorSpliceSitePos, tmpAnchorLength_acceptor);
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo;
			tmpFixNWDPinfo.doNWDP_withMismatchJumpCode(tmpAnchorSeq_acceptor, tmpAlterAnchorSeq_sharedDonor);
			int tmpPenalty = tmpFixNWDPinfo.getPenalty();
			tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor.push_back(tmpPenalty);
		}
	}	
	for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
		if(tmpAlterDonerSpliceSitePos != tmpStartPos)
		{
			tmpAlterSiteVec_sharedAcceptor.push_back(tmpAlterDonerSpliceSitePos);
			string tmpAlterAnchorSeq_sharedAcceptor = indexInfo->returnChromStrSubstr(tmpChrNameInt, 
				tmpAlterDonerSpliceSitePos - tmpAnchorLength_doner + 1, tmpAnchorLength_doner);		
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo;
			tmpFixNWDPinfo.doNWDP_withMismatchJumpCode(tmpAnchorSeq_doner, tmpAlterAnchorSeq_sharedAcceptor);
			int tmpPenalty = tmpFixNWDPinfo.getPenalty();
			tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor.push_back(tmpPenalty);
		}
	}
}

void alterSpliceAnchorSeqSimi_cmp2self(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo,
	int tmpAnchorLength_doner, int tmpAnchorLength_acceptor, SJhash_Info* SJ_self, AlignInferJunctionHash_Info* alignInferJuncHash_self,
	vector<int>& tmpAlterSiteVec_sharedDonor, vector<int>& tmpAlterSiteVec_sharedAcceptor,
	vector<int>& tmpJunc_alterSiteSupReadNum_sharedDonor, vector<int>& tmpJunc_alterSiteSupReadNum_sharedAcceptor,
	vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor, vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor)
{
	string tmpAnchorSeq_doner = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpStartPos - tmpAnchorLength_doner + 1, tmpAnchorLength_doner);
	string tmpAnchorSeq_acceptor = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpEndPos, tmpAnchorLength_acceptor);
	// generate tmpAlterOtherAcceptorSiteVec & tmpAlterOtherDonerSiteVec
	vector<int> tmpAlterAcceptorSpliceSitePosVec;
	SJ_self->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChrNameInt, tmpStartPos, tmpAlterAcceptorSpliceSitePosVec);
	vector<int> tmpAlterDonerSpliceSitePosVec;
	SJ_self->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChrNameInt, tmpEndPos, tmpAlterDonerSpliceSitePosVec);

	for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
		if(tmpAlterAcceptorSpliceSitePos != tmpEndPos)
		{
			tmpAlterSiteVec_sharedDonor.push_back(tmpAlterAcceptorSpliceSitePos);
			string tmpAlterAnchorSeq_sharedDonor = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpAlterAcceptorSpliceSitePos, tmpAnchorLength_acceptor);
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo;
			tmpFixNWDPinfo.doNWDP_withMismatchJumpCode(tmpAnchorSeq_acceptor, tmpAlterAnchorSeq_sharedDonor);
			int tmpPenalty = tmpFixNWDPinfo.getPenalty();
			tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor.push_back(tmpPenalty);
			int tmpSupReadNum; 
			bool tmpSearchSuccessBool = alignInferJuncHash_self->searchAndReturnSupNumInAlignInferJuncHash(tmpSupReadNum, tmpChrNameInt, tmpStartPos, tmpAlterAcceptorSpliceSitePos);
			if(!tmpSearchSuccessBool) // FIX ME, should not be false
				tmpSupReadNum = 0;
			tmpJunc_alterSiteSupReadNum_sharedDonor.push_back(tmpSupReadNum);
		}
	}	
	for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
		if(tmpAlterDonerSpliceSitePos != tmpStartPos)
		{
			tmpAlterSiteVec_sharedAcceptor.push_back(tmpAlterDonerSpliceSitePos);
			string tmpAlterAnchorSeq_sharedAcceptor = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpAlterDonerSpliceSitePos - tmpAnchorLength_doner + 1, tmpAnchorLength_doner);		
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo;
			tmpFixNWDPinfo.doNWDP_withMismatchJumpCode(tmpAnchorSeq_doner, tmpAlterAnchorSeq_sharedAcceptor);
			int tmpPenalty = tmpFixNWDPinfo.getPenalty();
			tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor.push_back(tmpPenalty);
			int tmpSupReadNum;
			bool tmpSearchSuccessBool = alignInferJuncHash_self->searchAndReturnSupNumInAlignInferJuncHash(tmpSupReadNum, tmpChrNameInt, tmpAlterDonerSpliceSitePos, tmpEndPos);
			if(!tmpSearchSuccessBool)
				tmpSupReadNum = 0;
			tmpJunc_alterSiteSupReadNum_sharedAcceptor.push_back(tmpSupReadNum);
		}
	}
}

int getTheMinimumPenalty_relativeLowSupport(int tmpSupNum, vector<int>& tmpAlterSiteSupNumVec, vector<int>& tmpAlterSitePenaltyVec)
{
	int minimumPenalty = 99;
	for(int tmp = 0; tmp < tmpAlterSiteSupNumVec.size(); tmp++)
	{
		int tmpAlterSiteSupNum = tmpAlterSiteSupNumVec[tmp];
		int tmpAlterSitePenalty = tmpAlterSitePenaltyVec[tmp];
		if((tmpSupNum <= tmpAlterSiteSupNum)&&(tmpAlterSitePenalty < minimumPenalty))
			minimumPenalty = tmpAlterSitePenalty;
	}
	return minimumPenalty;
}

int getTheMinimumPenalty(vector<int>& tmpPenaltyVec)
{
	int minimumPenalty = 99;
	for(int tmp = 0; tmp < tmpPenaltyVec.size(); tmp ++)
	{
		int tmpPenalty = tmpPenaltyVec[tmp];
		if(tmpPenalty < minimumPenalty)
			minimumPenalty = tmpPenalty;
	}
	return minimumPenalty;
}

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable inputIndexPath inputAnnSplice inputToFilterJunc defaultAnchorSeqLength " << endl;
		cout << " defaultSimilarAnchorSeqEditDistanceMax minSupNum_read minSupNum_sample outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[8];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

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
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;

	cout << endl << "[" << asctime(local) << "... start to initiate alignInferJunctionHashInfo_ref and SJhash_Info_ref" << endl;
	string inputRefSpliceFile = argv[2];
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_ref = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_ref->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_ref->insertJuncFromJuncFile_chrNamePos_strand_inclusive_fromGeneAnn(inputRefSpliceFile, indexInfo);
	SJhash_Info* SJ_ref = new SJhash_Info();
	SJ_ref->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	alignInferJunctionHashInfo_ref->convert2SJhashInfo(SJ_ref, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of initiating alignInferJunctionHashInfo_ref and SJhash_Info_ref" << endl;

	cout << endl << "[" << asctime(local) << "... start to initiate alignInferJunctionHashInfo_intropolis and SJ_intropolis" << endl;
	string inputToFilterJuncFile = argv[3];
	cout << "inputToFilterJunc: " << inputToFilterJuncFile << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_intropolis = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_intropolis->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_intropolis->insertJuncFromJuncFile_chrNamePos_readSupNum_intropolis(inputToFilterJuncFile, indexInfo);
	SJhash_Info* SJ_intropolis = new SJhash_Info();
	SJ_intropolis->initiateAreaAndStringHash(chromNum);
	alignInferJunctionHashInfo_intropolis->convert2SJhashInfo(SJ_intropolis, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of initiating alignInferJunctionHashInfo_intropolis and SJ_intropolis" << endl;

	cout << endl << "[" << asctime(local) << "... start to filter possible false positives from intropolis" << endl;
	string defaultAnchorSeqLengthStr = argv[4];
	int defaultAnchorSeqLength = atoi(defaultAnchorSeqLengthStr.c_str());
	string defaultSimilarAnchorSeqEditDistanceMaxStr = argv[5];
	int defaultSimilarAnchorSeqEditDistanceMax = atoi(defaultSimilarAnchorSeqEditDistanceMaxStr.c_str());
	double default_similarAnchorSeq_editDistanceRateMax = (double)defaultSimilarAnchorSeqEditDistanceMax/(double)defaultAnchorSeqLength;
	string minSupNum_read_str = argv[6];
	int	minSupNum_read = atoi(minSupNum_read_str.c_str());
	string minSupNum_sample_str = argv[7];
	int minSupNum_sample = atoi(minSupNum_sample_str.c_str());

	string outputJuncFile_kept = outputFolderStr + "kept.junc";
	string outputJuncFile_filteredOut = outputFolderStr + "filteredOut.junc";
	string outputJuncFile_invalid = outputFolderStr + "invalid.junc";
	string outputJuncFile_ori = outputFolderStr + "ori.junc";
	ofstream kept_ofs(outputJuncFile_kept.c_str());
	ofstream filteredOut_ofs(outputJuncFile_filteredOut.c_str());
	ofstream invalid_ofs(outputJuncFile_invalid.c_str());
	ofstream ori_ofs(outputJuncFile_ori.c_str());

	ifstream toFilterJunc_ifs(inputToFilterJuncFile.c_str());
	while(!toFilterJunc_ifs.eof())
	{
		string tmpStr;
		getline(toFilterJunc_ifs, tmpStr);
		cout << "tmpStr: " << tmpStr << endl;
		if(tmpStr == "")
			break;
		int tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, tmpJunc_sampleSupNum, tmpJunc_readSupNum;
		bool parse_success_bool = parseIntropolisJuncStr(tmpStr, tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, 
			tmpJunc_sampleSupNum, tmpJunc_readSupNum, indexInfo);
		if(!parse_success_bool)
		{
			invalid_ofs << tmpStr << endl;
			continue;
		}
		ori_ofs << tmpStr << endl;
		if((tmpJunc_sampleSupNum < minSupNum_sample)||(tmpJunc_readSupNum < minSupNum_read))
		{
			filteredOut_ofs << tmpStr << endl;
			continue;
		}
		// check tmpJunc_extensionSeqSimi_doner & tmpJunc_extensionSeqSimi_acceptor
		int tmpPenalty_matchThroughAtDoner, tmpPenalty_matchThroughAtAcceptor;
		extensionSeqSimi(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, defaultAnchorSeqLength, 
			defaultAnchorSeqLength, tmpPenalty_matchThroughAtDoner, tmpPenalty_matchThroughAtAcceptor);
		double tmpPenalty_matchThroughAtDoner_perc = (double)tmpPenalty_matchThroughAtDoner / (double)defaultAnchorSeqLength;
		double tmpPenalty_matchThroughAtAcceptor_perc = (double)tmpPenalty_matchThroughAtAcceptor / (double)defaultAnchorSeqLength;
		if((tmpPenalty_matchThroughAtDoner_perc <= default_similarAnchorSeq_editDistanceRateMax)
			||(tmpPenalty_matchThroughAtAcceptor_perc <= default_similarAnchorSeq_editDistanceRateMax))
		{
			filteredOut_ofs << tmpStr << endl;
			continue;
		}

		vector<int>	tmpJunc_alterSiteVec_sharedDonor_cmp2ref, tmpJunc_alterSiteVec_sharedAcceptor_cmp2ref,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2ref, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2ref;
		alterSpliceAnchorSeqSimi_cmp2ref(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, defaultAnchorSeqLength, defaultAnchorSeqLength, SJ_ref,
			tmpJunc_alterSiteVec_sharedDonor_cmp2ref, tmpJunc_alterSiteVec_sharedAcceptor_cmp2ref,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2ref, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2ref);
		int tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref = getTheMinimumPenalty(tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2ref);
		int tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref = getTheMinimumPenalty(tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2ref);
		double tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref_perc = (double)tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref / (double)defaultAnchorSeqLength;
		double tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref_perc = (double)tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref / (double)defaultAnchorSeqLength;
		if((tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref_perc <= default_similarAnchorSeq_editDistanceRateMax)
			||(tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref_perc <= default_similarAnchorSeq_editDistanceRateMax))
		{
			filteredOut_ofs << tmpStr << endl;
			continue;
		}	

		vector<int>	tmpJunc_alterSiteVec_sharedDonor_cmp2self, tmpJunc_alterSiteVec_sharedAcceptor_cmp2self,
			tmpJunc_alterSiteSupReadNum_sharedDonor_cmp2self, tmpJunc_alterSiteSupReadNum_sharedAcceptor_cmp2self,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2self, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2self;		
		alterSpliceAnchorSeqSimi_cmp2self(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, defaultAnchorSeqLength, defaultAnchorSeqLength, SJ_intropolis, alignInferJunctionHashInfo_intropolis,
			tmpJunc_alterSiteVec_sharedDonor_cmp2self, tmpJunc_alterSiteVec_sharedAcceptor_cmp2self,
			tmpJunc_alterSiteSupReadNum_sharedDonor_cmp2self, tmpJunc_alterSiteSupReadNum_sharedAcceptor_cmp2self,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2self, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2self);
		int tmpPenalty_alterSpliceSite_sharedDonor_cmp2self = getTheMinimumPenalty_relativeLowSupport(
			tmpJunc_readSupNum, tmpJunc_alterSiteVec_sharedDonor_cmp2self, tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor_cmp2self);
		int tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2self = getTheMinimumPenalty_relativeLowSupport(
			tmpJunc_readSupNum, tmpJunc_alterSiteVec_sharedAcceptor_cmp2self, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor_cmp2self);
		double tmpPenalty_alterSpliceSite_sharedDonor_cmp2self_perc = (double)tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref / (double)defaultAnchorSeqLength;
		double tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2self_perc = (double)tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref / (double)defaultAnchorSeqLength;
		if((tmpPenalty_alterSpliceSite_sharedDonor_cmp2ref_perc <= default_similarAnchorSeq_editDistanceRateMax)
			||(tmpPenalty_alterSpliceSite_sharedAcceptor_cmp2ref_perc <= default_similarAnchorSeq_editDistanceRateMax))
		{
			filteredOut_ofs << tmpStr << endl;
			continue;
		}

		kept_ofs << tmpStr << endl;
	}
	toFilterJunc_ifs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of filtering possible false positives from intropolis" << endl;	


	kept_ofs.close();
	filteredOut_ofs.close();
	invalid_ofs.close();
	ori_ofs.close();

	delete alignInferJunctionHashInfo_intropolis;
	delete SJ_intropolis;
	delete alignInferJunctionHashInfo_ref;
	delete SJ_ref;
	free(chrom);
	delete indexInfo;
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	log_ofs.close();
	return 0;
}