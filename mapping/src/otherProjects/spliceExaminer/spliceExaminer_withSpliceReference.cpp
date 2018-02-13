// classify junctions into 3 categories.
// category -- true: annotated
// category -- false: junc with reverse strands, 0-editDistance anchorSeqs
// category -- to classify junctions (others)
// input: 
// 1. junc -- chrName, startPos, endPos, supNum, 
//    (anchorLength_doner, anchorLength_acceptor, flankStringScore (6--can, 2--semi, 0--non))
// 2. reference genome index
// 3. spliceJunctionReference
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

bool parseJuncStr(string& tmpStr, int& tmpJunc_chrNameInt, int& tmpJunc_startPos, int& tmpJunc_endPos, 
	int& tmpJunc_supNum, int& tmpJunc_donerAnchorLength, int& tmpJunc_acceptorAnchorLength, 
	Index_Info* indexInfo, bool anchorLengthProvidedOrNot_bool, int default_anchor_length)
{
	vector<string> tmpJuncFieldVec;
	int tmpStartLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
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
	string tmpJunc_startPosStr = tmpJuncFieldVec[1];
	string tmpJunc_endPosStr = tmpJuncFieldVec[2];
	string tmpJunc_supNumStr = tmpJuncFieldVec[4];
	tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
	if(tmpJunc_chrNameInt < 0)
		return false;
	tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
	tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());	
	tmpJunc_supNum = atoi(tmpJunc_supNumStr.c_str());
	if(anchorLengthProvidedOrNot_bool)
	{
		string tmpJunc_donerAnchorLengthStr = tmpJuncFieldVec[5];
		string tmpJunc_acceptorAnchorLengthStr = tmpJuncFieldVec[6];
		tmpJunc_donerAnchorLength = atoi(tmpJunc_donerAnchorLengthStr.c_str());
		tmpJunc_acceptorAnchorLength = atoi(tmpJunc_acceptorAnchorLengthStr.c_str());
	}
	else
	{
		tmpJunc_donerAnchorLength = default_anchor_length;
		tmpJunc_acceptorAnchorLength = default_anchor_length;
	}
	return true;
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

void alterSpliceAnchorSeqSimi(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, Index_Info* indexInfo,
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

int returnMinPenalty(vector<int>& tmpPenaltyVec)
{
	int tmpMinPenalty = 9999;
	for(int tmp = 0; tmp < tmpPenaltyVec.size(); tmp++)
	{
		int tmpPenalty = tmpPenaltyVec[tmp];
		if(tmpPenalty < tmpMinPenalty)
			tmpMinPenalty = tmpPenalty;
	}
	return tmpMinPenalty;
}

double returnMinPenaltyFloat(vector<double>& tmpPenaltyFloatVec)
{
	double tmpMinPenaltyFloat = 999.0;
	for(int tmp = 0; tmp < tmpPenaltyFloatVec.size(); tmp++)
	{
		double tmpPenaltyFloat = tmpPenaltyFloatVec[tmp];
		if(tmpPenaltyFloat < tmpMinPenaltyFloat)
			tmpMinPenaltyFloat = tmpPenaltyFloat;
	}
	return tmpMinPenaltyFloat;
}

// bool return_canBeExtendedOrAlterSplicedPerfectly_bool(int tmpPenalty_matchThroughAtDoner, int tmpPenalty_matchThroughAtAcceptor,
// 	vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor, vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor)
// {
// 	if((tmpPenalty_matchThroughAtDoner == 0)||(tmpPenalty_matchThroughAtAcceptor == 0))
// 		return true;
// 	else
// 	{
// 		for(int tmp = 0; tmp < tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor.size(); tmp++)
// 		{
// 			if(tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor[tmp] == 0)
// 				return true;
// 		}
// 		for(int tmp = 0; tmp < tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor.size(); tmp++)
// 		{
// 			if(tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor[tmp] == 0)
// 				return true;
// 		}
// 	}
// 	return false;
// }

double return_minPenaltyFloat_extension_alterSplice_seqSimi(
	int tmpJunc_donerAnchorLength, int tmpJunc_acceptorAnchorLength,
	int tmpPenalty_matchThroughAtDoner, int tmpPenalty_matchThroughAtAcceptor,
	vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor, vector<int>& tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor)
{
	double tmpMinPenaltyFloat_matchThroughAtAcceptor = (double)tmpPenalty_matchThroughAtAcceptor / (double)tmpJunc_donerAnchorLength;
	double tmpMinPenaltyFloat_matchThroughAtDoner = (double)tmpPenalty_matchThroughAtDoner / (double)tmpJunc_acceptorAnchorLength;
	int tmpMinPenalty_alterSplice_sharedDoner = returnMinPenalty(tmpAlterSpliceAnchorSeqSimiPenalty_sharedDonor);
	double tmpMinPenaltyFloat_alterSplice_sharedDoner = (double)tmpMinPenalty_alterSplice_sharedDoner / (double)tmpJunc_acceptorAnchorLength;
	int tmpMinPenalty_alterSplice_sharedAcceptor = returnMinPenalty(tmpAlterSpliceAnchorSeqSimiPenalty_sharedAcceptor);
	double tmpMinPenaltyFloat_alterSplice_sharedAcceptor = (double)tmpMinPenalty_alterSplice_sharedAcceptor / (double)tmpJunc_donerAnchorLength;
	vector<double> tmpPenaltyFloatVec;
	tmpPenaltyFloatVec.push_back(tmpMinPenaltyFloat_matchThroughAtAcceptor);
	tmpPenaltyFloatVec.push_back(tmpMinPenaltyFloat_matchThroughAtDoner);
	tmpPenaltyFloatVec.push_back(tmpMinPenaltyFloat_alterSplice_sharedDoner);
	tmpPenaltyFloatVec.push_back(tmpMinPenaltyFloat_alterSplice_sharedAcceptor);
	double tmpMinPenaltyFloat = returnMinPenaltyFloat(tmpPenaltyFloatVec);
	return tmpMinPenaltyFloat;
}

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexPath inputRefSplice inputJunc outputFolder threadNum anchorLengthProvidedOrNotBool" << endl;
		exit(1);
	}
	string threadNumStr = argv[5];
	int threadNum = atoi(threadNumStr.c_str());

	int default_anchor_length = 15;
	double trueJunc_minPenaltyFloat_min = 0.15;
	string anchorLengthProvidedOrNotBoolStr = argv[6];
	bool anchorLengthProvidedOrNot_bool;
	if((anchorLengthProvidedOrNotBoolStr == "true")||(anchorLengthProvidedOrNotBoolStr == "True")
		||(anchorLengthProvidedOrNotBoolStr == "TRUE"))
		anchorLengthProvidedOrNot_bool = true;
	else if((anchorLengthProvidedOrNotBoolStr == "false")||(anchorLengthProvidedOrNotBoolStr == "False")
		||(anchorLengthProvidedOrNotBoolStr == "FALSE"))
		anchorLengthProvidedOrNot_bool = false;
	else
	{
		cout << "anchorLengthProvidedOrNotBoolStr: " << anchorLengthProvidedOrNotBoolStr << endl;
		cout << "should be true or false" << endl;
		exit(1);
	}

	string outputFolderStr = argv[4];
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
	cout << endl << "[" << asctime(local) << "... start to initiate alignInferJunctionHashInfo_ref" << endl;
	string inputRefSpliceFile = argv[2];
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_ref = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_ref->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_ref->insertJuncFromJuncFile_chrNamePos_strand_inclusive_fromGeneAnn(inputRefSpliceFile, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate SJhashInfo_ref from alignInferJunctionHashInfo_ref" << endl;
	SJhash_Info* SJ_ref = new SJhash_Info();
	SJ_ref->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	alignInferJunctionHashInfo_ref->convert2SJhashInfo(SJ_ref, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to classify alignInferJunctionHashInfo_result" << endl;
    string trueJunc_file = outputFolderStr + "true.junc";
	string falseJunc_file = outputFolderStr + "false.junc";
	string toDetermineJunc_file = outputFolderStr + "toDetermine.junc";
	string detailInfoJunc_file = outputFolderStr + "detailInfo.junc";
	ofstream trueJunc_ofs(trueJunc_file.c_str());
	ofstream falseJunc_ofs(falseJunc_file.c_str());
	ofstream toDetermineJunc_ofs(toDetermineJunc_file.c_str());
	ofstream detailInfoJunc_ofs(detailInfoJunc_file.c_str());
	string inputJuncFile = argv[3];
	ifstream inputJunc_ifs(inputJuncFile.c_str());
	int tmpJuncNum = 0;
	while(!inputJunc_ifs.eof())
	{
		string tmpStr;
		getline(inputJunc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpJuncNum ++;
		int tmpThousandIndex = tmpJuncNum / 10000;
		if(tmpJuncNum == tmpThousandIndex * 10000)
			cout << "Processed Junc #: " << tmpJuncNum << endl;
		int tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, tmpJunc_supNum, 
			tmpJunc_donerAnchorLength, tmpJunc_acceptorAnchorLength;
		bool parseJuncStr_bool = parseJuncStr(tmpStr, tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos,
				tmpJunc_supNum, tmpJunc_donerAnchorLength, tmpJunc_acceptorAnchorLength, indexInfo,
				anchorLengthProvidedOrNot_bool, default_anchor_length);
		if(!parseJuncStr_bool)
			continue;
		// check tmpJunc_extensionSeqSimi_doner & tmpJunc_extensionSeqSimi_acceptor
		int tmpPenalty_matchThroughAtDoner, tmpPenalty_matchThroughAtAcceptor;
		extensionSeqSimi(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo,
			tmpJunc_donerAnchorLength, tmpJunc_acceptorAnchorLength,
			tmpPenalty_matchThroughAtDoner, tmpPenalty_matchThroughAtAcceptor);
		// check tmpAlterSiteVec_sharedDonor & tmpAlterSiteVec_sharedAcceptor 
		// 		& tmpAlterSpliceAnchorSeqSimi_sharedDonor & tmpAlterSpliceAnchorSeqSimi_sharedAcceptor
		vector<int>	tmpJunc_alterSiteVec_sharedDonor, tmpJunc_alterSiteVec_sharedAcceptor,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor;
		alterSpliceAnchorSeqSimi(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo,
			tmpJunc_donerAnchorLength, tmpJunc_acceptorAnchorLength, SJ_ref,
			tmpJunc_alterSiteVec_sharedDonor, tmpJunc_alterSiteVec_sharedAcceptor,
			tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor);

		string tmpJuncDetailedInfoStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt) + "\t" 
			+ int_to_str(tmpJunc_startPos) + "\t" + int_to_str(tmpJunc_endPos) + "\tJUNC\t" + int_to_str(tmpJunc_supNum) + "\t"
			+ int_to_str(tmpJunc_donerAnchorLength) + "\t" + int_to_str(tmpJunc_acceptorAnchorLength) + "\t"
			+ int_to_str(tmpPenalty_matchThroughAtDoner) + "\t" + int_to_str(tmpPenalty_matchThroughAtAcceptor) + "\t";

		int tmpJunc_alterSiteVec_sharedDonor_size = tmpJunc_alterSiteVec_sharedDonor.size();
		int tmpJunc_alterSiteVec_sharedAcceptor_size = tmpJunc_alterSiteVec_sharedAcceptor.size();
		// alterSite_sharedDonor
		if(tmpJunc_alterSiteVec_sharedDonor_size == 0)
			tmpJuncDetailedInfoStr += "NONE";
		else
		{	
			for(int tmp = 0; tmp < tmpJunc_alterSiteVec_sharedDonor_size; tmp++)
				tmpJuncDetailedInfoStr = tmpJuncDetailedInfoStr + int_to_str(tmpJunc_alterSiteVec_sharedDonor[tmp]) + ",";
		}
		tmpJuncDetailedInfoStr += "\t";
		// alterSite_sharedAcceptor
		if(tmpJunc_alterSiteVec_sharedAcceptor_size == 0)
			tmpJuncDetailedInfoStr += "NONE";
		else
		{
			for(int tmp = 0; tmp < tmpJunc_alterSiteVec_sharedAcceptor_size; tmp++)
				tmpJuncDetailedInfoStr = tmpJuncDetailedInfoStr + int_to_str(tmpJunc_alterSiteVec_sharedAcceptor[tmp]) + ",";
		}
		tmpJuncDetailedInfoStr += "\t";
		// alterSite_anchorSeqSimi_sharedDonor
		if(tmpJunc_alterSiteVec_sharedDonor_size == 0)
			tmpJuncDetailedInfoStr += "NONE";
		else
		{	
			for(int tmp = 0; tmp < tmpJunc_alterSiteVec_sharedDonor_size; tmp++)
				tmpJuncDetailedInfoStr = tmpJuncDetailedInfoStr + int_to_str(tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor[tmp]) + ",";
		}
		tmpJuncDetailedInfoStr += "\t";
		// alterSite_anchorSeqSimi_sharedAcceptor
		if(tmpJunc_alterSiteVec_sharedAcceptor_size == 0)
			tmpJuncDetailedInfoStr += "NONE";
		else
		{
			for(int tmp = 0; tmp < tmpJunc_alterSiteVec_sharedAcceptor_size; tmp++)
				tmpJuncDetailedInfoStr = tmpJuncDetailedInfoStr + int_to_str(tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor[tmp]) + ",";
		}

		// categorize splice junctions
		bool tmpJunc_annotated_bool = alignInferJunctionHashInfo_ref->SJexistInAlignInferJuncHash(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
		if(tmpJunc_annotated_bool)
			trueJunc_ofs << tmpJuncDetailedInfoStr << endl;
		else
		{
			string tmpJunc_annotatedSJstrand = alignInferJunctionHashInfo_ref->returnRegionOverlappedSJsStrand(
				tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, SJ_ref);
			string tmpJunc_flankString_strand = indexInfo->returnFlankStringStrand(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
			// bool canBeExtendedOrAlterSplicedPerfectly_bool = return_canBeExtendedOrAlterSplicedPerfectly_bool(tmpPenalty_matchThroughAtDoner, 
			// 	tmpPenalty_matchThroughAtAcceptor, tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor);
			// if(canBeExtendedOrAlterSplicedPerfectly_bool)
			// 	falseJunc_ofs << tmpJuncDetailedInfoStr << endl;
			double tmpJunc_minPenaltyFloat = return_minPenaltyFloat_extension_alterSplice_seqSimi(
				tmpJunc_donerAnchorLength, tmpJunc_acceptorAnchorLength, tmpPenalty_matchThroughAtDoner, tmpPenalty_matchThroughAtAcceptor,
				tmpJunc_alterSpliceAnchorSeqSimi_sharedDonor, tmpJunc_alterSpliceAnchorSeqSimi_sharedAcceptor);
			if(tmpJunc_minPenaltyFloat <= trueJunc_minPenaltyFloat_min)
				falseJunc_ofs << tmpJuncDetailedInfoStr << endl;
			else
			{	
				if(tmpJunc_annotatedSJstrand == "N")
					toDetermineJunc_ofs << tmpJuncDetailedInfoStr << endl;
				else if(tmpJunc_annotatedSJstrand == "+")
				{
					if(tmpJunc_flankString_strand == "-")
						falseJunc_ofs << tmpJuncDetailedInfoStr << endl;
					else
						toDetermineJunc_ofs << tmpJuncDetailedInfoStr << endl;
				}	
				else // tmpJunc_annotatedSJstrand == "-"
				{
					if(tmpJunc_flankString_strand == "+")
						falseJunc_ofs << tmpJuncDetailedInfoStr << endl;
					else
						toDetermineJunc_ofs << tmpJuncDetailedInfoStr << endl;
				}
			}
		}
		detailInfoJunc_ofs << tmpJuncDetailedInfoStr << endl;
	}
	detailInfoJunc_ofs.close();
	inputJunc_ifs.close();
	trueJunc_ofs.close();
	falseJunc_ofs.close();
	toDetermineJunc_ofs.close();
	log_ofs.close();
	delete SJ_ref;
	delete alignInferJunctionHashInfo_ref;
	free(chrom);
	delete indexInfo;
	return 0;
}