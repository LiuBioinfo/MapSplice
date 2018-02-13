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
	int defaultAnchorSeqLen_3 = 20;
		
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

	string tmpJunc_anchorSeq_doner_3 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos - defaultAnchorSeqLen_3 + 1, defaultAnchorSeqLen_3);
	string tmpJunc_anchorSeq_acceptor_3 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos, defaultAnchorSeqLen_3);
	string tmpJunc_extension_atDoner_3 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_startPos + 1, defaultAnchorSeqLen_3);
	string tmpJunc_extension_atAcceptor_3 = indexInfo->returnChromStrSubstr(tmpJunc_chrNameInt, tmpJunc_endPos - defaultAnchorSeqLen_3, defaultAnchorSeqLen_3);				
		
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

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtDoner_3;
	tmpFixNWDPinfo_matchThroughAtDoner_3.doNWDP_withMismatchJumpCode(tmpJunc_extension_atDoner_3, tmpJunc_anchorSeq_acceptor_3);
	int penalty_matchThrough_doner_3 = tmpFixNWDPinfo_matchThroughAtDoner_3.getPenalty();

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtAcceptor_3;
	tmpFixNWDPinfo_matchThroughAtAcceptor_3.doNWDP_withMismatchJumpCode(tmpJunc_extension_atAcceptor_3, tmpJunc_anchorSeq_doner_3);
	int penalty_matchThrough_acceptor_3 = tmpFixNWDPinfo_matchThroughAtAcceptor_3.getPenalty();

	double penalty_matchThrough_doner_1_ratio = (double)penalty_matchThrough_doner_1 / (double)defaultAnchorSeqLen_1;
	double penalty_matchThrough_doner_2_ratio = (double)penalty_matchThrough_doner_2 / (double)defaultAnchorSeqLen_2;
	double penalty_matchThrough_doner_3_ratio = (double)penalty_matchThrough_doner_3 / (double)defaultAnchorSeqLen_3;
	double penalty_matchThrough_acceptor_1_ratio = (double)penalty_matchThrough_acceptor_1 / (double)defaultAnchorSeqLen_1;
	double penalty_matchThrough_acceptor_2_ratio = (double)penalty_matchThrough_acceptor_2 / (double)defaultAnchorSeqLen_2;
	double penalty_matchThrough_acceptor_3_ratio = (double)penalty_matchThrough_acceptor_3 / (double)defaultAnchorSeqLen_3;

	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_1);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_1);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_2);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_2);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_3);
	tmpAnchorLengthVec.push_back(defaultAnchorSeqLen_3);
	tmpPenaltyVec.push_back(penalty_matchThrough_doner_1);
	tmpPenaltyVec.push_back(penalty_matchThrough_acceptor_1);
	tmpPenaltyVec.push_back(penalty_matchThrough_doner_2);
	tmpPenaltyVec.push_back(penalty_matchThrough_acceptor_2);
	tmpPenaltyVec.push_back(penalty_matchThrough_doner_3);
	tmpPenaltyVec.push_back(penalty_matchThrough_acceptor_3);

	vector<double> tmpPenaltyRatioVec;
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_doner_1_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_doner_2_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_doner_3_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_acceptor_1_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_acceptor_2_ratio);
	tmpPenaltyRatioVec.push_back(penalty_matchThrough_acceptor_3_ratio);

	double tmpPenaltyRatio_least = getLeastPenaltyRatio(tmpPenaltyRatioVec);
	return tmpPenaltyRatio_least;
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
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;

	string inputIntropolisJuncFile = argv[2];
	ifstream intropolis_ifs(inputIntropolisJuncFile.c_str());
	string outputFileWithExtensionSeqSimi = argv[3];
	ofstream extensionSeqSimi_ofs(outputFileWithExtensionSeqSimi.c_str());
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
		vector<int> tmpAnchorLengthVec; 
		vector<int> tmpPenaltyVec;
		double tmpJunc_extensionSeq_penaltyRatio_least = extensionSeqSimi_returnLeastPenaltyRatio(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos, indexInfo, tmpAnchorLengthVec, tmpPenaltyVec);
		extensionSeqSimi_ofs << tmpStr << "\t" << tmpJunc_extensionSeq_penaltyRatio_least; 
		int tmpVecSize = tmpAnchorLengthVec.size();
		for(int tmp = 0; tmp < tmpVecSize; tmp ++)
			extensionSeqSimi_ofs << "\t" << tmpAnchorLengthVec[tmp] << ":" << tmpPenaltyVec[tmp];
		extensionSeqSimi_ofs << endl;		
	}	
	extensionSeqSimi_ofs.close();
	delete indexInfo;
	free(chrom);
	intropolis_ifs.close();
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}	