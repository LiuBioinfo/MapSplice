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
using namespace std;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputFusionBedPe outputFolder" << endl;
		exit(1);
	}
	int fusionSiteAnchorSeqLength = 12;

	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[3];
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

	string querySequenceFile = outputFolderStr + "querySeq.txt";
	string queryIdFile = outputFolderStr + "queryId.txt";
	ofstream querySeq_ofs(querySequenceFile.c_str());
	ofstream queryId_ofs(queryIdFile.c_str());
	string inputFusionBedPeFile = argv[2];
	ifstream fusionBedPe_ifs(inputFusionBedPeFile.c_str());
	while(!fusionBedPe_ifs.eof())
	{
		string tmpStr;
		getline(fusionBedPe_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpBedPeFieldVec;		
		parseStr2fieldVec(tmpBedPeFieldVec, tmpStr);
		string tmpChrName_1 = tmpBedPeFieldVec[0];
		int tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrName_1);
		string tmpChrPosStr_1 = tmpBedPeFieldVec[1];
		int tmpChrPos_1 = atoi(tmpChrPosStr_1.c_str());
		string tmpChrName_2 = tmpBedPeFieldVec[2];
		int tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrName_2);
		string tmpChrPosStr_2 = tmpBedPeFieldVec[3];
		int tmpChrPos_2 = atoi(tmpChrPosStr_2.c_str());
		string tmpStrand_1 = tmpBedPeFieldVec[6];
		string tmpStrand_2 = tmpBedPeFieldVec[7];
		
		string tmpFusionSiteAnchorSeq_1, tmpFusionSiteAnchorSeq_2;
		if(tmpStrand_1 == "+") 
			tmpFusionSiteAnchorSeq_1 = indexInfo->returnChromStrSubstr(tmpChrNameInt_1, 
				tmpChrPos_1 - fusionSiteAnchorSeqLength + 1, fusionSiteAnchorSeqLength);
		else
			tmpFusionSiteAnchorSeq_1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				tmpChrNameInt_1, tmpChrPos_1, fusionSiteAnchorSeqLength));
		if(tmpStrand_2 == "+")
			tmpFusionSiteAnchorSeq_2 = indexInfo->returnChromStrSubstr(tmpChrNameInt_2, 
				tmpChrPos_2, fusionSiteAnchorSeqLength);
		else
			tmpFusionSiteAnchorSeq_2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
				tmpChrNameInt_2, tmpChrPos_2 - fusionSiteAnchorSeqLength + 1, fusionSiteAnchorSeqLength));
		string regularFusionSeq = tmpFusionSiteAnchorSeq_1 + tmpFusionSiteAnchorSeq_2;
		string oneBaseOutAtDonerFusionSeq = tmpFusionSiteAnchorSeq_1.substr(0, tmpFusionSiteAnchorSeq_1.length() - 1) + tmpFusionSiteAnchorSeq_2;
		string oneBaseOutAtAcceptorFusionSeq = tmpFusionSiteAnchorSeq_1 + tmpFusionSiteAnchorSeq_2.substr(1, tmpFusionSiteAnchorSeq_2.length() - 1);
		querySeq_ofs << regularFusionSeq << endl << oneBaseOutAtDonerFusionSeq << endl << oneBaseOutAtAcceptorFusionSeq << endl;
		queryId_ofs << tmpStr << "\t" << "REGULAR" << endl << tmpStr << "\tDROP_ONE_BASE_AT_DONOR" << endl << tmpStr << "\tDROP_ONE_BASE_AT_ACCEPTOR" << endl;
	}
	fusionBedPe_ifs.close();

	querySeq_ofs.close();
	queryId_ofs.close();
	log_ofs.close();
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	return 0;
}