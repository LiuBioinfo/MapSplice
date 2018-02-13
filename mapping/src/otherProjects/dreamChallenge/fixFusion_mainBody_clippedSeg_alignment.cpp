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
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
#include "../../general/read_info.h"
#include "general/geneAnnEntryHash.h"
#include "general/fixFusion_peRead.h"

using namespace std;

bool end1orEnd2(int tmpFlag)
{
	return (tmpFlag & 0x40);
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

bool parse_sam_clippedSeg_mainBody(string& tmpClippedSegSam, string& tmpMainBodySam_1, string& tmpMainBodySam_2, Index_Info* indexInfo,
	int& tmp_chrNameInt_mainBody, int& tmp_chrPos_mainBody, int& tmp_chrNameInt_clippedSeg, int& tmp_chrPos_clippedSeg,
	vector<Jump_Code>& tmp_jumpCodeVec_mainBody, vector<Jump_Code>& tmp_jumpCodeVec_clippedSeg,
	bool& tmp_mainBody_forOrRevMap_bool, bool& tmp_clippedSegAt_downstreamOrUpstream_end_bool, bool& tmp_clippedSeg_forOrRevMap_bool, 
	PE_Read_Info& tmpPeReadInfo)
{
	cout << "parse_sam_clippedSeg_mainBody starts ..." << endl;
	vector<string> tmpClippedSegSam_fieldVec;
	parseStr2fieldVec(tmpClippedSegSam_fieldVec, tmpClippedSegSam);
	int chrNameInt_clippedSeg = indexInfo->convertStringToInt(tmpClippedSegSam_fieldVec[1]); //cout << "chrNameInt_clippedSeg: " << chrNameInt_clippedSeg << endl;
	int chrPos_clippedSeg = atoi(tmpClippedSegSam_fieldVec[2].c_str()); //cout << "chrPos_clippedSeg: " << chrPos_clippedSeg << endl;
	string cigarString_clippedSeg = tmpClippedSegSam_fieldVec[3]; //cout << "cigarString_clippedSeg: " << cigarString_clippedSeg << endl;
	int clippedSeg_mapped_read_end_NO = atoi(tmpClippedSegSam_fieldVec[4].c_str());
	int clippedSeg_mapStrandCompared2mainBody = atoi(tmpClippedSegSam_fieldVec[5].c_str()); // 1 -- the same strand cmpTo mainBody; 0 -- reverse strand cmpTo mainBody
	cout << "clippedSeg_mapped_read_end_NO: " << clippedSeg_mapped_read_end_NO << endl;
	cout << "clippedSeg_mapStrandCompared2mainBody: " << clippedSeg_mapStrandCompared2mainBody << endl;

	vector<string> tmpMainBodySam_1_fieldVec;
	cout << "tmpMainBodySam_1: " << endl << tmpMainBodySam_1 << endl;
	parseStr2fieldVec(tmpMainBodySam_1_fieldVec, tmpMainBodySam_1);
	cout << "tmpMainBodySam_1_fieldVec.size(): " << tmpMainBodySam_1_fieldVec.size() << endl;
	int flag_1 = atoi(tmpMainBodySam_1_fieldVec[1].c_str());
	int chrNameInt_1 = indexInfo->convertStringToInt(tmpMainBodySam_1_fieldVec[2]);
	cout << "flag_1: " << flag_1 << endl;
	cout << "chrNameInt_1: " << chrNameInt_1 << endl;
	int chrPos_1 = atoi(tmpMainBodySam_1_fieldVec[3].c_str());
	cout << "chrPos_1: " << chrPos_1 << endl;
	string cigarString_1 = tmpMainBodySam_1_fieldVec[5];
	string readSeq_1 = tmpMainBodySam_1_fieldVec[9];
	cout << "cigarString_1: " << cigarString_1 << endl;
	cout << "readSeq_1: " << readSeq_1 << endl;

	vector<string> tmpMainBodySam_2_fieldVec;
	parseStr2fieldVec(tmpMainBodySam_2_fieldVec, tmpMainBodySam_2);
	int flag_2 = atoi(tmpMainBodySam_2_fieldVec[1].c_str());
	int chrNameInt_2 = indexInfo->convertStringToInt(tmpMainBodySam_2_fieldVec[2]);
	cout << "flag_2: " << flag_1 << endl;
	cout << "chrNameInt_2: " << chrNameInt_2 << endl;
	int chrPos_2 = atoi(tmpMainBodySam_2_fieldVec[3].c_str());
	cout << "chrPos_2: " << chrPos_2 << endl;	
	string cigarString_2 = tmpMainBodySam_2_fieldVec[5];
	string readSeq_2 = tmpMainBodySam_2_fieldVec[9];
	cout << "cigarString_2: " << cigarString_2 << endl;
	cout << "readSeq_2: " << readSeq_2 << endl;

	tmp_chrNameInt_clippedSeg = chrNameInt_clippedSeg; 
	tmp_chrPos_clippedSeg = chrPos_clippedSeg;
	cout << "tmp_chrNameInt_clippedSeg: " << tmp_chrNameInt_clippedSeg << endl;
	cout << "tmp_chrPos_clippedSeg: " << tmp_chrPos_clippedSeg << endl;
	cigarString2jumpCodeVec(cigarString_clippedSeg, tmp_jumpCodeVec_clippedSeg);
	tmp_clippedSegAt_downstreamOrUpstream_end_bool = (clippedSeg_mapped_read_end_NO == 2); // 1 -- upstream end; 2 -- downstream end
	tmp_clippedSeg_forOrRevMap_bool = (clippedSeg_mapStrandCompared2mainBody == 1);
	cout << "tmp_clippedSegAt_downstreamOrUpstream_end_bool: " << tmp_clippedSegAt_downstreamOrUpstream_end_bool << endl;
	cout << "tmp_clippedSeg_forOrRevMap_bool: " << tmp_clippedSeg_forOrRevMap_bool << endl;
	if(tmp_clippedSegAt_downstreamOrUpstream_end_bool)
	{
		tmp_chrNameInt_mainBody = chrNameInt_2;
		tmp_chrPos_mainBody = chrPos_2;
		cigarString2jumpCodeVec(cigarString_2, tmp_jumpCodeVec_mainBody);
	}
	else
	{
		tmp_chrNameInt_mainBody = chrNameInt_1;
		tmp_chrPos_mainBody = chrPos_1;
		cigarString2jumpCodeVec(cigarString_1, tmp_jumpCodeVec_mainBody);
	}
	bool seq1_1st_end_bool = end1orEnd2(flag_1);
	//bool seq2_1st_end_bool = end1orEnd2(flag_2);	
	string rawSeq_1, rawSeq_2;
	cout << "readSeq_1.length(): " << readSeq_1.length() << endl << readSeq_1 << endl;
	//cout << "the LastBase: " << readSeq_1.at(100) << endl;
	cout << "readSeq_2.length(): " << readSeq_2.length() << endl << readSeq_2 << endl;
	//cout << "the LastBase: " << readSeq_2.at(100) << endl;
	if(seq1_1st_end_bool)
	{
		tmp_mainBody_forOrRevMap_bool = true;
		rawSeq_1 = readSeq_1;
		rawSeq_2 = convertStringToReverseComplement(readSeq_2);
	}
	else
	{
		tmp_mainBody_forOrRevMap_bool = false;
		rawSeq_1 = convertStringToReverseComplement(readSeq_2);
		rawSeq_2 = readSeq_1;
	}
	cout << "rawSeq_1: " << endl << rawSeq_1 << endl;
	cout << "rawSeq_2: " << endl << rawSeq_2 << endl;
	// if(seq1_for_or_rcm_bool)
	// 	rawSeq_1 = readSeq_1;
	// else
	// 	rawSeq_1 = convertStringToReverseComplement(readSeq_1);
	// if(seq2_for_or_rcm_bool)
	// 	rawSeq_2 = readSeq_2;
	// else
	// 	rawSeq_2 = convertStringToReverseComplement(readSeq_2);	
	tmpPeReadInfo.initiateReadInfo(tmpMainBodySam_1_fieldVec[0], tmpMainBodySam_2_fieldVec[0], rawSeq_1, rawSeq_2, "*", "*", true, false);
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputGeneAnnEntryFile inputSam_mainBody_clipppedSeg outputFolder" << endl;
		exit(1);
	}
	int maximum_allowed_mismatchNum = 5; 
	int maximum_allowed_insertionLength = 10;
	bool geneAnnEntryBoundaryOnly_bool = false;
	bool geneAnnIncorporated_bool = true; 
	bool insertionAllowed_bool = true;

	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[4];
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

	string inputGeneAnnEntryFile = argv[2];
	cout << "start to load simplified gene annotation file" << endl;
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	cout << "start to initiate GeneAnnEntryArea2infoIndexMapVec" << endl;
	//cout << "FAILED to load geneAnn !!!!!!!!" << endl;
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);

	cout << "start to fixDoubleAnchorFusion" << endl;
	string inputSam_mainBody_clipppedSeg = argv[3];
	ifstream sam_mainBody_clippedSeg_ifs(inputSam_mainBody_clipppedSeg.c_str());
	while(!sam_mainBody_clippedSeg_ifs.eof())
	{
		string tmpClippedSegSam, tmpMainBodySam_1, tmpMainBodySam_2;
		getline(sam_mainBody_clippedSeg_ifs, tmpClippedSegSam);
		if(tmpClippedSegSam == "")
			break;
		cout << "tmpClippedSegSam: " << tmpClippedSegSam << endl;
		getline(sam_mainBody_clippedSeg_ifs, tmpMainBodySam_1);
		getline(sam_mainBody_clippedSeg_ifs, tmpMainBodySam_2);
		cout << "tmpMainBodySam_1: " << tmpMainBodySam_1 << endl;
		cout << "tmpMainBodySam_2: " << tmpMainBodySam_2 << endl;

		int tmp_chrNameInt_mainBody, tmp_chrPos_mainBody, tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg;
		vector<Jump_Code> tmp_jumpCodeVec_mainBody, tmp_jumpCodeVec_clippedSeg;
		bool tmp_mainBody_forOrRevMap_bool, tmp_clippedSegAt_downstreamOrUpstream_end_bool, tmp_clippedSeg_forOrRevMap_bool;
		PE_Read_Info tmpPeReadInfo;

		bool parse_success_bool = parse_sam_clippedSeg_mainBody(tmpClippedSegSam, tmpMainBodySam_1, tmpMainBodySam_2, indexInfo,
			tmp_chrNameInt_mainBody, tmp_chrPos_mainBody, tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg,
			tmp_jumpCodeVec_mainBody, tmp_jumpCodeVec_clippedSeg,
			tmp_mainBody_forOrRevMap_bool, tmp_clippedSegAt_downstreamOrUpstream_end_bool, tmp_clippedSeg_forOrRevMap_bool, tmpPeReadInfo);
		//cout << "tmp_mainBody_forOrRevMap_bool: " << tmp_mainBody_forOrRevMap_bool << endl;
		//cout << "tmp_clippedSegAt_downstreamOrUpstream_end_bool: " << tmp_clippedSegAt_downstreamOrUpstream_end_bool << endl;
		//cout << "tmp_clippedSeg_forOrRevMap_bool: " << tmp_clippedSeg_forOrRevMap_bool << endl;

		if(parse_success_bool)
		{	
			FixFusion_PeRead tmpFixFusionPeReadInfo;
			tmpFixFusionPeReadInfo.initiate(tmp_chrNameInt_mainBody, tmp_chrPos_mainBody, tmp_jumpCodeVec_mainBody, 
				tmp_mainBody_forOrRevMap_bool, tmp_clippedSegAt_downstreamOrUpstream_end_bool, tmp_clippedSeg_forOrRevMap_bool, 
				tmp_chrNameInt_clippedSeg, tmp_chrPos_clippedSeg, tmp_jumpCodeVec_clippedSeg, tmpPeReadInfo);
			tmpFixFusionPeReadInfo.fixFusion(tmpPeReadInfo, tmpGeneAnnHashInfo, indexInfo, 
				geneAnnEntryBoundaryOnly_bool, geneAnnIncorporated_bool, insertionAllowed_bool,
				maximum_allowed_mismatchNum, maximum_allowed_insertionLength);
			bool tmpFusionFixed_success_bool = tmpFixFusionPeReadInfo.return_fusionFixed_bool();
			if(tmpFusionFixed_success_bool)
			{	
				int tmpFusionFixed_donerAnchorLength = tmpFixFusionPeReadInfo.return_fusionFixed_donerAnchorLength();
				int tmpFusionFixed_foreignInsertionLength = tmpFixFusionPeReadInfo.return_fusionFixed_foreignInsertionLength();
				int tmpFusionFixed_mismatchNum = tmpFixFusionPeReadInfo.return_fusionFixed_mismatchNum();
				cout << "tmpFusionFixed_donerAnchorLength: " << tmpFusionFixed_donerAnchorLength << endl;
				cout << "tmpFusionFixed_foreignInsertionLength: " << tmpFusionFixed_foreignInsertionLength << endl;
			}
		}
	}

	log_ofs.close();
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	return 0;
}