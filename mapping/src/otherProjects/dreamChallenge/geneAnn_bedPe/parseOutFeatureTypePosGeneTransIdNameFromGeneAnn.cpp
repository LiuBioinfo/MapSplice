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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexPath inputGeneAnnGTFfile outputFilePrefix" << endl;
		exit(1);
	}	
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputGeneAnnGTFfile = argv[2];
	string outputFilePrefix = argv[3];
	string outputFilePrefix_valid = outputFilePrefix + ".geneAnnEntry.valid";
	string outputFilePrefix_invalid = outputFilePrefix + ".geneAnnEntry.invalid";
	ofstream geneAnnEntry_valid_ofs(outputFilePrefix_valid.c_str());
	ofstream geneAnnEntry_invalid_ofs(outputFilePrefix_invalid.c_str());
	ifstream gtf_ifs(inputGeneAnnGTFfile.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpGTFstr;
		getline(gtf_ifs, tmpGTFstr);
		if(tmpGTFstr == "")
			break;
		if(tmpGTFstr.at(0) == '#')
		{
			geneAnnEntry_invalid_ofs << tmpGTFstr << endl;
			continue;
		}
		vector<string> tmpGTFfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 8; tmp++)
		{
			int tabLoc = tmpGTFstr.find("\t", startLoc);
			string tmpGTFfield = tmpGTFstr.substr(startLoc, tabLoc-startLoc);
			tmpGTFfieldVec.push_back(tmpGTFfield);
			startLoc = tabLoc + 1;
		}
		string tmpGTFlastField = tmpGTFstr.substr(startLoc);
		tmpGTFfieldVec.push_back(tmpGTFlastField);
		string tmpChrName = tmpGTFfieldVec[0];
		if(tmpChrName.at(0) != 'c')
			tmpChrName = "chr" + tmpChrName;
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			geneAnnEntry_invalid_ofs << tmpGTFstr << endl;
			continue;
		}
		string tmpSource = tmpGTFfieldVec[1];
		string tmpFeatureType = tmpGTFfieldVec[2];
		string tmpStartPosStr = tmpGTFfieldVec[3];
		string tmpEndPosStr = tmpGTFfieldVec[4];
		// int tmpStartPos = atoi(tmpStartPosStr.c_str());
		// int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpGTFstrand = tmpGTFfieldVec[6];
		string tmpOtherStr = tmpGTFfieldVec[8];
		string tmpGeneId, tmpGeneName, tmpTranscriptId, tmpTranscriptName;
		// gene id
		int tmpGeneId_loc = tmpGTFstr.find("gene_id");
		if(tmpGeneId_loc == string::npos)
			tmpGeneId == "NULL";
		else
		{
			int next1stQuota_loc_geneId = tmpGTFstr.find("\"", tmpGeneId_loc + 1);
			int next2ndQuota_loc_geneId = tmpGTFstr.find("\"", next1stQuota_loc_geneId + 1);
			tmpGeneId = tmpGTFstr.substr(next1stQuota_loc_geneId + 1, next2ndQuota_loc_geneId - next1stQuota_loc_geneId - 1);
		}		
		// gene name
		int tmpGeneName_loc = tmpGTFstr.find("gene_name");
		if(tmpGeneName_loc == string::npos)
			tmpGeneName == "NULL";
		else
		{
			int next1stQuota_loc_geneName = tmpGTFstr.find("\"", tmpGeneName_loc + 1);
			int next2ndQuota_loc_geneName = tmpGTFstr.find("\"", next1stQuota_loc_geneName + 1);
			tmpGeneName = tmpGTFstr.substr(next1stQuota_loc_geneName + 1, next2ndQuota_loc_geneName - next1stQuota_loc_geneName - 1);
		}
		// transcript id
		int tmpTranscriptId_loc = tmpGTFstr.find("transcript_id");
		if(tmpTranscriptId_loc == string::npos)
			tmpTranscriptId = "NULL";
		else
		{
			int next1stQuota_loc_transcriptId = tmpGTFstr.find("\"", tmpTranscriptId_loc + 1);
			int next2ndQuota_loc_transcriptId = tmpGTFstr.find("\"", next1stQuota_loc_transcriptId + 1);
			tmpTranscriptId = tmpGTFstr.substr(next1stQuota_loc_transcriptId + 1, next2ndQuota_loc_transcriptId - next1stQuota_loc_transcriptId - 1);
		}
		// transcript name
		int tmpTranscriptName_loc = tmpGTFstr.find("transcript_name");
		if(tmpTranscriptName_loc == string::npos)
			tmpTranscriptName = "NULL";
		else
		{
			int next1stQuota_loc_transcriptName = tmpGTFstr.find("\"", tmpTranscriptName_loc + 1);
			int next2ndQuota_loc_transcriptName = tmpGTFstr.find("\"", next1stQuota_loc_transcriptName + 1);
			tmpTranscriptName = tmpGTFstr.substr(next1stQuota_loc_transcriptName + 1, next2ndQuota_loc_transcriptName - next1stQuota_loc_transcriptName - 1);		
		}
		geneAnnEntry_valid_ofs << tmpChrName << "\t" << tmpStartPosStr << "\t" << tmpEndPosStr << "\t" << tmpGTFstrand << "\t" << tmpFeatureType 
			<< "\t" << tmpGeneId << "\t" << tmpGeneName << "\t" << tmpTranscriptId << "\t" << tmpTranscriptName << "\t" << tmpSource  << endl;
	}
	gtf_ifs.close();
	geneAnnEntry_valid_ofs.close();
	geneAnnEntry_invalid_ofs.close();
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	return 0;
}