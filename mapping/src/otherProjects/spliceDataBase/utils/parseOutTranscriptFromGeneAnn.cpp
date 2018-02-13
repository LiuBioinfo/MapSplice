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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputGTF outputTranscriptFile" << endl;
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
	Index_Info* indexInfo = new Index_Info(indexParameterFileStr);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	chrom_bit_file_ifs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;
	cout << endl << "[" << asctime(local) << "... start to parse gene annotation file" << endl;

	string outputTranscriptFile = argv[3];
	// string outputTranscriptFile_tmp = outputTranscriptFile + ".tmp";
	// ofstream transcript_ofs_tmp(outputTranscriptFile_tmp.c_str());
	vector<string> chrNameVec_exon;
	vector<int> startPosVec_exon;
	vector<int> endPosVec_exon;
	vector<string> strandVec_exon;
	vector<string> transcriptIdVec_exon;
	vector<string> geneNameVec_exon;
	string inputGTFfile = argv[2];
	ifstream gtf_ifs(inputGTFfile.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == "#")
			continue;
		vector<string> tmpGtfFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc + 1);
			if(tabLoc == string::npos)
				break;
			string tmpGtfFieldStr = tmpStr.substr(startLoc, tabLoc - startLoc);
			tmpGtfFieldStrVec.push_back(tmpGtfFieldStr);
			startLoc = tabLoc + 1;
		}
		if(tmpGtfFieldStrVec.size() < 7)
			continue;
		string tmpChrNameStr = tmpGtfFieldStrVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		string tmpTypeStr = tmpGtfFieldStrVec[2];
		if(tmpTypeStr != "exon")
			continue;
		string tmpStartPosStr = tmpGtfFieldStrVec[3];
		string tmpEndPosStr = tmpGtfFieldStrVec[4];
		string tmpStrandStr = tmpGtfFieldStrVec[6];
		if((tmpChrNameInt < 0)||((tmpStrandStr != "+")&&(tmpStrandStr != "-")))
			continue;
		int tmpGeneName_loc = tmpStr.find("gene_name");
		if(tmpGeneName_loc == string::npos)
			continue;
		int next1stQuota_loc_geneName = tmpStr.find("\"", tmpGeneName_loc + 1);
		int next2ndQuota_loc_geneName = tmpStr.find("\"", next1stQuota_loc_geneName + 1);
		string tmpGeneNameStr = tmpStr.substr(next1stQuota_loc_geneName + 1, 
			next2ndQuota_loc_geneName - next1stQuota_loc_geneName - 1);
		int tmpTranscriptId_loc = tmpStr.find("transcript_id");
		if(tmpTranscriptId_loc == string::npos)
			continue;
		int next1stQuota_loc_transcriptId = tmpStr.find("\"", tmpTranscriptId_loc + 1);
		int next2ndQuota_loc_transcriptId = tmpStr.find("\"", next1stQuota_loc_transcriptId + 1);
		string tmpTranscriptIdStr = tmpStr.substr(next1stQuota_loc_transcriptId + 1, 
			next2ndQuota_loc_transcriptId - next1stQuota_loc_transcriptId - 1);
		//transcript_ofs_tmp << tmpChrNameStr << "\t" << tmpStartPosStr << "\t" << tmpEndPosStr << "\t" 
		//	<< tmpStrandStr << "\t" << tmpTranscriptIdStr << "\t" << tmpGeneNameStr << endl;
		chrNameVec_exon.push_back(tmpChrNameStr);
		int tmpStartPosInt = atoi(tmpStartPosStr.c_str());
		startPosVec_exon.push_back(tmpStartPosInt);
		int tmpEndPosInt = atoi(tmpEndPosStr.c_str());
		endPosVec_exon.push_back(tmpEndPosInt);
		strandVec_exon.push_back(tmpStrandStr);
		transcriptIdVec_exon.push_back(tmpTranscriptIdStr);
		geneNameVec_exon.push_back(tmpGeneNameStr);	
	}
	gtf_ifs.close();

	ofstream transcript_ofs(outputTranscriptFile.c_str());
	int totalExonNum = chrNameVec_exon.size();
	string tmpChrName_currentTranscript = chrNameVec_exon[0];
	vector<int> tmpStartPosVec_currentTranscript;
	tmpStartPosVec_currentTranscript.push_back(startPosVec_exon[0]);
	vector<int> tmpEndPosVec_currentTranscript;
	tmpEndPosVec_currentTranscript.push_back(endPosVec_exon[0]);
	string tmpStrand_currentTranscript = strandVec_exon[0];
	string tmpTranscriptId_currentTranscript = transcriptIdVec_exon[0];
	string tmpGeneName_currentTranscript = geneNameVec_exon[0];
	for(int tmp = 1; tmp < totalExonNum; tmp ++)
	{
		//vector<string> chrNameVec_exon;
		//vector<int> startPosVec_exon;
		//vector<int> endPosVec_exon;
		//vector<string> strandVec_exon;
		//vector<string> transcriptIdVec_exon;
		//vector<string> exonNameVec_exon;
		string tmpChrName_tmpExon = chrNameVec_exon[tmp];
		int tmpStartPos_tmpExon = startPosVec_exon[tmp];
		int tmpEndPos_tmpExon = endPosVec_exon[tmp];
		string tmpStrand_tmpExon = strandVec_exon[tmp];
		string tmpTranscriptId_tmpExon = transcriptIdVec_exon[tmp];
		string tmpGeneName_tmpExon = geneNameVec_exon[tmp];
		if(tmpTranscriptId_tmpExon != tmpTranscriptId_currentTranscript) // new transcript 
		{
			// output last transcript
			sort(tmpStartPosVec_currentTranscript.begin(), tmpStartPosVec_currentTranscript.end());
			sort(tmpEndPosVec_currentTranscript.begin(), tmpEndPosVec_currentTranscript.end());
			int tmpExonNum_currentTranscript = tmpStartPosVec_currentTranscript.size();
			transcript_ofs << tmpChrName_currentTranscript << "\t" << tmpExonNum_currentTranscript << "\t";
			for(int tmpExonIndex = 0; tmpExonIndex < tmpExonNum_currentTranscript; tmpExonIndex ++)
				transcript_ofs << tmpStartPosVec_currentTranscript[tmpExonIndex] << ",";
			transcript_ofs << "\t";
			for(int tmpExonIndex = 0; tmpExonIndex < tmpExonNum_currentTranscript; tmpExonIndex ++)
				transcript_ofs << tmpEndPosVec_currentTranscript[tmpExonIndex] << ",";
			transcript_ofs << "\t" << tmpStrand_currentTranscript << "\t" << tmpTranscriptId_currentTranscript
				<< "\t" <<  tmpGeneName_currentTranscript << endl;

			tmpChrName_currentTranscript = tmpChrName_tmpExon;
			vector<int>().swap(tmpStartPosVec_currentTranscript);
			tmpStartPosVec_currentTranscript.push_back(tmpStartPos_tmpExon);
			vector<int>().swap(tmpEndPosVec_currentTranscript);
			tmpEndPosVec_currentTranscript.push_back(tmpEndPos_tmpExon);
			tmpStrand_currentTranscript = tmpStrand_tmpExon;
			tmpTranscriptId_currentTranscript = tmpTranscriptId_tmpExon;
			tmpGeneName_currentTranscript = tmpGeneName_tmpExon;
		}
		else // the same transcript as the current one
		{
			if((tmpChrName_currentTranscript != tmpChrName_tmpExon)
				||(tmpGeneName_currentTranscript != tmpGeneName_tmpExon)
				||(tmpStrand_currentTranscript != tmpStrand_tmpExon))
			{
				cout << "tmpChrName_currentTranscript: " << tmpChrName_currentTranscript << endl;
				cout << "tmpChrName_tmpExon: " << tmpChrName_tmpExon << endl;
				cout << "tmpGeneName_currentTranscript: " << tmpGeneName_currentTranscript << endl;
				cout << "tmpGeneName_tmpExon: " << tmpGeneName_tmpExon << endl;
				cout << "tmpStrand_currentTranscript: " << tmpStrand_currentTranscript << endl;
				cout << "tmpStrand_tmpExon: " << tmpStrand_tmpExon << endl;
				exit(1);
			}
			tmpStartPosVec_currentTranscript.push_back(tmpStartPos_tmpExon);
			tmpEndPosVec_currentTranscript.push_back(tmpEndPos_tmpExon);
		}
		if(tmp == totalExonNum - 1) // the last transcript
		{
			sort(tmpStartPosVec_currentTranscript.begin(), tmpStartPosVec_currentTranscript.end());
			sort(tmpEndPosVec_currentTranscript.begin(), tmpEndPosVec_currentTranscript.end());
			int tmpExonNum_currentTranscript = tmpStartPosVec_currentTranscript.size();
			transcript_ofs << tmpChrName_currentTranscript << "\t" << tmpExonNum_currentTranscript << "\t";
			for(int tmpExonIndex = 0; tmpExonIndex < tmpExonNum_currentTranscript; tmpExonIndex ++)
				transcript_ofs << tmpStartPosVec_currentTranscript[tmpExonIndex] << ",";
			transcript_ofs << "\t";
			for(int tmpExonIndex = 0; tmpExonIndex < tmpExonNum_currentTranscript; tmpExonIndex ++)
				transcript_ofs << tmpEndPosVec_currentTranscript[tmpExonIndex] << ",";
			transcript_ofs << "\t" << tmpStrand_currentTranscript << "\t" << tmpTranscriptId_currentTranscript
				<< "\t" <<  tmpGeneName_currentTranscript << endl;			
		}
	}
	transcript_ofs.close();
	delete indexInfo;
	free(chrom);
	return 0;
}