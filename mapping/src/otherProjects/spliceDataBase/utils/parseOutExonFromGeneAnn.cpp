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
		cout << "Executable inputIndexFolderPath inputGTF outputExonFile" << endl;
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

	string outputExonFile = argv[3];
	ofstream exon_ofs(outputExonFile.c_str());
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
		int next1stQuota_loc = tmpStr.find("\"", tmpGeneName_loc + 1);
		int next2ndQuota_loc = tmpStr.find("\"", next1stQuota_loc + 1);
		string tmpGeneNameStr = tmpStr.substr(next1stQuota_loc + 1, next2ndQuota_loc - next1stQuota_loc - 1);
		exon_ofs << tmpChrNameStr << "\t" << tmpStartPosStr << "\t" << tmpEndPosStr << "\t" << tmpStrandStr << "\t" << tmpGeneNameStr << endl;
	}
	gtf_ifs.close();
	exon_ofs.close();
	delete  indexInfo;
	free(chrom);
	return 0;
}