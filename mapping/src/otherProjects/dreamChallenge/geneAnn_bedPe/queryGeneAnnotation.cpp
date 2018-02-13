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
#include "../general/geneAnnEntryHash.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputGeneAnnEntryFile queryChrName queryChrPos" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	//string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	//ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;
	//char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	//chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	//indexInfo->readGenome(chrom);
	indexInfo->initiate_withoutLoadingSeq();
	//cout << "start to initiateChrNameIndexArray " << endl;
	//indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputGeneAnnEntryFile = argv[2];
	cout << "start to load simplified gene annotation file" << endl;
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	cout << "start to initiate GeneAnnEntryArea2infoIndexMapVec" << endl;
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);

	cout << "start to do query" << endl;
	string quertChrName = argv[3];
	string queryChrPosStr = argv[4];
	int queryChrPos = atoi(queryChrPosStr.c_str());
	vector<string> foundGeneAnnEntryVec;
	tmpGeneAnnHashInfo.searchAndReturnGeneAnnEntryStrVec(foundGeneAnnEntryVec, quertChrName, queryChrPos, indexInfo);

	int foundGeneAnnEntryVecSize = foundGeneAnnEntryVec.size();
	cout << "In total, " << foundGeneAnnEntryVecSize << " entries are found in gene annotation file." << endl;
	for(int tmp = 0; tmp < foundGeneAnnEntryVecSize; tmp++)
	{
		string tmpEntry = foundGeneAnnEntryVec[tmp];
		cout << "Entry " << tmp + 1 << ": " << tmpEntry << endl;
	}

	delete indexInfo;
	parameter_ifs.close();
	//chrom_bit_file_ifs.close();
	//free(chrom);
	return 0;
}