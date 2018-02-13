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
	if((argc != 7)&&(argc != 6))
	{
		cout << "Executable inputIndexPath inputGeneAnnEntryFile inputGeneAnnEntryIndexFile queryChrName queryChrPos Mode_id(optional)" << endl;
		cout << "Allowed Mode_id: 0 -- within; 1 -- boundaryOnly" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;
	indexInfo->initiate_withoutLoadingSeq();
	cout << "end of initiating indexInfo" << endl;

	string inputGeneAnnEntryFile = argv[2];
	string inputGeneAnnEntryIndexFile = argv[3];
	cout << "start to load simplified gene annotation file" << endl;
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	cout << "start to initiate GeneAnnEntryInfoVec" << endl;
	tmpGeneAnnHashInfo.initiateGeneAnnEntryInfoVec(inputGeneAnnEntryFile, indexInfo);
	cout << "start to initiate GeneAnnEntryArea2infoIndexMapVec" << endl;
	tmpGeneAnnHashInfo.initiateGeneAnnEntryArea2infoIndexMapVec(inputGeneAnnEntryIndexFile, indexInfo);

	cout << "start to do query" << endl;
	string quertChrName = argv[4];
	string queryChrPosStr = argv[5];
	string modeIdStr;
	int modeId;
	if(argc == 7)
	{	
		modeIdStr = argv[6];
		modeId = atoi(modeIdStr.c_str());
	}
	else
		modeId = 0;

	int queryChrPos = atoi(queryChrPosStr.c_str());
	vector<string> foundGeneAnnEntryVec;
	tmpGeneAnnHashInfo.searchAndReturnGeneAnnEntryStrVec(foundGeneAnnEntryVec, quertChrName, queryChrPos, indexInfo, modeId);

	int foundGeneAnnEntryVecSize = foundGeneAnnEntryVec.size();
	cout << "In total, " << foundGeneAnnEntryVecSize << " entries are found in gene annotation file." << endl;
	for(int tmp = 0; tmp < foundGeneAnnEntryVecSize; tmp++)
	{
		string tmpEntry = foundGeneAnnEntryVec[tmp];
		cout << "Entry " << tmp + 1 << ": " << tmpEntry << endl;
	}

	delete indexInfo;
	parameter_ifs.close();
	return 0;
}