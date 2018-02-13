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
#include "general/alu_entry_hash.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputAluAnnFile queryChrName queryChrPos" << endl;
		exit(1);
	}
	cout << "initiate indexInfo ..." << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;

	string inputAluAnnEntryFile = argv[2];
	cout << "start to load salu annotation file" << endl;
	Alu_Entry_Hash tmpAluAnnHashInfo;
	tmpAluAnnHashInfo.initiate_aluAnnEntryArea2infoIndexMapVec(chromNum);
	tmpAluAnnHashInfo.loadAluAnn(inputAluAnnEntryFile, indexInfo);

	cout << "start to do query" << endl;
	string quertChrName = argv[3];
	string queryChrPosStr = argv[4];
	int queryChrPos = atoi(queryChrPosStr.c_str());
	vector<string> foundAluAnnEntryVec;
	tmpAluAnnHashInfo.searchAndReturnAluAnnEntryStrVec(foundAluAnnEntryVec, quertChrName, queryChrPos, indexInfo);

	int foundAluAnnEntryVecSize = foundAluAnnEntryVec.size();
	cout << "In total, " << foundAluAnnEntryVecSize << " entries are found in alu annotation file." << endl;
	for(int tmp = 0; tmp < foundAluAnnEntryVecSize; tmp++)
	{
		string tmpEntry = foundAluAnnEntryVec[tmp];
		cout << "Entry " << tmp + 1 << ": " << tmpEntry << endl;
	}

	delete indexInfo;
	parameter_ifs.close();
	return 0;
}