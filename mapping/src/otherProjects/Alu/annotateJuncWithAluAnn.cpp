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
		cout << "Executable inputIndexPath inputAluAnnFile inputJuncFile outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	string outputDirStr = outputFolderStr + "/";
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputJuncFile_AnnotatedWithAluAnn = outputDirStr + "/junc_annotatedWithAlu.txt";
	string alu_valid_file = outputDirStr + "alu_valid.txt";
	string alu_invalid_file = outputDirStr + "alu_invalid.txt";
	//ofstream alu_valid_ofs(alu_valid_file.c_str());
	//ofstream alu_invalid_ofs(alu_invalid_file.c_str());
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
	tmpAluAnnHashInfo.loadAluAnn(inputAluAnnEntryFile, indexInfo, alu_valid_file, alu_invalid_file);

	cout << "start to scan junc file with alu ann" << endl;
	string inputJuncFile = argv[3];
	//string outputJuncFile_AnnotatedWithAluAnn = argv[4];
	ifstream junc_ifs(inputJuncFile.c_str());
	ofstream annotatedJunc_ofs(outputJuncFile_AnnotatedWithAluAnn.c_str());
	while(!junc_ifs.eof())
	{
		string tmpStr;
		getline(junc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//cout<<  tmpStr << endl;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpJunc_chrName = tmpStr.substr(0, tabLoc_1);
		string tmpJunc_startPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpJunc_endPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		// cout << "tmpJunc_chrName: " << tmpJunc_chrName << endl;
		// cout << "tmpJunc_startPos: " << tmpJunc_startPos << endl;
		// cout << "tmpJunc_endPos: " << tmpJunc_endPos << endl;
		vector<string> tmpJunc_startPos_aluVec;
		vector<string> tmpJunc_endPos_aluVec;
		tmpAluAnnHashInfo.searchAndReturnAluGeneIdVec(tmpJunc_startPos_aluVec, tmpJunc_chrName, 
			tmpJunc_startPos, indexInfo);
		tmpAluAnnHashInfo.searchAndReturnAluGeneIdVec(tmpJunc_endPos_aluVec, tmpJunc_chrName, 
			tmpJunc_endPos, indexInfo);
		//cout << "tmpJunc_startPos_aluVec.size(): " << tmpJunc_startPos_aluVec.size() << endl;
		//cout << "tmpJunc_endPos_aluVec.size(): " << tmpJunc_endPos_aluVec.size() << endl;
		if((tmpJunc_startPos_aluVec.size() > 1)||(tmpJunc_endPos_aluVec.size() > 1))
		{
			annotatedJunc_ofs << tmpStr;		
			if(tmpJunc_startPos_aluVec.size() > 1)
				annotatedJunc_ofs << "\t" << tmpJunc_startPos_aluVec[0];
			else
				annotatedJunc_ofs << "\tNULL";
			if(tmpJunc_endPos_aluVec.size() > 1)
				annotatedJunc_ofs << "\t" << tmpJunc_endPos_aluVec[0] << endl;
			else
				annotatedJunc_ofs << "\tNULL" << endl;
		}
	}

	// string quertChrName = argv[3];
	// string queryChrPosStr = argv[4];
	// int queryChrPos = atoi(queryChrPosStr.c_str());
	// vector<string> foundAluAnnEntryVec;
	// tmpAluAnnHashInfo.searchAndReturnAluAnnEntryStrVec(foundAluAnnEntryVec, quertChrName, queryChrPos, indexInfo);

	// int foundAluAnnEntryVecSize = foundAluAnnEntryVec.size();
	// cout << "In total, " << foundAluAnnEntryVecSize << " entries are found in alu annotation file." << endl;
	// for(int tmp = 0; tmp < foundAluAnnEntryVecSize; tmp++)
	// {
	// 	string tmpEntry = foundAluAnnEntryVec[tmp];
	// 	cout << "Entry " << tmp + 1 << ": " << tmpEntry << endl;
	// }
	//alu_valid_ofs.close();
	//alu_invalid_ofs.close();
	junc_ifs.close();
	annotatedJunc_ofs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}