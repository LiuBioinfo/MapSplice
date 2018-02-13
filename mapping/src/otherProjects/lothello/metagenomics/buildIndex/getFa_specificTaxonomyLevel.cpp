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
		cout << "Executable fa2taxInfoFile outputDir taxonomyLevel" << endl;
		exit(1);
	}
	int taxonomyLevelInt = -1; // 3--phylum; 4--class; 5--order; 6-- family; 7--genus; 8--species 
	string taxonomyLevelStr = argv[3];
	if((taxonomyLevelStr == "Phylum")||(taxonomyLevelStr == "PHYLUM")||(taxonomyLevelStr == "phylum"))
		taxonomyLevelInt = 3;
	else if((taxonomyLevelStr == "Genus")||(taxonomyLevelStr == "GENUS")||(taxonomyLevelStr == "genus"))
		taxonomyLevelInt = 7;
	else if((taxonomyLevelStr == "Species")||(taxonomyLevelStr == "SPECIES")||(taxonomyLevelStr == "species"))
		taxonomyLevelInt = 8;
	else
	{
		cout << "invalid taxonomyLevel: " << taxonomyLevelStr << endl;
		exit(1);
	}

	string outputDir = argv[2];
	outputDir += "/";
	string cmd_mkdir_outputDir = "mkdir " + outputDir;
	system(cmd_mkdir_outputDir.c_str());

	vector<string> taxonomyIdVec;
	vector< vector<string> > taxonomyFaPathVecVec;

	int taxonomyIdFieldIndex;
	if(taxonomyLevelInt == 3) // phylum
		taxonomyIdFieldIndex = 7;
	else if(taxonomyLevelInt == 7) // genus
		taxonomyIdFieldIndex = 3;
	else if(taxonomyLevelInt == 8) // species
		taxonomyIdFieldIndex = 2;
	else
	{
		cout << "invalid taxonomyLevelInt: " << taxonomyLevelInt << endl;
		exit(1);
	}

	vector<string> rawFaPathVec;
	vector<string> rawTaxonomyIdVec;
	string fa2taxInfoFile = argv[1];
	ifstream fa2taxInfo_ifs(fa2taxInfoFile.c_str());
	while(!fa2taxInfo_ifs.eof())
	{
		string tmpStr;
		getline(fa2taxInfo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpFaPath = tmpFieldVec[0];
		string tmpTaxoId = tmpFieldVec[taxonomyIdFieldIndex];
		rawFaPathVec.push_back(tmpFaPath);
		rawTaxonomyIdVec.push_back(tmpTaxoId);
	}
	fa2taxInfo_ifs.close();
	cout << "rawFa #: " << rawFaPathVec.size() << endl;

	for(int tmp = 0; tmp < rawFaPathVec.size(); tmp++)
	{
		string tmpFaPath = rawFaPathVec[tmp];
		string tmpTaxoId = rawTaxonomyIdVec[tmp];
		int currentTaxonomyIdVecSize = taxonomyIdVec.size();
		bool newTaxoId_bool = true;
		for(int tmp2 = 0; tmp2 < currentTaxonomyIdVecSize; tmp2++)
		{
			string tmpExistingTaxoId = taxonomyIdVec[tmp2];
			if(tmpExistingTaxoId == tmpTaxoId)
			{
				taxonomyFaPathVecVec[tmp2].push_back(tmpFaPath);
				newTaxoId_bool = false;
				break;
			}
		}
		if(newTaxoId_bool)
		{
			taxonomyIdVec.push_back(tmpTaxoId);
			vector<string> tmpFaPathVec;
			tmpFaPathVec.push_back(tmpFaPath);
			taxonomyFaPathVecVec.push_back(tmpFaPathVec);
		}
	}
	cout << "taxonomyIdVec.size(): " << taxonomyIdVec.size() << endl;
	// merge Fa vec at specific taxo level
	for(int tmp = 0; tmp < taxonomyIdVec.size(); tmp++)
	{
		string cmd_cat = "cat";
		for(int tmp2 = 0; tmp2 < taxonomyFaPathVecVec[tmp].size(); tmp2 ++)
		{
			cmd_cat += " ";
			cmd_cat += (taxonomyFaPathVecVec[tmp])[tmp2];
		}
		cmd_cat += " > ";
		cmd_cat += outputDir;
		cmd_cat += taxonomyIdVec[tmp];
		cmd_cat += ".fa";
		cout << "tmp merge fa cmd: " << endl << cmd_cat << endl;
		system(cmd_cat.c_str());
	}
	return 0;
}