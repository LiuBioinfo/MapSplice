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

time_t nowtime;
struct tm *local;

bool reformatChrNameWithAccessionIdMap(string& raw_chrName, string& new_chrName, vector< pair<string, string> >& accessionIdMapVec)
{
	int accessionIdMapVecSize = accessionIdMapVec.size();
	for(int tmp = 0; tmp < accessionIdMapVecSize; tmp++)
	{
		string tmpAccessionId = accessionIdMapVec[tmp].first;
		if(raw_chrName == tmpAccessionId)
		{
			new_chrName = accessionIdMapVec[tmp].second;
			return true;
		}
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputAccessionIdMapFile inputGTF outputReformattedGTF" << endl;
		exit(1);
	}
	string inputAccessionIdMapFile = argv[1];
	vector< pair<string, string> > accessionIdMapVec; // <chrName, accessionId>
	ifstream accessionIdMap_ifs(inputAccessionIdMapFile.c_str());
	while(!accessionIdMap_ifs.eof())
	{
		string tmpStr;
		getline(accessionIdMap_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1st = tmpStr.find("\t");
		if(tabLoc_1st == string::npos)
			continue;
		string tmpChrName = tmpStr.substr(0, tabLoc_1st);
		string tmpAccessionId = tmpStr.substr(tabLoc_1st + 1);
		accessionIdMapVec.push_back(pair<string,string>(tmpAccessionId, tmpChrName));
	}
	accessionIdMap_ifs.close();
	string inputGTF = argv[2];
	string outputReformattedGTF = argv[3];
	ifstream gtf_ifs(inputGTF.c_str());
	ofstream reformattedGtf_ofs(outputReformattedGTF.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1st = tmpStr.find("\t");
		if(tabLoc_1st == string::npos)
			continue;
		string tmpChrStr_raw = tmpStr.substr(0, tabLoc_1st);
		string tmpGtfStr_others_raw = tmpStr.substr(tabLoc_1st + 1);
		string tmpChrStr_new;
		bool reformatChrName_success_bool = reformatChrNameWithAccessionIdMap(tmpChrStr_raw, tmpChrStr_new, accessionIdMapVec);
		if(reformatChrName_success_bool)
			reformattedGtf_ofs << tmpChrStr_new << "\t" << tmpGtfStr_others_raw << endl;
	}
	gtf_ifs.close();
	reformattedGtf_ofs.close();
	return 0;
}