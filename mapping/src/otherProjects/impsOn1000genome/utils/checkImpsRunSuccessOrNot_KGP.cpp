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

using namespace std;

void getIdVec(string& dataListFile, vector<string>& idVec)
{
	ifstream id_ifs(dataListFile.c_str());
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		if(tabLoc == string::npos)
			idVec.push_back(tmpStr);
		else
			idVec.push_back(tmpStr.substr(0, tabLoc));
	}
	id_ifs.close();
}

bool impsRunSuccessOrNot(string& tmpLogFile)
{
	bool runSuccessOrNotBool = false;
	ifstream tmpLog_ifs(tmpLogFile.c_str());
	if((!tmpLog_ifs.is_open())||(tmpLog_ifs.fail()))
		return false;
	vector<string> tmpStrVec;
	while(!tmpLog_ifs.eof())
	{
		string tmpStr;
		getline(tmpLog_ifs, tmpStr);
		tmpStrVec.push_back(tmpStr);
	}
	tmpLog_ifs.close();
	cout << "tmpStrVec.size(): " << tmpStrVec.size() << endl;
	int tmpStrvecSize = tmpStrVec.size();
	if(tmpStrvecSize < 4)
		return false;
	string tmpLastLine_1 = tmpStrVec[tmpStrvecSize-1];
	string tmpLastLine_2 = tmpStrVec[tmpStrvecSize-2];
	string tmpLastLine_3 = tmpStrVec[tmpStrvecSize-3];
	string tmpLastLine_4 = tmpStrVec[tmpStrvecSize-4];
	cout << "tmpLastLine_1: " << tmpLastLine_1 << endl;
	cout << "tmpLastLine_2: " << tmpLastLine_2 << endl;
	cout << "tmpLastLine_3: " << tmpLastLine_3 << endl;
	cout << "tmpLastLine_4: " << tmpLastLine_4 << endl;
	if((tmpLastLine_1.find("all jobs done") != string::npos)
		||(tmpLastLine_2.find("all jobs done") != string::npos)
		||(tmpLastLine_3.find("all jobs done") != string::npos)
		||(tmpLastLine_4.find("all jobs done") != string::npos))
		return true;
	else
		return false;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 dataList" << endl;
		cout << "#2 outputResults" << endl;
		exit(1);
	}
	string dataList = argv[1];
	string outputResults = argv[2];

	vector<string> idVec;
	getIdVec(dataList, idVec);

	ofstream results_ofs(outputResults.c_str());
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		results_ofs << tmpId << "\t";
		cout << "tmpId: " << tmpId << endl;
		string tmpOutputDir = "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/imps_twoHap_results/" + tmpId + "_impsTwoHapSingleIndex";
		string tmpLogFile = tmpOutputDir + "/runtime.log";
		cout << "tmpLogFile: " << tmpLogFile << endl;
		if(impsRunSuccessOrNot(tmpLogFile))
			results_ofs << "Yes" << endl;
		else
			results_ofs << "No" << endl;
	}
	results_ofs.close();
	return 0;
}