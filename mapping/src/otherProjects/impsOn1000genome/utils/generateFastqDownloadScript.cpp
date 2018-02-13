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

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 individualIdList" << endl;
		cout << "#2 individualId2fastqSiteList" << endl;
		cout << "#3 fastqDataDir" << endl;
		cout << "#4 downloadScript" << endl;
		exit(1);
	}
	string toDownloadIndividualIdList = argv[1];
	string individualId2fastqSiteList = argv[2];
	string fastqDataDir = argv[3];
	string downloadScript = argv[4];

	vector<string> toDownloadIdVec;
	ifstream toDownloadId_ifs(toDownloadIndividualIdList.c_str());
	while(!toDownloadId_ifs.eof())
	{
		string tmpStr;
		getline(toDownloadId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		toDownloadIdVec.push_back(tmpStr);
	}
	toDownloadId_ifs.close();
	
	vector<string> idVec;
	vector<string> fqVec_1;
	vector<string> fqVec_2;
	ifstream id2fq_ifs(individualId2fastqSiteList.c_str());
	while(!id2fq_ifs.eof())
	{
		string tmpStr_1;
		getline(id2fq_ifs, tmpStr_1);
		if(tmpStr_1 == "")
			break;
		string tmpStr_2;
		getline(id2fq_ifs, tmpStr_2);
		if(tmpStr_2 == "")
			break;		
		int tabLoc_1 = tmpStr_1.find("\t");
		int tabLoc_2 = tmpStr_2.find("\t");
		string tmpId_1 = tmpStr_1.substr(0, tabLoc_1);
		string tmpFq_1 = tmpStr_1.substr(tabLoc_1 + 1);
		string tmpId_2 = tmpStr_2.substr(0, tabLoc_2);
		string tmpFq_2 = tmpStr_2.substr(tabLoc_2 + 1);		
		if(tmpId_1 != tmpId_2)
		{
			cout << "tmpId_1 != tmpId_2" << endl;
			exit(1);
		}
		idVec.push_back(tmpId_1);
		fqVec_1.push_back(tmpFq_1);
		fqVec_2.push_back(tmpFq_2);
	}	
	id2fq_ifs.close();

	ofstream downloadScript_ofs(downloadScript.c_str());
	for(int tmp = 0; tmp < toDownloadIdVec.size(); tmp++)
	{
		string tmpToDownloadId = toDownloadIdVec[tmp];
		string tmpToDownload_fq_1, tmpToDownload_fq_2;
		bool tmpIdFound_bool = false;
		for(int tmp2 = 0; tmp2 < idVec.size(); tmp2++)
		{
			if(idVec[tmp2] == tmpToDownloadId)
			{
				tmpToDownload_fq_1 = fqVec_1[tmp2];
				tmpToDownload_fq_2 = fqVec_2[tmp2];
				tmpIdFound_bool = true;
				break;
			}
		}
		if(!tmpIdFound_bool)
		{
			cout << "!tmpIdFound_bool" << endl;
			exit(1);
		}
		downloadScript_ofs << "cd " << fastqDataDir << "/" << tmpToDownloadId << "/" << endl;
		downloadScript_ofs << "rm " << fastqDataDir << "/" << tmpToDownloadId << "/" << tmpToDownloadId << "_1.fastq" << endl;
		downloadScript_ofs << "rm " << fastqDataDir << "/" << tmpToDownloadId << "/" << tmpToDownloadId << "_2.fastq" << endl;
		downloadScript_ofs << "wget " << tmpToDownload_fq_1 << endl;
		downloadScript_ofs << "wget " << tmpToDownload_fq_2 << endl;
		downloadScript_ofs << endl;		
	}
	downloadScript_ofs.close();
	return 0;
}