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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../general/flkStrChangedSJ_info.h"

using namespace std;

bool isCanSJ(string& tmpFlkStr)	
{
	if((tmpFlkStr == "GTAG")||(tmpFlkStr == "CTAC")
		||(tmpFlkStr == "ATAC")||(tmpFlkStr == "GTAT")
		||(tmpFlkStr == "GCAG")||(tmpFlkStr == "CTGC"))
		return true;
	else
		return false;
}

int canSJhapNum(string& tmpFlkStr_hap1, string& tmpFlkStr_hap2)
{
	int tmpCanSJhapNum = 0;
	if(isCanSJ(tmpFlkStr_hap1))
		tmpCanSJhapNum ++;
	if(isCanSJ(tmpFlkStr_hap2))
		tmpCanSJhapNum ++;	
	return tmpCanSJhapNum;
}

void summarizeCanSJhapSampleSupNumPair(vector< pair<int,int> >& tmpCanSJhapNumAndSupNumPairVec, 
	int& tmpCanSJ_nonHap_sampleNum, int& tmpCanSJ_nonHap_supNum, 
	int& tmpCanSJ_existHap_sampleNum, int& tmpCanSJ_existHap_supNum, 
	int& tmpCanSJ_singleHap_sampleNum, int& tmpCanSJ_singleHap_supNum,
	int& tmpCanSJ_bothHap_sampleNum, int& tmpCanSJ_bothHap_supNum)
{
	tmpCanSJ_nonHap_sampleNum = 0; 
	tmpCanSJ_nonHap_supNum = 0;
	tmpCanSJ_existHap_sampleNum = 0;
	tmpCanSJ_existHap_supNum = 0;
	tmpCanSJ_singleHap_sampleNum = 0;
	tmpCanSJ_singleHap_supNum = 0;
	tmpCanSJ_bothHap_sampleNum = 0;
	tmpCanSJ_bothHap_supNum	= 0;
	for(int tmp = 0; tmp < tmpCanSJhapNumAndSupNumPairVec.size(); tmp++)
	{
		int tmpCanHapNum = tmpCanSJhapNumAndSupNumPairVec[tmp].first;
		int tmpSupNum = tmpCanSJhapNumAndSupNumPairVec[tmp].second;
		if(tmpCanHapNum == 0)
		{
			tmpCanSJ_nonHap_sampleNum ++;
			tmpCanSJ_nonHap_supNum += tmpSupNum;
		}
		else if(tmpCanHapNum == 1)
		{
			tmpCanSJ_singleHap_sampleNum ++;
			tmpCanSJ_singleHap_supNum += tmpSupNum;
			tmpCanSJ_existHap_sampleNum ++;
			tmpCanSJ_existHap_supNum += tmpSupNum;
		}
		else if(tmpCanHapNum == 2)
		{
			tmpCanSJ_bothHap_sampleNum ++;
			tmpCanSJ_bothHap_supNum += tmpSupNum;
			tmpCanSJ_existHap_sampleNum ++;
			tmpCanSJ_existHap_supNum += tmpSupNum;
		}
		else
		{
			cout << "invalid tmpCanHapNum: " << tmpCanHapNum << endl;
			exit(1);
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout <<  "#1 index_info" << endl;
		cout << "#2 changedSJlist" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();
	parameter_ifs.close();

	string changedSJfileList = argv[2];
	string outputDir = argv[3]; outputDir += "/";
	string mkdir_output = "mkdir -p " + outputDir;
	system(mkdir_output.c_str());
	string log_file = outputDir + "log";
	ofstream log_ofs(log_file.c_str());

	cout << "start to read changedSJfileVec..." << endl;
	vector<string> changedSJfileVec;
	ifstream changedSJfileList_ifs(changedSJfileList.c_str());
	while(!changedSJfileList_ifs.eof())
	{
		string tmpStr;
		getline(changedSJfileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		changedSJfileVec.push_back(tmpStr);
	}
	changedSJfileList_ifs.close();

	int fileNum = changedSJfileVec.size();
	cout << "fileNum: " << fileNum << endl;
	vector<FlkStrChangedSJ_Hash_Info> flkStrChangedSJhashInfoVec;
	FlkStrChangedSJ_Hash_Info flkStrChangedSJhashInfo_merged;
	flkStrChangedSJhashInfo_merged.initiate_chrNum(chromNum);
	cout<< "start to load each flkStrChangedSJ_file" << endl;
	for(int tmp = 0; tmp < fileNum; tmp ++)
	{
		flkStrChangedSJhashInfo_merged.loadFlkStrChangedSJfile(changedSJfileVec[tmp], indexInfo);
		FlkStrChangedSJ_Hash_Info tmpFlkStrChangedSJhashInfo;
		tmpFlkStrChangedSJhashInfo.initiate_chrNum(chromNum);
		tmpFlkStrChangedSJhashInfo.loadFlkStrChangedSJfile(changedSJfileVec[tmp], indexInfo);
		flkStrChangedSJhashInfoVec.push_back(tmpFlkStrChangedSJhashInfo);
	}

	string output_non2canSJ = outputDir + "non2can.txt";
	ofstream non2can_ofs(output_non2canSJ.c_str());
	int totalFlkStrChangedSJnum = flkStrChangedSJhashInfo_merged.return_flkStrChangedSJnum();
	cout << "totalFlkStrChangedSJnum: " << totalFlkStrChangedSJnum << endl;
	for(int tmp = 0; tmp < totalFlkStrChangedSJnum; tmp++)
	{
		int tmpChrNameInt = flkStrChangedSJhashInfo_merged.return_chrNameInt_vecIndex(tmp);
		string tmpChrName = indexInfo->returnChrNameStr(tmpChrNameInt);
		int tmpDonerEndPos = flkStrChangedSJhashInfo_merged.return_donerEndPos_vecIndex(tmp);
		int tmpAcceptorStartPos = flkStrChangedSJhashInfo_merged.return_acceptorStartPos_vecIndex(tmp);
		int tmpSupNum = flkStrChangedSJhashInfo_merged.return_supNum_vecIndex(tmp);
		string tmpFlkStr_ref = flkStrChangedSJhashInfo_merged.return_flkStr_ref_vecIndex(tmp);
		if(isCanSJ(tmpFlkStr_ref))
			continue;
		vector< pair<int,int> > tmpCanSJhapNumAndSupNumPairVec;
		for(int tmp2 = 0; tmp2 < fileNum; tmp2++)
		{
			int tmpCanSJhapNum, tmpSupNum;
			flkStrChangedSJhashInfoVec[tmp2].search_and_return_canSJhapNum_supNum(
				tmpCanSJhapNum, tmpSupNum, tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
			tmpCanSJhapNumAndSupNumPairVec.push_back(pair<int,int>(tmpCanSJhapNum, tmpSupNum));
		}
		int tmpCanSJ_nonHap_sampleNum = 0, tmpCanSJ_existHap_sampleNum = 0, tmpCanSJ_singleHap_sampleNum = 0, tmpCanSJ_bothHap_sampleNum;
		int tmpCanSJ_nonHap_supNum = 0, tmpCanSJ_existHap_supNum = 0, tmpCanSJ_singleHap_supNum = 0, tmpCanSJ_bothHap_supNum;
		summarizeCanSJhapSampleSupNumPair(tmpCanSJhapNumAndSupNumPairVec, 
			tmpCanSJ_nonHap_sampleNum, tmpCanSJ_nonHap_supNum, 
			tmpCanSJ_existHap_sampleNum, tmpCanSJ_existHap_supNum, 
			tmpCanSJ_singleHap_sampleNum, tmpCanSJ_singleHap_supNum,
			tmpCanSJ_bothHap_sampleNum, tmpCanSJ_bothHap_supNum);
		if(tmpCanSJ_existHap_sampleNum > 0)
			non2can_ofs << tmpChrName << "\t" << tmpDonerEndPos << "\t" << tmpAcceptorStartPos << "\t" << tmpFlkStr_ref 
				<< "\t" << tmpCanSJ_nonHap_sampleNum << "\t" << tmpCanSJ_nonHap_supNum 
				<< "\t" << tmpCanSJ_existHap_sampleNum << "\t" << tmpCanSJ_existHap_supNum
				<< "\t" << tmpCanSJ_singleHap_sampleNum << "\t" << tmpCanSJ_singleHap_supNum 
				<< "\t" << tmpCanSJ_bothHap_sampleNum << "\t" << tmpCanSJ_bothHap_supNum << endl;
	}
	non2can_ofs.close();
	delete indexInfo;
	log_ofs.close();
	return 0;
}