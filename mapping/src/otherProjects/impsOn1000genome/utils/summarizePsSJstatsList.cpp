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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

void get3valueFromStr(string& tmpStr, int& val, int& val_min3, int& val_min5)
{
	int tabLoc_1 = tmpStr.find("\t");
	int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
	string tmpValStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpValStr_min3 = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	string tmpValStr_min5 = tmpStr.substr(tabLoc_3 + 1);
	val = atoi(tmpValStr.c_str());
	val_min3 = atoi(tmpValStr_min3.c_str());
	val_min5 = atoi(tmpValStr_min5.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputPsSJstatsFileList" << endl;
		cout << "#2 outputSummarizedResults" << endl;
		exit(1);
	}
	string inputPsSJstatsFileList = argv[1];
	string outputSummarizedResults = argv[2];

	vector<string> psSJstatsFileVec;
	ifstream psSJstatsFileList_ifs(inputPsSJstatsFileList.c_str());
	while(!psSJstatsFileList_ifs.eof())
	{
		string tmpStr;
		getline(psSJstatsFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		psSJstatsFileVec.push_back(tmpStr);
	}
	psSJstatsFileList_ifs.close();

	vector<int> totalJuncNumVec;
	vector<int> canJuncNumVec;
	vector<int> semJuncNumVec;
	vector<int> nonJuncNumVec;
	vector<int> keptJuncNumVec;
	vector<int> changedJuncNumVec;
	vector<int> changedJuncNumVec_hom;
	vector<int> changedJuncNumVec_het;
	vector<int> changedJuncNumVec_can2can;
	vector<int> changedJuncNumVec_can2sem;
	vector<int> changedJuncNumVec_can2non;
	vector<int> changedJuncNumVec_sem2can;
	vector<int> changedJuncNumVec_sem2sem;
	vector<int> changedJuncNumVec_sem2non;
	vector<int> changedJuncNumVec_non2can;
	vector<int> changedJuncNumVec_non2sem;
	vector<int> changedJuncNumVec_non2non;
	vector<int> changedJuncNumVec_can2canCan;
	vector<int> changedJuncNumVec_can2canSem;
	vector<int> changedJuncNumVec_can2canNon;
	vector<int> changedJuncNumVec_can2semSem;
	vector<int> changedJuncNumVec_can2semNon;
	vector<int> changedJuncNumVec_can2nonNon;	
	vector<int> changedJuncNumVec_sem2canCan;
	vector<int> changedJuncNumVec_sem2canSem;
	vector<int> changedJuncNumVec_sem2canNon;
	vector<int> changedJuncNumVec_sem2semSem;
	vector<int> changedJuncNumVec_sem2semNon;
	vector<int> changedJuncNumVec_sem2nonNon;
	vector<int> changedJuncNumVec_non2canCan;
	vector<int> changedJuncNumVec_non2canSem;
	vector<int> changedJuncNumVec_non2canNon;
	vector<int> changedJuncNumVec_non2semSem;
	vector<int> changedJuncNumVec_non2semNon;
	vector<int> changedJuncNumVec_non2nonNon;

	vector<int> totalJuncNumVec_min3;
	vector<int> canJuncNumVec_min3;
	vector<int> semJuncNumVec_min3;
	vector<int> nonJuncNumVec_min3;
	vector<int> keptJuncNumVec_min3;
	vector<int> changedJuncNumVec_min3;
	vector<int> changedJuncNumVec_hom_min3;
	vector<int> changedJuncNumVec_het_min3;
	vector<int> changedJuncNumVec_can2can_min3;
	vector<int> changedJuncNumVec_can2sem_min3;
	vector<int> changedJuncNumVec_can2non_min3;
	vector<int> changedJuncNumVec_sem2can_min3;
	vector<int> changedJuncNumVec_sem2sem_min3;
	vector<int> changedJuncNumVec_sem2non_min3;
	vector<int> changedJuncNumVec_non2can_min3;
	vector<int> changedJuncNumVec_non2sem_min3;
	vector<int> changedJuncNumVec_non2non_min3;
	vector<int> changedJuncNumVec_can2canCan_min3;
	vector<int> changedJuncNumVec_can2canSem_min3;
	vector<int> changedJuncNumVec_can2canNon_min3;
	vector<int> changedJuncNumVec_can2semSem_min3;
	vector<int> changedJuncNumVec_can2semNon_min3;
	vector<int> changedJuncNumVec_can2nonNon_min3;	
	vector<int> changedJuncNumVec_sem2canCan_min3;
	vector<int> changedJuncNumVec_sem2canSem_min3;
	vector<int> changedJuncNumVec_sem2canNon_min3;
	vector<int> changedJuncNumVec_sem2semSem_min3;
	vector<int> changedJuncNumVec_sem2semNon_min3;
	vector<int> changedJuncNumVec_sem2nonNon_min3;
	vector<int> changedJuncNumVec_non2canCan_min3;
	vector<int> changedJuncNumVec_non2canSem_min3;
	vector<int> changedJuncNumVec_non2canNon_min3;
	vector<int> changedJuncNumVec_non2semSem_min3;
	vector<int> changedJuncNumVec_non2semNon_min3;
	vector<int> changedJuncNumVec_non2nonNon_min3;

	vector<int> totalJuncNumVec_min5;
	vector<int> canJuncNumVec_min5;
	vector<int> semJuncNumVec_min5;
	vector<int> nonJuncNumVec_min5;
	vector<int> keptJuncNumVec_min5;
	vector<int> changedJuncNumVec_min5;
	vector<int> changedJuncNumVec_hom_min5;
	vector<int> changedJuncNumVec_het_min5;
	vector<int> changedJuncNumVec_can2can_min5;
	vector<int> changedJuncNumVec_can2sem_min5;
	vector<int> changedJuncNumVec_can2non_min5;
	vector<int> changedJuncNumVec_sem2can_min5;
	vector<int> changedJuncNumVec_sem2sem_min5;
	vector<int> changedJuncNumVec_sem2non_min5;
	vector<int> changedJuncNumVec_non2can_min5;
	vector<int> changedJuncNumVec_non2sem_min5;
	vector<int> changedJuncNumVec_non2non_min5;
	vector<int> changedJuncNumVec_can2canCan_min5;
	vector<int> changedJuncNumVec_can2canSem_min5;
	vector<int> changedJuncNumVec_can2canNon_min5;
	vector<int> changedJuncNumVec_can2semSem_min5;
	vector<int> changedJuncNumVec_can2semNon_min5;
	vector<int> changedJuncNumVec_can2nonNon_min5;	
	vector<int> changedJuncNumVec_sem2canCan_min5;
	vector<int> changedJuncNumVec_sem2canSem_min5;
	vector<int> changedJuncNumVec_sem2canNon_min5;
	vector<int> changedJuncNumVec_sem2semSem_min5;
	vector<int> changedJuncNumVec_sem2semNon_min5;
	vector<int> changedJuncNumVec_sem2nonNon_min5;
	vector<int> changedJuncNumVec_non2canCan_min5;
	vector<int> changedJuncNumVec_non2canSem_min5;
	vector<int> changedJuncNumVec_non2canNon_min5;
	vector<int> changedJuncNumVec_non2semSem_min5;
	vector<int> changedJuncNumVec_non2semNon_min5;
	vector<int> changedJuncNumVec_non2nonNon_min5;

	for(int tmp = 0; tmp < psSJstatsFileVec.size(); tmp++)
	{
		string tmpStatsFile = psSJstatsFileVec[tmp];
		ifstream tmpStats_ifs(tmpStatsFile.c_str());
		int tmpLineNum = 0;
		vector<string> tmpLineVec;
		while(!tmpStats_ifs.eof())
		{
			tmpLineNum ++;
			if(tmpLineNum > 38)
				break;
			string tmpStr;
			getline(tmpStats_ifs, tmpStr);
			tmpLineVec.push_back(tmpStr);
		}
		tmpStats_ifs.close();
		string tmpStr_totalJunc = tmpLineVec[0];
		string tmpStr_canJunc = tmpLineVec[2];
		string tmpStr_semJunc = tmpLineVec[3];
		string tmpStr_nonJunc = tmpLineVec[4];
		string tmpStr_keptJunc = tmpLineVec[6];
		string tmpStr_changedJunc = tmpLineVec[8];
		string tmpStr_changedJunc_hom = tmpLineVec[9];
		string tmpStr_changedJunc_het = tmpLineVec[10];
		string tmpStr_changedJunc_can2can = tmpLineVec[11];
		string tmpStr_changedJunc_can2sem = tmpLineVec[12];
		string tmpStr_changedJunc_can2non = tmpLineVec[13];
		string tmpStr_changedJunc_sem2can = tmpLineVec[14];
		string tmpStr_changedJunc_sem2sem = tmpLineVec[15];
		string tmpStr_changedJunc_sem2non = tmpLineVec[16];
		string tmpStr_changedJunc_non2can = tmpLineVec[17];
		string tmpStr_changedJunc_non2sem = tmpLineVec[18];
		string tmpStr_changedJunc_non2non = tmpLineVec[19];
		string tmpStr_changedJunc_can2canCan = tmpLineVec[20];
		string tmpStr_changedJunc_can2canSem = tmpLineVec[21];
		string tmpStr_changedJunc_can2canNon = tmpLineVec[22];
		string tmpStr_changedJunc_can2semSem = tmpLineVec[23];
		string tmpStr_changedJunc_can2semNon = tmpLineVec[24];
		string tmpStr_changedJunc_can2nonNon = tmpLineVec[25];
		string tmpStr_changedJunc_sem2canCan = tmpLineVec[26];
		string tmpStr_changedJunc_sem2canSem = tmpLineVec[27];
		string tmpStr_changedJunc_sem2canNon = tmpLineVec[28];
		string tmpStr_changedJunc_sem2semSem = tmpLineVec[29];
		string tmpStr_changedJunc_sem2semNon = tmpLineVec[30];
		string tmpStr_changedJunc_sem2nonNon = tmpLineVec[31];
		string tmpStr_changedJunc_non2canCan = tmpLineVec[32];
		string tmpStr_changedJunc_non2canSem = tmpLineVec[33];
		string tmpStr_changedJunc_non2canNon = tmpLineVec[34];
		string tmpStr_changedJunc_non2semSem = tmpLineVec[35];
		string tmpStr_changedJunc_non2semNon = tmpLineVec[36];
		string tmpStr_changedJunc_non2nonNon = tmpLineVec[37];

		int tmp_totalJuncNum, tmp_canJuncNum, tmp_semJuncNum, tmp_nonJuncNum, tmp_keptJuncNum,
			tmp_changedJuncNum, tmp_changedJuncNum_hom, tmp_changedJuncNum_het,
			tmp_changedJuncNum_can2can, tmp_changedJuncNum_can2sem, tmp_changedJuncNum_can2non,
			tmp_changedJuncNum_sem2can, tmp_changedJuncNum_sem2sem, tmp_changedJuncNum_sem2non,
			tmp_changedJuncNum_non2can, tmp_changedJuncNum_non2sem, tmp_changedJuncNum_non2non,
			tmp_changedJuncNum_can2canCan, tmp_changedJuncNum_can2canSem, tmp_changedJuncNum_can2canNon,
			tmp_changedJuncNum_can2semSem, tmp_changedJuncNum_can2semNon, tmp_changedJuncNum_can2nonNon,
			tmp_changedJuncNum_sem2canCan, tmp_changedJuncNum_sem2canSem, tmp_changedJuncNum_sem2canNon,
			tmp_changedJuncNum_sem2semSem, tmp_changedJuncNum_sem2semNon, tmp_changedJuncNum_sem2nonNon,
			tmp_changedJuncNum_non2canCan, tmp_changedJuncNum_non2canSem, tmp_changedJuncNum_non2canNon,
			tmp_changedJuncNum_non2semSem, tmp_changedJuncNum_non2semNon, tmp_changedJuncNum_non2nonNon;
		int tmp_totalJuncNum_min3, tmp_canJuncNum_min3, tmp_semJuncNum_min3, tmp_nonJuncNum_min3, tmp_keptJuncNum_min3,
			tmp_changedJuncNum_min3, tmp_changedJuncNum_hom_min3, tmp_changedJuncNum_het_min3,
			tmp_changedJuncNum_can2can_min3, tmp_changedJuncNum_can2sem_min3, tmp_changedJuncNum_can2non_min3,
			tmp_changedJuncNum_sem2can_min3, tmp_changedJuncNum_sem2sem_min3, tmp_changedJuncNum_sem2non_min3,
			tmp_changedJuncNum_non2can_min3, tmp_changedJuncNum_non2sem_min3, tmp_changedJuncNum_non2non_min3,
			tmp_changedJuncNum_can2canCan_min3, tmp_changedJuncNum_can2canSem_min3, tmp_changedJuncNum_can2canNon_min3,
			tmp_changedJuncNum_can2semSem_min3, tmp_changedJuncNum_can2semNon_min3, tmp_changedJuncNum_can2nonNon_min3,
			tmp_changedJuncNum_sem2canCan_min3, tmp_changedJuncNum_sem2canSem_min3, tmp_changedJuncNum_sem2canNon_min3,
			tmp_changedJuncNum_sem2semSem_min3, tmp_changedJuncNum_sem2semNon_min3, tmp_changedJuncNum_sem2nonNon_min3,
			tmp_changedJuncNum_non2canCan_min3, tmp_changedJuncNum_non2canSem_min3, tmp_changedJuncNum_non2canNon_min3,
			tmp_changedJuncNum_non2semSem_min3, tmp_changedJuncNum_non2semNon_min3, tmp_changedJuncNum_non2nonNon_min3;
		int tmp_totalJuncNum_min5, tmp_canJuncNum_min5, tmp_semJuncNum_min5, tmp_nonJuncNum_min5, tmp_keptJuncNum_min5,
			tmp_changedJuncNum_min5, tmp_changedJuncNum_hom_min5, tmp_changedJuncNum_het_min5,
			tmp_changedJuncNum_can2can_min5, tmp_changedJuncNum_can2sem_min5, tmp_changedJuncNum_can2non_min5,
			tmp_changedJuncNum_sem2can_min5, tmp_changedJuncNum_sem2sem_min5, tmp_changedJuncNum_sem2non_min5,
			tmp_changedJuncNum_non2can_min5, tmp_changedJuncNum_non2sem_min5, tmp_changedJuncNum_non2non_min5,
			tmp_changedJuncNum_can2canCan_min5, tmp_changedJuncNum_can2canSem_min5, tmp_changedJuncNum_can2canNon_min5,
			tmp_changedJuncNum_can2semSem_min5, tmp_changedJuncNum_can2semNon_min5, tmp_changedJuncNum_can2nonNon_min5,
			tmp_changedJuncNum_sem2canCan_min5, tmp_changedJuncNum_sem2canSem_min5, tmp_changedJuncNum_sem2canNon_min5,
			tmp_changedJuncNum_sem2semSem_min5, tmp_changedJuncNum_sem2semNon_min5, tmp_changedJuncNum_sem2nonNon_min5,
			tmp_changedJuncNum_non2canCan_min5, tmp_changedJuncNum_non2canSem_min5, tmp_changedJuncNum_non2canNon_min5,
			tmp_changedJuncNum_non2semSem_min5, tmp_changedJuncNum_non2semNon_min5, tmp_changedJuncNum_non2nonNon_min5;					
		get3valueFromStr(tmpStr_totalJunc, tmp_totalJuncNum, tmp_totalJuncNum_min3, tmp_totalJuncNum_min5);
		get3valueFromStr(tmpStr_canJunc, tmp_canJuncNum, tmp_canJuncNum_min3, tmp_canJuncNum_min5);
		get3valueFromStr(tmpStr_semJunc, tmp_semJuncNum, tmp_semJuncNum_min3, tmp_semJuncNum_min5);
		get3valueFromStr(tmpStr_nonJunc, tmp_nonJuncNum, tmp_nonJuncNum_min3, tmp_nonJuncNum_min5);
		get3valueFromStr(tmpStr_keptJunc, tmp_keptJuncNum, tmp_keptJuncNum_min3, tmp_keptJuncNum_min5);
		get3valueFromStr(tmpStr_changedJunc, tmp_changedJuncNum, tmp_changedJuncNum_min3, tmp_changedJuncNum_min5);
		get3valueFromStr(tmpStr_changedJunc_hom, tmp_changedJuncNum_hom, tmp_changedJuncNum_hom_min3, tmp_changedJuncNum_hom_min5);
		get3valueFromStr(tmpStr_changedJunc_het, tmp_changedJuncNum_het, tmp_changedJuncNum_het_min3, tmp_changedJuncNum_het_min5);
		get3valueFromStr(tmpStr_changedJunc_can2can, tmp_changedJuncNum_can2can, tmp_changedJuncNum_can2can_min3, tmp_changedJuncNum_can2can_min5);
		get3valueFromStr(tmpStr_changedJunc_can2sem, tmp_changedJuncNum_can2sem, tmp_changedJuncNum_can2sem_min3, tmp_changedJuncNum_can2sem_min5);
		get3valueFromStr(tmpStr_changedJunc_can2non, tmp_changedJuncNum_can2non, tmp_changedJuncNum_can2non_min3, tmp_changedJuncNum_can2non_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2can, tmp_changedJuncNum_sem2can, tmp_changedJuncNum_sem2can_min3, tmp_changedJuncNum_sem2can_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2sem, tmp_changedJuncNum_sem2sem, tmp_changedJuncNum_sem2sem_min3, tmp_changedJuncNum_sem2sem_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2non, tmp_changedJuncNum_sem2non, tmp_changedJuncNum_sem2non_min3, tmp_changedJuncNum_sem2non_min5);
		get3valueFromStr(tmpStr_changedJunc_non2can, tmp_changedJuncNum_non2can, tmp_changedJuncNum_non2can_min3, tmp_changedJuncNum_non2can_min5);
		get3valueFromStr(tmpStr_changedJunc_non2sem, tmp_changedJuncNum_non2sem, tmp_changedJuncNum_non2sem_min3, tmp_changedJuncNum_non2sem_min5);
		get3valueFromStr(tmpStr_changedJunc_non2non, tmp_changedJuncNum_non2non, tmp_changedJuncNum_non2non_min3, tmp_changedJuncNum_non2non_min5);
		get3valueFromStr(tmpStr_changedJunc_can2canCan, tmp_changedJuncNum_can2canCan, tmp_changedJuncNum_can2canCan_min3, tmp_changedJuncNum_can2canCan_min5);
		get3valueFromStr(tmpStr_changedJunc_can2canSem, tmp_changedJuncNum_can2canSem, tmp_changedJuncNum_can2canSem_min3, tmp_changedJuncNum_can2canSem_min5);
		get3valueFromStr(tmpStr_changedJunc_can2canNon, tmp_changedJuncNum_can2canNon, tmp_changedJuncNum_can2canNon_min3, tmp_changedJuncNum_can2canNon_min5);
		get3valueFromStr(tmpStr_changedJunc_can2semSem, tmp_changedJuncNum_can2semSem, tmp_changedJuncNum_can2semSem_min3, tmp_changedJuncNum_can2semSem_min5);
		get3valueFromStr(tmpStr_changedJunc_can2semNon, tmp_changedJuncNum_can2semNon, tmp_changedJuncNum_can2semNon_min3, tmp_changedJuncNum_can2semNon_min5);
		get3valueFromStr(tmpStr_changedJunc_can2nonNon, tmp_changedJuncNum_can2nonNon, tmp_changedJuncNum_can2nonNon_min3, tmp_changedJuncNum_can2nonNon_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2canCan, tmp_changedJuncNum_sem2canCan, tmp_changedJuncNum_sem2canCan_min3, tmp_changedJuncNum_sem2canCan_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2canSem, tmp_changedJuncNum_sem2canSem, tmp_changedJuncNum_sem2canSem_min3, tmp_changedJuncNum_sem2canSem_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2canNon, tmp_changedJuncNum_sem2canNon, tmp_changedJuncNum_sem2canNon_min3, tmp_changedJuncNum_sem2canNon_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2semSem, tmp_changedJuncNum_sem2semSem, tmp_changedJuncNum_sem2semSem_min3, tmp_changedJuncNum_sem2semSem_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2semNon, tmp_changedJuncNum_sem2semNon, tmp_changedJuncNum_sem2semNon_min3, tmp_changedJuncNum_sem2semNon_min5);
		get3valueFromStr(tmpStr_changedJunc_sem2nonNon, tmp_changedJuncNum_sem2nonNon, tmp_changedJuncNum_sem2nonNon_min3, tmp_changedJuncNum_sem2nonNon_min5);
		get3valueFromStr(tmpStr_changedJunc_non2canCan, tmp_changedJuncNum_non2canCan, tmp_changedJuncNum_non2canCan_min3, tmp_changedJuncNum_non2canCan_min5);
		get3valueFromStr(tmpStr_changedJunc_non2canSem, tmp_changedJuncNum_non2canSem, tmp_changedJuncNum_non2canSem_min3, tmp_changedJuncNum_non2canSem_min5);
		get3valueFromStr(tmpStr_changedJunc_non2canNon, tmp_changedJuncNum_non2canNon, tmp_changedJuncNum_non2canNon_min3, tmp_changedJuncNum_non2canNon_min5);
		get3valueFromStr(tmpStr_changedJunc_non2semSem, tmp_changedJuncNum_non2semSem, tmp_changedJuncNum_non2semSem_min3, tmp_changedJuncNum_non2semSem_min5);
		get3valueFromStr(tmpStr_changedJunc_non2semNon, tmp_changedJuncNum_non2semNon, tmp_changedJuncNum_non2semNon_min3, tmp_changedJuncNum_non2semNon_min5);
		get3valueFromStr(tmpStr_changedJunc_non2nonNon, tmp_changedJuncNum_non2nonNon, tmp_changedJuncNum_non2nonNon_min3, tmp_changedJuncNum_non2nonNon_min5);

		totalJuncNumVec.push_back(tmp_totalJuncNum);
		canJuncNumVec.push_back(tmp_canJuncNum);
		semJuncNumVec.push_back(tmp_semJuncNum);
		nonJuncNumVec.push_back(tmp_nonJuncNum);
		keptJuncNumVec.push_back(tmp_keptJuncNum);
		changedJuncNumVec.push_back(tmp_changedJuncNum);
		changedJuncNumVec_hom.push_back(tmp_changedJuncNum_hom);
		changedJuncNumVec_het.push_back(tmp_changedJuncNum_het);
		changedJuncNumVec_can2can.push_back(tmp_changedJuncNum_can2can);
		changedJuncNumVec_can2sem.push_back(tmp_changedJuncNum_can2sem);
		changedJuncNumVec_can2non.push_back(tmp_changedJuncNum_can2non);
		changedJuncNumVec_sem2can.push_back(tmp_changedJuncNum_sem2can);
		changedJuncNumVec_sem2sem.push_back(tmp_changedJuncNum_sem2sem);
		changedJuncNumVec_sem2non.push_back(tmp_changedJuncNum_sem2non);
		changedJuncNumVec_non2can.push_back(tmp_changedJuncNum_non2can);
		changedJuncNumVec_non2sem.push_back(tmp_changedJuncNum_non2sem);
		changedJuncNumVec_non2non.push_back(tmp_changedJuncNum_non2non);
		changedJuncNumVec_can2canCan.push_back(tmp_changedJuncNum_can2canCan);
		changedJuncNumVec_can2canSem.push_back(tmp_changedJuncNum_can2canSem);
		changedJuncNumVec_can2canNon.push_back(tmp_changedJuncNum_can2canNon);
		changedJuncNumVec_can2semSem.push_back(tmp_changedJuncNum_can2semSem);
		changedJuncNumVec_can2semNon.push_back(tmp_changedJuncNum_can2semNon);
		changedJuncNumVec_can2nonNon.push_back(tmp_changedJuncNum_can2nonNon);
		changedJuncNumVec_sem2canCan.push_back(tmp_changedJuncNum_sem2canCan);
		changedJuncNumVec_sem2canSem.push_back(tmp_changedJuncNum_sem2canSem);
		changedJuncNumVec_sem2canNon.push_back(tmp_changedJuncNum_sem2canNon);
		changedJuncNumVec_sem2semSem.push_back(tmp_changedJuncNum_sem2semSem);
		changedJuncNumVec_sem2semNon.push_back(tmp_changedJuncNum_sem2semNon);
		changedJuncNumVec_sem2nonNon.push_back(tmp_changedJuncNum_sem2nonNon);
		changedJuncNumVec_non2canCan.push_back(tmp_changedJuncNum_non2canCan);
		changedJuncNumVec_non2canSem.push_back(tmp_changedJuncNum_non2canSem);
		changedJuncNumVec_non2canNon.push_back(tmp_changedJuncNum_non2canNon);
		changedJuncNumVec_non2semSem.push_back(tmp_changedJuncNum_non2semSem);
		changedJuncNumVec_non2semNon.push_back(tmp_changedJuncNum_non2semNon);
		changedJuncNumVec_non2nonNon.push_back(tmp_changedJuncNum_non2nonNon);

		totalJuncNumVec_min3.push_back(tmp_totalJuncNum_min3);
		canJuncNumVec_min3.push_back(tmp_canJuncNum_min3);
		semJuncNumVec_min3.push_back(tmp_semJuncNum_min3);
		nonJuncNumVec_min3.push_back(tmp_nonJuncNum_min3);
		keptJuncNumVec_min3.push_back(tmp_keptJuncNum_min3);
		changedJuncNumVec_min3.push_back(tmp_changedJuncNum_min3);
		changedJuncNumVec_hom_min3.push_back(tmp_changedJuncNum_hom_min3);
		changedJuncNumVec_het_min3.push_back(tmp_changedJuncNum_het_min3);
		changedJuncNumVec_can2can_min3.push_back(tmp_changedJuncNum_can2can_min3);
		changedJuncNumVec_can2sem_min3.push_back(tmp_changedJuncNum_can2sem_min3);
		changedJuncNumVec_can2non_min3.push_back(tmp_changedJuncNum_can2non_min3);
		changedJuncNumVec_sem2can_min3.push_back(tmp_changedJuncNum_sem2can_min3);
		changedJuncNumVec_sem2sem_min3.push_back(tmp_changedJuncNum_sem2sem_min3);
		changedJuncNumVec_sem2non_min3.push_back(tmp_changedJuncNum_sem2non_min3);
		changedJuncNumVec_non2can_min3.push_back(tmp_changedJuncNum_non2can_min3);
		changedJuncNumVec_non2sem_min3.push_back(tmp_changedJuncNum_non2sem_min3);
		changedJuncNumVec_non2non_min3.push_back(tmp_changedJuncNum_non2non_min3);
		changedJuncNumVec_can2canCan_min3.push_back(tmp_changedJuncNum_can2canCan_min3);
		changedJuncNumVec_can2canSem_min3.push_back(tmp_changedJuncNum_can2canSem_min3);
		changedJuncNumVec_can2canNon_min3.push_back(tmp_changedJuncNum_can2canNon_min3);
		changedJuncNumVec_can2semSem_min3.push_back(tmp_changedJuncNum_can2semSem_min3);
		changedJuncNumVec_can2semNon_min3.push_back(tmp_changedJuncNum_can2semNon_min3);
		changedJuncNumVec_can2nonNon_min3.push_back(tmp_changedJuncNum_can2nonNon_min3);
		changedJuncNumVec_sem2canCan_min3.push_back(tmp_changedJuncNum_sem2canCan_min3);
		changedJuncNumVec_sem2canSem_min3.push_back(tmp_changedJuncNum_sem2canSem_min3);
		changedJuncNumVec_sem2canNon_min3.push_back(tmp_changedJuncNum_sem2canNon_min3);
		changedJuncNumVec_sem2semSem_min3.push_back(tmp_changedJuncNum_sem2semSem_min3);
		changedJuncNumVec_sem2semNon_min3.push_back(tmp_changedJuncNum_sem2semNon_min3);
		changedJuncNumVec_sem2nonNon_min3.push_back(tmp_changedJuncNum_sem2nonNon_min3);
		changedJuncNumVec_non2canCan_min3.push_back(tmp_changedJuncNum_non2canCan_min3);
		changedJuncNumVec_non2canSem_min3.push_back(tmp_changedJuncNum_non2canSem_min3);
		changedJuncNumVec_non2canNon_min3.push_back(tmp_changedJuncNum_non2canNon_min3);
		changedJuncNumVec_non2semSem_min3.push_back(tmp_changedJuncNum_non2semSem_min3);
		changedJuncNumVec_non2semNon_min3.push_back(tmp_changedJuncNum_non2semNon_min3);
		changedJuncNumVec_non2nonNon_min3.push_back(tmp_changedJuncNum_non2nonNon_min3);

		totalJuncNumVec_min5.push_back(tmp_totalJuncNum_min5);
		canJuncNumVec_min5.push_back(tmp_canJuncNum_min5);
		semJuncNumVec_min5.push_back(tmp_semJuncNum_min5);
		nonJuncNumVec_min5.push_back(tmp_nonJuncNum_min5);
		keptJuncNumVec_min5.push_back(tmp_keptJuncNum_min5);
		changedJuncNumVec_min5.push_back(tmp_changedJuncNum_min5);
		changedJuncNumVec_hom_min5.push_back(tmp_changedJuncNum_hom_min5);
		changedJuncNumVec_het_min5.push_back(tmp_changedJuncNum_het_min5);
		changedJuncNumVec_can2can_min5.push_back(tmp_changedJuncNum_can2can_min5);
		changedJuncNumVec_can2sem_min5.push_back(tmp_changedJuncNum_can2sem_min5);
		changedJuncNumVec_can2non_min5.push_back(tmp_changedJuncNum_can2non_min5);
		changedJuncNumVec_sem2can_min5.push_back(tmp_changedJuncNum_sem2can_min5);
		changedJuncNumVec_sem2sem_min5.push_back(tmp_changedJuncNum_sem2sem_min5);
		changedJuncNumVec_sem2non_min5.push_back(tmp_changedJuncNum_sem2non_min5);
		changedJuncNumVec_non2can_min5.push_back(tmp_changedJuncNum_non2can_min5);
		changedJuncNumVec_non2sem_min5.push_back(tmp_changedJuncNum_non2sem_min5);
		changedJuncNumVec_non2non_min5.push_back(tmp_changedJuncNum_non2non_min5);
		changedJuncNumVec_can2canCan_min5.push_back(tmp_changedJuncNum_can2canCan_min5);
		changedJuncNumVec_can2canSem_min5.push_back(tmp_changedJuncNum_can2canSem_min5);
		changedJuncNumVec_can2canNon_min5.push_back(tmp_changedJuncNum_can2canNon_min5);
		changedJuncNumVec_can2semSem_min5.push_back(tmp_changedJuncNum_can2semSem_min5);
		changedJuncNumVec_can2semNon_min5.push_back(tmp_changedJuncNum_can2semNon_min5);
		changedJuncNumVec_can2nonNon_min5.push_back(tmp_changedJuncNum_can2nonNon_min5);
		changedJuncNumVec_sem2canCan_min5.push_back(tmp_changedJuncNum_sem2canCan_min5);
		changedJuncNumVec_sem2canSem_min5.push_back(tmp_changedJuncNum_sem2canSem_min5);
		changedJuncNumVec_sem2canNon_min5.push_back(tmp_changedJuncNum_sem2canNon_min5);
		changedJuncNumVec_sem2semSem_min5.push_back(tmp_changedJuncNum_sem2semSem_min5);
		changedJuncNumVec_sem2semNon_min5.push_back(tmp_changedJuncNum_sem2semNon_min5);
		changedJuncNumVec_sem2nonNon_min5.push_back(tmp_changedJuncNum_sem2nonNon_min5);
		changedJuncNumVec_non2canCan_min5.push_back(tmp_changedJuncNum_non2canCan_min5);
		changedJuncNumVec_non2canSem_min5.push_back(tmp_changedJuncNum_non2canSem_min5);
		changedJuncNumVec_non2canNon_min5.push_back(tmp_changedJuncNum_non2canNon_min5);
		changedJuncNumVec_non2semSem_min5.push_back(tmp_changedJuncNum_non2semSem_min5);
		changedJuncNumVec_non2semNon_min5.push_back(tmp_changedJuncNum_non2semNon_min5);
		changedJuncNumVec_non2nonNon_min5.push_back(tmp_changedJuncNum_non2nonNon_min5);		
	}

	int totalJuncNum_sum = 0, canJuncNum_sum = 0, semJuncNum_sum = 0, nonJuncNum_sum = 0, keptJuncNum_sum = 0,
		changedJuncNum_sum = 0, changedJuncNum_hom_sum = 0, changedJuncNum_het_sum = 0,
		changedJuncNum_can2can_sum = 0, changedJuncNum_can2sem_sum = 0, changedJuncNum_can2non_sum = 0,
		changedJuncNum_sem2can_sum = 0, changedJuncNum_sem2sem_sum = 0, changedJuncNum_sem2non_sum = 0,
		changedJuncNum_non2can_sum = 0, changedJuncNum_non2sem_sum = 0, changedJuncNum_non2non_sum = 0,
		changedJuncNum_can2canCan_sum = 0, changedJuncNum_can2canSem_sum = 0, changedJuncNum_can2canNon_sum = 0,
		changedJuncNum_can2semSem_sum = 0, changedJuncNum_can2semNon_sum = 0, changedJuncNum_can2nonNon_sum = 0,
		changedJuncNum_sem2canCan_sum = 0, changedJuncNum_sem2canSem_sum = 0, changedJuncNum_sem2canNon_sum = 0,
		changedJuncNum_sem2semSem_sum = 0, changedJuncNum_sem2semNon_sum = 0, changedJuncNum_sem2nonNon_sum = 0,
		changedJuncNum_non2canCan_sum = 0, changedJuncNum_non2canSem_sum = 0, changedJuncNum_non2canNon_sum = 0,
		changedJuncNum_non2semSem_sum = 0, changedJuncNum_non2semNon_sum = 0, changedJuncNum_non2nonNon_sum = 0;
	int totalJuncNum_min3_sum = 0, canJuncNum_min3_sum = 0, semJuncNum_min3_sum = 0, nonJuncNum_min3_sum = 0, keptJuncNum_min3_sum = 0,
		changedJuncNum_min3_sum = 0, changedJuncNum_hom_min3_sum = 0, changedJuncNum_het_min3_sum = 0,
		changedJuncNum_can2can_min3_sum = 0, changedJuncNum_can2sem_min3_sum = 0, changedJuncNum_can2non_min3_sum = 0,
		changedJuncNum_sem2can_min3_sum = 0, changedJuncNum_sem2sem_min3_sum = 0, changedJuncNum_sem2non_min3_sum = 0,
		changedJuncNum_non2can_min3_sum = 0, changedJuncNum_non2sem_min3_sum = 0, changedJuncNum_non2non_min3_sum = 0,
		changedJuncNum_can2canCan_min3_sum = 0, changedJuncNum_can2canSem_min3_sum = 0, changedJuncNum_can2canNon_min3_sum = 0,
		changedJuncNum_can2semSem_min3_sum = 0, changedJuncNum_can2semNon_min3_sum = 0, changedJuncNum_can2nonNon_min3_sum = 0,
		changedJuncNum_sem2canCan_min3_sum = 0, changedJuncNum_sem2canSem_min3_sum = 0, changedJuncNum_sem2canNon_min3_sum = 0,
		changedJuncNum_sem2semSem_min3_sum = 0, changedJuncNum_sem2semNon_min3_sum = 0, changedJuncNum_sem2nonNon_min3_sum = 0,
		changedJuncNum_non2canCan_min3_sum = 0, changedJuncNum_non2canSem_min3_sum = 0, changedJuncNum_non2canNon_min3_sum = 0,
		changedJuncNum_non2semSem_min3_sum = 0, changedJuncNum_non2semNon_min3_sum = 0, changedJuncNum_non2nonNon_min3_sum = 0;
	int totalJuncNum_min5_sum = 0, canJuncNum_min5_sum = 0, semJuncNum_min5_sum = 0, nonJuncNum_min5_sum = 0, keptJuncNum_min5_sum = 0,
		changedJuncNum_min5_sum = 0, changedJuncNum_hom_min5_sum = 0, changedJuncNum_het_min5_sum = 0,
		changedJuncNum_can2can_min5_sum = 0, changedJuncNum_can2sem_min5_sum = 0, changedJuncNum_can2non_min5_sum = 0,
		changedJuncNum_sem2can_min5_sum = 0, changedJuncNum_sem2sem_min5_sum = 0, changedJuncNum_sem2non_min5_sum = 0,
		changedJuncNum_non2can_min5_sum = 0, changedJuncNum_non2sem_min5_sum = 0, changedJuncNum_non2non_min5_sum = 0,
		changedJuncNum_can2canCan_min5_sum = 0, changedJuncNum_can2canSem_min5_sum = 0, changedJuncNum_can2canNon_min5_sum = 0,
		changedJuncNum_can2semSem_min5_sum = 0, changedJuncNum_can2semNon_min5_sum = 0, changedJuncNum_can2nonNon_min5_sum = 0,
		changedJuncNum_sem2canCan_min5_sum = 0, changedJuncNum_sem2canSem_min5_sum = 0, changedJuncNum_sem2canNon_min5_sum = 0,
		changedJuncNum_sem2semSem_min5_sum = 0, changedJuncNum_sem2semNon_min5_sum = 0, changedJuncNum_sem2nonNon_min5_sum = 0,
		changedJuncNum_non2canCan_min5_sum = 0, changedJuncNum_non2canSem_min5_sum = 0, changedJuncNum_non2canNon_min5_sum = 0,
		changedJuncNum_non2semSem_min5_sum = 0, changedJuncNum_non2semNon_min5_sum = 0, changedJuncNum_non2nonNon_min5_sum = 0;					
	
	for(int tmp = 0; tmp < psSJstatsFileVec.size(); tmp++)
	{
		totalJuncNum_sum += totalJuncNumVec[tmp];
		canJuncNum_sum += canJuncNumVec[tmp];
		semJuncNum_sum += semJuncNumVec[tmp];
		nonJuncNum_sum += nonJuncNumVec[tmp];
		keptJuncNum_sum += keptJuncNumVec[tmp];
		changedJuncNum_sum += changedJuncNumVec[tmp];
		changedJuncNum_hom_sum += changedJuncNumVec_hom[tmp];
		changedJuncNum_het_sum += changedJuncNumVec_het[tmp];
		changedJuncNum_can2can_sum += changedJuncNumVec_can2can[tmp];
		changedJuncNum_can2sem_sum += changedJuncNumVec_can2sem[tmp];
		changedJuncNum_can2non_sum += changedJuncNumVec_can2non[tmp];
		changedJuncNum_sem2can_sum += changedJuncNumVec_sem2can[tmp];
		changedJuncNum_sem2sem_sum += changedJuncNumVec_sem2sem[tmp];
		changedJuncNum_sem2non_sum += changedJuncNumVec_sem2non[tmp];
		changedJuncNum_non2can_sum += changedJuncNumVec_non2can[tmp];
		changedJuncNum_non2sem_sum += changedJuncNumVec_non2sem[tmp];
		changedJuncNum_non2non_sum += changedJuncNumVec_non2non[tmp];
		changedJuncNum_can2canCan_sum += changedJuncNumVec_can2canCan[tmp];
		changedJuncNum_can2canSem_sum += changedJuncNumVec_can2canSem[tmp];
		changedJuncNum_can2canNon_sum += changedJuncNumVec_can2canNon[tmp];
		changedJuncNum_can2semSem_sum += changedJuncNumVec_can2semSem[tmp];
		changedJuncNum_can2semNon_sum += changedJuncNumVec_can2semNon[tmp];
		changedJuncNum_can2nonNon_sum += changedJuncNumVec_can2nonNon[tmp];	
		changedJuncNum_sem2canCan_sum += changedJuncNumVec_sem2canCan[tmp];
		changedJuncNum_sem2canSem_sum += changedJuncNumVec_sem2canSem[tmp];
		changedJuncNum_sem2canNon_sum += changedJuncNumVec_sem2canNon[tmp];
		changedJuncNum_sem2semSem_sum += changedJuncNumVec_sem2semSem[tmp];
		changedJuncNum_sem2semNon_sum += changedJuncNumVec_sem2semNon[tmp];
		changedJuncNum_sem2nonNon_sum += changedJuncNumVec_sem2nonNon[tmp];
		changedJuncNum_non2canCan_sum += changedJuncNumVec_non2canCan[tmp];
		changedJuncNum_non2canSem_sum += changedJuncNumVec_non2canSem[tmp];
		changedJuncNum_non2canNon_sum += changedJuncNumVec_non2canNon[tmp];
		changedJuncNum_non2semSem_sum += changedJuncNumVec_non2semSem[tmp];
		changedJuncNum_non2semNon_sum += changedJuncNumVec_non2semNon[tmp];
		changedJuncNum_non2nonNon_sum += changedJuncNumVec_non2nonNon[tmp];

		totalJuncNum_min3_sum += totalJuncNumVec_min3[tmp];
		canJuncNum_min3_sum += canJuncNumVec_min3[tmp];
		semJuncNum_min3_sum += semJuncNumVec_min3[tmp];
		nonJuncNum_min3_sum += nonJuncNumVec_min3[tmp];
		keptJuncNum_min3_sum += keptJuncNumVec_min3[tmp];
		changedJuncNum_min3_sum += changedJuncNumVec_min3[tmp];
		changedJuncNum_hom_min3_sum += changedJuncNumVec_hom_min3[tmp];
		changedJuncNum_het_min3_sum += changedJuncNumVec_het_min3[tmp];
		changedJuncNum_can2can_min3_sum += changedJuncNumVec_can2can_min3[tmp];
		changedJuncNum_can2sem_min3_sum += changedJuncNumVec_can2sem_min3[tmp];
		changedJuncNum_can2non_min3_sum += changedJuncNumVec_can2non_min3[tmp];
		changedJuncNum_sem2can_min3_sum += changedJuncNumVec_sem2can_min3[tmp];
		changedJuncNum_sem2sem_min3_sum += changedJuncNumVec_sem2sem_min3[tmp];
		changedJuncNum_sem2non_min3_sum += changedJuncNumVec_sem2non_min3[tmp];
		changedJuncNum_non2can_min3_sum += changedJuncNumVec_non2can_min3[tmp];
		changedJuncNum_non2sem_min3_sum += changedJuncNumVec_non2sem_min3[tmp];
		changedJuncNum_non2non_min3_sum += changedJuncNumVec_non2non_min3[tmp];
		changedJuncNum_can2canCan_min3_sum += changedJuncNumVec_can2canCan_min3[tmp];
		changedJuncNum_can2canSem_min3_sum += changedJuncNumVec_can2canSem_min3[tmp];
		changedJuncNum_can2canNon_min3_sum += changedJuncNumVec_can2canNon_min3[tmp];
		changedJuncNum_can2semSem_min3_sum += changedJuncNumVec_can2semSem_min3[tmp];
		changedJuncNum_can2semNon_min3_sum += changedJuncNumVec_can2semNon_min3[tmp];
		changedJuncNum_can2nonNon_min3_sum += changedJuncNumVec_can2nonNon_min3[tmp];	
		changedJuncNum_sem2canCan_min3_sum += changedJuncNumVec_sem2canCan_min3[tmp];
		changedJuncNum_sem2canSem_min3_sum += changedJuncNumVec_sem2canSem_min3[tmp];
		changedJuncNum_sem2canNon_min3_sum += changedJuncNumVec_sem2canNon_min3[tmp];
		changedJuncNum_sem2semSem_min3_sum += changedJuncNumVec_sem2semSem_min3[tmp];
		changedJuncNum_sem2semNon_min3_sum += changedJuncNumVec_sem2semNon_min3[tmp];
		changedJuncNum_sem2nonNon_min3_sum += changedJuncNumVec_sem2nonNon_min3[tmp];
		changedJuncNum_non2canCan_min3_sum += changedJuncNumVec_non2canCan_min3[tmp];
		changedJuncNum_non2canSem_min3_sum += changedJuncNumVec_non2canSem_min3[tmp];
		changedJuncNum_non2canNon_min3_sum += changedJuncNumVec_non2canNon_min3[tmp];
		changedJuncNum_non2semSem_min3_sum += changedJuncNumVec_non2semSem_min3[tmp];
		changedJuncNum_non2semNon_min3_sum += changedJuncNumVec_non2semNon_min3[tmp];
		changedJuncNum_non2nonNon_min3_sum += changedJuncNumVec_non2nonNon_min3[tmp];

		totalJuncNum_min5_sum += totalJuncNumVec_min5[tmp];
		canJuncNum_min5_sum += canJuncNumVec_min5[tmp];
		semJuncNum_min5_sum += semJuncNumVec_min5[tmp];
		nonJuncNum_min5_sum += nonJuncNumVec_min5[tmp];
		keptJuncNum_min5_sum += keptJuncNumVec_min5[tmp];
		changedJuncNum_min5_sum += changedJuncNumVec_min5[tmp];
		changedJuncNum_hom_min5_sum += changedJuncNumVec_hom_min5[tmp];
		changedJuncNum_het_min5_sum += changedJuncNumVec_het_min5[tmp];
		changedJuncNum_can2can_min5_sum += changedJuncNumVec_can2can_min5[tmp];
		changedJuncNum_can2sem_min5_sum += changedJuncNumVec_can2sem_min5[tmp];
		changedJuncNum_can2non_min5_sum += changedJuncNumVec_can2non_min5[tmp];
		changedJuncNum_sem2can_min5_sum += changedJuncNumVec_sem2can_min5[tmp];
		changedJuncNum_sem2sem_min5_sum += changedJuncNumVec_sem2sem_min5[tmp];
		changedJuncNum_sem2non_min5_sum += changedJuncNumVec_sem2non_min5[tmp];
		changedJuncNum_non2can_min5_sum += changedJuncNumVec_non2can_min5[tmp];
		changedJuncNum_non2sem_min5_sum += changedJuncNumVec_non2sem_min5[tmp];
		changedJuncNum_non2non_min5_sum += changedJuncNumVec_non2non_min5[tmp];
		changedJuncNum_can2canCan_min5_sum += changedJuncNumVec_can2canCan_min5[tmp];
		changedJuncNum_can2canSem_min5_sum += changedJuncNumVec_can2canSem_min5[tmp];
		changedJuncNum_can2canNon_min5_sum += changedJuncNumVec_can2canNon_min5[tmp];
		changedJuncNum_can2semSem_min5_sum += changedJuncNumVec_can2semSem_min5[tmp];
		changedJuncNum_can2semNon_min5_sum += changedJuncNumVec_can2semNon_min5[tmp];
		changedJuncNum_can2nonNon_min5_sum += changedJuncNumVec_can2nonNon_min5[tmp];	
		changedJuncNum_sem2canCan_min5_sum += changedJuncNumVec_sem2canCan_min5[tmp];
		changedJuncNum_sem2canSem_min5_sum += changedJuncNumVec_sem2canSem_min5[tmp];
		changedJuncNum_sem2canNon_min5_sum += changedJuncNumVec_sem2canNon_min5[tmp];
		changedJuncNum_sem2semSem_min5_sum += changedJuncNumVec_sem2semSem_min5[tmp];
		changedJuncNum_sem2semNon_min5_sum += changedJuncNumVec_sem2semNon_min5[tmp];
		changedJuncNum_sem2nonNon_min5_sum += changedJuncNumVec_sem2nonNon_min5[tmp];
		changedJuncNum_non2canCan_min5_sum += changedJuncNumVec_non2canCan_min5[tmp];
		changedJuncNum_non2canSem_min5_sum += changedJuncNumVec_non2canSem_min5[tmp];
		changedJuncNum_non2canNon_min5_sum += changedJuncNumVec_non2canNon_min5[tmp];
		changedJuncNum_non2semSem_min5_sum += changedJuncNumVec_non2semSem_min5[tmp];
		changedJuncNum_non2semNon_min5_sum += changedJuncNumVec_non2semNon_min5[tmp];
		changedJuncNum_non2nonNon_min5_sum += changedJuncNumVec_non2nonNon_min5[tmp];	
	}

	ofstream summary_ofs(outputSummarizedResults.c_str());
	summary_ofs << "Junction Total #:\t" << totalJuncNum_sum << "\t" << totalJuncNum_min3_sum << "\t" << totalJuncNum_min5_sum << endl << endl;

	summary_ofs << "Junction Can #:\t" << canJuncNum_sum << "\t" << canJuncNum_min3_sum << "\t" << canJuncNum_min5_sum << endl;
	summary_ofs << "Junction Sem #:\t" << semJuncNum_sum << "\t" << semJuncNum_min3_sum << "\t" << semJuncNum_min5_sum << endl;
	summary_ofs << "Junction Non #:\t" << nonJuncNum_sum << "\t" << nonJuncNum_min3_sum << "\t" << nonJuncNum_min5_sum << endl << endl;

	summary_ofs << "Junction FlkStr Kept #:\t" << keptJuncNum_sum << "\t" << keptJuncNum_min3_sum << "\t" << keptJuncNum_min5_sum << endl << endl;
	summary_ofs << "Junction FlkStr Changed #:\t" << changedJuncNum_sum << "\t" << changedJuncNum_min3_sum << "\t" << changedJuncNum_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Hom #:\t" << changedJuncNum_hom_sum << "\t" << changedJuncNum_hom_min3_sum << "\t" << changedJuncNum_hom_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Het #:\t" << changedJuncNum_het_sum << "\t" << changedJuncNum_het_min3_sum << "\t" << changedJuncNum_het_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Can #:\t" << changedJuncNum_can2can_sum << "\t" << changedJuncNum_can2can_min3_sum << "\t" << changedJuncNum_can2can_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Sem #:\t" << changedJuncNum_can2sem_sum << "\t" << changedJuncNum_can2sem_min3_sum << "\t" << changedJuncNum_can2sem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Non #:\t" << changedJuncNum_can2non_sum << "\t" << changedJuncNum_can2non_min3_sum << "\t" << changedJuncNum_can2non_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Can #:\t" << changedJuncNum_sem2can_sum << "\t" << changedJuncNum_sem2can_min3_sum << "\t" << changedJuncNum_sem2can_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Sem #:\t" << changedJuncNum_sem2sem_sum << "\t" << changedJuncNum_sem2sem_min3_sum << "\t" << changedJuncNum_sem2sem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Non #:\t" << changedJuncNum_sem2non_sum << "\t" << changedJuncNum_sem2non_min3_sum << "\t" << changedJuncNum_sem2non_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Can #:\t" << changedJuncNum_non2can_sum << "\t" << changedJuncNum_non2can_min3_sum << "\t" << changedJuncNum_non2can_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Sem #:\t" << changedJuncNum_non2sem_sum << "\t" << changedJuncNum_non2sem_min3_sum << "\t" << changedJuncNum_non2sem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Non #:\t" << changedJuncNum_non2non_sum << "\t" << changedJuncNum_non2non_min3_sum << "\t" << changedJuncNum_non2non_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Can/Can #:\t" << changedJuncNum_can2canCan_sum << "\t" << changedJuncNum_can2canCan_min3_sum << "\t" << changedJuncNum_can2canCan_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Can/Sem #:\t" << changedJuncNum_can2canSem_sum << "\t" << changedJuncNum_can2canSem_min3_sum << "\t" << changedJuncNum_can2canSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Can/Non #:\t" << changedJuncNum_can2canNon_sum << "\t" << changedJuncNum_can2canNon_min3_sum << "\t" << changedJuncNum_can2canNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Sem/Sem #:\t" << changedJuncNum_can2semSem_sum << "\t" << changedJuncNum_can2semSem_min3_sum << "\t" << changedJuncNum_can2semSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Sem/Non #:\t" << changedJuncNum_can2semNon_sum << "\t" << changedJuncNum_can2semNon_min3_sum << "\t" << changedJuncNum_can2semNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Can to Non/Non #:\t" << changedJuncNum_can2nonNon_sum << "\t" << changedJuncNum_can2nonNon_min3_sum << "\t" << changedJuncNum_can2nonNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Can/Can #:\t" << changedJuncNum_sem2canCan_sum << "\t" << changedJuncNum_sem2canCan_min3_sum << "\t" << changedJuncNum_sem2canCan_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Can/Sem #:\t" << changedJuncNum_sem2canSem_sum << "\t" << changedJuncNum_sem2canSem_min3_sum << "\t" << changedJuncNum_sem2canSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Can/Non #:\t" << changedJuncNum_sem2canNon_sum << "\t" << changedJuncNum_sem2canNon_min3_sum << "\t" << changedJuncNum_sem2canNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Sem/Sem #:\t" << changedJuncNum_sem2semSem_sum << "\t" << changedJuncNum_sem2semSem_min3_sum << "\t" << changedJuncNum_sem2semSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Sem/Non #:\t" << changedJuncNum_sem2semNon_sum << "\t" << changedJuncNum_sem2semNon_min3_sum << "\t" << changedJuncNum_sem2semNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Sem to Non/Non #:\t" << changedJuncNum_sem2nonNon_sum << "\t" << changedJuncNum_sem2nonNon_min3_sum << "\t" << changedJuncNum_sem2nonNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Can/Can #:\t" << changedJuncNum_non2canCan_sum << "\t" << changedJuncNum_non2canCan_min3_sum << "\t" << changedJuncNum_non2canCan_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Can/Sem #:\t" << changedJuncNum_non2canSem_sum << "\t" << changedJuncNum_non2canSem_min3_sum << "\t" << changedJuncNum_non2canSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Can/Non #:\t" << changedJuncNum_non2canNon_sum << "\t" << changedJuncNum_non2canNon_min3_sum << "\t" << changedJuncNum_non2canNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Sem/Sem #:\t" << changedJuncNum_non2semSem_sum << "\t" << changedJuncNum_non2semSem_min3_sum << "\t" << changedJuncNum_non2semSem_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Sem/Non #:\t" << changedJuncNum_non2semNon_sum << "\t" << changedJuncNum_non2semNon_min3_sum << "\t" << changedJuncNum_non2semNon_min5_sum << endl;
	summary_ofs << "Junction FlkStr Changed Non to Non/Non #:\t" << changedJuncNum_non2nonNon_sum << "\t" << changedJuncNum_non2nonNon_min3_sum << "\t" << changedJuncNum_non2nonNon_min5_sum << endl;	
	summary_ofs.close();
	return 0;
}

