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
#include "../../../../general/read_block_test.h"
#include "../../../../general/otherFunc.h"
#include "../../../../general/index_info.h"

using namespace std;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("_", startLoc);
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
	if(argc != 5)
	{
		cout << "Executable inputRegionAssignedReadFa_1 inputRegionAssignedReadFa_2 totalRegionNum outputReadFaPrefix" << endl;;
		exit(1);
	}
	string inputRegionAssignedReadFa_1 = argv[1];
	string inputRegionAssignedReadFa_2 = argv[2];
	string totalRegionNumStr = argv[3];
	int totalRegionNum = atoi(totalRegionNumStr.c_str());
	string outputReadFaPrefix = argv[4];
	cout << "totalRegionNum: " << totalRegionNum << endl;
	cout << "outputReadFaPrefix: " << outputReadFaPrefix << endl;
	//vector<string> outputReadFaFileVec;
	vector<ofstream*> readFaOfsVec_1;
	vector<ofstream*> readFaOfsVec_2;
	for(int tmp = 1; tmp <= totalRegionNum; tmp++)
	{
		string tmpRegionFaFile_1 = outputReadFaPrefix + "." + int_to_str(tmp) + ".1.fa";
		string tmpRegionFaFile_2 = outputReadFaPrefix + "." + int_to_str(tmp) + ".2.fa";
		cout << "tmpRegionFaFile_1: " << tmpRegionFaFile_1 << endl;
		cout << "tmpRegionFaFile_2: " << tmpRegionFaFile_2 << endl;
		//outputReadFaFileVec.push_back(tmpRegionFaFile);
		ofstream* tmpOfs_1 = new ofstream(tmpRegionFaFile_1.c_str());
		ofstream* tmpOfs_2 = new ofstream(tmpRegionFaFile_2.c_str());
		readFaOfsVec_1.push_back(tmpOfs_1);
		readFaOfsVec_2.push_back(tmpOfs_2);
	}

	ifstream rawFa_ifs_1(inputRegionAssignedReadFa_1.c_str());
	ifstream rawFa_ifs_2(inputRegionAssignedReadFa_2.c_str());
	while((!rawFa_ifs_1.eof())||(!rawFa_ifs_2.eof()))
	{
		string tmpStr_1_id;
		getline(rawFa_ifs_1, tmpStr_1_id);
		string tmpStr_2_id;
		getline(rawFa_ifs_2, tmpStr_2_id);
		if((tmpStr_1_id == "")||(tmpStr_2_id == ""))
			break;
		string tmpStr_1_seq;
		getline(rawFa_ifs_1, tmpStr_1_seq);
		string tmpStr_2_seq;
		getline(rawFa_ifs_2, tmpStr_2_seq);
		vector<string> idFieldVec_1;
		//vector<string> idFieldVec_2;
		parseStr2fieldVec(idFieldVec_1, tmpStr_1_id);
		//parseStr2fieldVec(idFieldVec_2, tmpStr_2_seq);
		string tmpRegionIdStr = idFieldVec_1[4];
		int tmpRegionId = atoi(tmpRegionIdStr.c_str());
		int tmpRegionIndex = tmpRegionId - 1;
		*(readFaOfsVec_1[tmpRegionIndex]) << tmpStr_1_id << endl << tmpStr_1_seq << endl;
		*(readFaOfsVec_2[tmpRegionIndex]) << tmpStr_2_id << endl << tmpStr_2_seq << endl;
	}

	for(int tmp = 1; tmp <= totalRegionNum; tmp++)
	{
		(*readFaOfsVec_1[tmp-1]).close();
		(*readFaOfsVec_2[tmp-1]).close();
		delete readFaOfsVec_1[tmp-1];
		delete readFaOfsVec_2[tmp-1];
	}
	rawFa_ifs_1.close();
	rawFa_ifs_2.close();
	return 0;
}