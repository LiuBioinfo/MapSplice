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

void loadSRRlist(string& tmpSRRlist, vector<string>& tmpSRRvec)
{
	ifstream tmpSRR_ifs(tmpSRRlist.c_str());
	while(!tmpSRR_ifs.eof())
	{
		string tmpStr;
		getline(tmpSRR_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int SRRloc_start = tmpStr.rfind("SRR");
		int SRRloc_end = tmpStr.rfind(".");
		string tmpSRR;
		if(SRRloc_end != string::npos)
			tmpSRR = tmpStr.substr(SRRloc_start, SRRloc_start - SRRloc_end);
		else
			tmpSRR = tmpStr.substr(SRRloc_start);
		tmpSRRvec.push_back(tmpSRR);	
	}
	tmpSRR_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 SRRlist_1" << endl;
		cout << "#2 SRRlist_2" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string SRRlist_1 = argv[1];
	string SRRlist_2 = argv[2];
	string outputDir = argv[3];
	outputDir += "/";
	string cmd_mkdir = "mkdir -p " + outputDir;
	system(cmd_mkdir.c_str());

	vector<string> SRRvec_1;
	vector<string> SRRvec_2;
	loadSRRlist(SRRlist_1, SRRvec_1);
	loadSRRlist(SRRlist_2, SRRvec_2);

	string output_1_ori = outputDir + "ori_1.txt";
	string output_2_ori = outputDir + "ori_2.txt";
	ofstream ori_1_ofs(output_1_ori.c_str());
	ofstream ori_2_ofs(output_2_ori.c_str());
	for(int tmp = 0; tmp < SRRvec_1.size(); tmp++)
		ori_1_ofs << SRRvec_1[tmp] << endl;
	for(int tmp = 0; tmp < SRRvec_2.size(); tmp++)
		ori_2_ofs << SRRvec_2[tmp] << endl;	
	ori_1_ofs.close();
	ori_2_ofs.close();

	string output_bothShare = outputDir + "both_share.txt";
	string output_1_not_2 = outputDir + "1_2.txt";
	string output_2_not_1 = outputDir + "2_1.txt";
	vector<string> bothShareVec;
	vector<string> oneTwoVec;
	vector<string> twoOneVec;
	for(int tmp = 0; tmp < SRRvec_1.size(); tmp++)
	{
		string tmpSRR_1 = SRRvec_1[tmp];
		bool alsoInList2_bool = false;
		for(int tmp2 = 0; tmp2 < SRRvec_2.size(); tmp2++)
		{
			string tmpSRR_2 = SRRvec_2[tmp2];
			if(tmpSRR_1 == tmpSRR_2)
			{
				alsoInList2_bool = true;
				break;
			}
		}
		if(alsoInList2_bool)
			bothShareVec.push_back(tmpSRR_1);
		else
			oneTwoVec.push_back(tmpSRR_1);
	}
	for(int tmp = 0; tmp < SRRvec_2.size(); tmp++)
	{
		string tmpSRR_2 = SRRvec_2[tmp];
		bool alsoInList1_bool = false;
		for(int tmp1 = 0; tmp1 < SRRvec_1.size(); tmp1++)
		{
			string tmpSRR_1 = SRRvec_1[tmp1];
			if(tmpSRR_2 == tmpSRR_1)
			{
				alsoInList1_bool = true;
				break;
			}
		}
		if(!alsoInList1_bool)
			twoOneVec.push_back(tmpSRR_2);
	}

	ofstream bothShare_ofs(output_bothShare.c_str());
	ofstream one_two_ofs(output_1_not_2.c_str());
	ofstream two_one_ofs(output_2_not_1.c_str());
	for(int tmp = 0; tmp < bothShareVec.size(); tmp++)
		bothShare_ofs << bothShareVec[tmp] << endl;
	for(int tmp = 0; tmp < oneTwoVec.size(); tmp++)
		one_two_ofs << oneTwoVec[tmp] << endl;
	for(int tmp = 0; tmp < twoOneVec.size(); tmp++)
		two_one_ofs << twoOneVec[tmp] << endl;
	bothShare_ofs.close();
	one_two_ofs.close();
	two_one_ofs.close();
	return 0;
}