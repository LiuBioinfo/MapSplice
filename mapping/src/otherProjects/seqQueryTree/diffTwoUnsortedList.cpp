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

void loadList(string& listFile, vector<string>& tmpVec)
{
	ifstream list_ifs(listFile.c_str());
	while(!list_ifs.eof())
	{
		string tmpStr;
		getline(list_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpVec.push_back(tmpStr);
	}
	list_ifs.close();
}


int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 list_1" << endl;
		cout << "#2 list_2" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string list_1 = argv[1];
	string list_2 = argv[2];
	string outputDir = argv[3];
	outputDir += "/";
	string cmd_mkdir = "mkdir -p " + outputDir;
	system(cmd_mkdir.c_str());

	vector<string> vec_1;
	vector<string> vec_2;
	loadList(list_1, vec_1);
	loadList(list_2, vec_2);

	string output_1_ori = outputDir + "ori_1.txt";
	string output_2_ori = outputDir + "ori_2.txt";
	ofstream ori_1_ofs(output_1_ori.c_str());
	ofstream ori_2_ofs(output_2_ori.c_str());
	for(int tmp = 0; tmp < vec_1.size(); tmp++)
		ori_1_ofs << vec_1[tmp] << endl;
	for(int tmp = 0; tmp < vec_2.size(); tmp++)
		ori_2_ofs << vec_2[tmp] << endl;	
	ori_1_ofs.close();
	ori_2_ofs.close();

	string output_bothShare = outputDir + "both_share.txt";
	string output_1_not_2 = outputDir + "1_2.txt";
	string output_2_not_1 = outputDir + "2_1.txt";
	vector<string> bothShareVec;
	vector<string> oneTwoVec;
	vector<string> twoOneVec;
	for(int tmp = 0; tmp < vec_1.size(); tmp++)
	{
		string tmpId_1 = vec_1[tmp];
		bool alsoInList2_bool = false;
		for(int tmp2 = 0; tmp2 < vec_2.size(); tmp2++)
		{
			string tmpId_2 = vec_2[tmp2];
			if(tmpId_1 == tmpId_2)
			{
				alsoInList2_bool = true;
				break;
			}
		}
		if(alsoInList2_bool)
			bothShareVec.push_back(tmpId_1);
		else
			oneTwoVec.push_back(tmpId_1);
	}
	for(int tmp = 0; tmp < vec_2.size(); tmp++)
	{
		string tmpId_2 = vec_2[tmp];
		bool alsoInList1_bool = false;
		for(int tmp1 = 0; tmp1 < vec_1.size(); tmp1++)
		{
			string tmpId_1 = vec_1[tmp1];
			if(tmpId_2 == tmpId_1)
			{
				alsoInList1_bool = true;
				break;
			}
		}
		if(!alsoInList1_bool)
			twoOneVec.push_back(tmpId_2);
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