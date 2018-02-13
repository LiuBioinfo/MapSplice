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

int searchCategoryIdInMap(string& tmpSRR, map<string, int>& srr2sizeCategoryMap)
{
	map<string, int>::iterator tmpIter = srr2sizeCategoryMap.find(tmpSRR);
	if(tmpIter != srr2sizeCategoryMap.end()) // found
		return (tmpIter->second);
	else
	{
		cout << "tmpSRR not found!" << endl;
		exit(1);
	}
}

int categoryStr2int(string& categoryStr)
{
	if(categoryStr == "0~300MB") // min_cov = 1
		return 1;
	else if(categoryStr == "300~500MB") // min_cov = 3
		return 2;
	else if(categoryStr == "500MB~1GB") // min_cov = 10
		return 3;
	else if(categoryStr == "1~3GB") // min_cov = 20
		return 4;
	else if(categoryStr == ">3GB") // min_cov = 50
		return 5;
	else
	{
		cout << "Error! invalid categoryStr: " << categoryStr << endl;
		exit(1);
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 SRR_list_sizeCategory" << endl;
		cout << "#2 input_diff_results" << endl;
		cout << "#3 results" << endl;
		exit(1);
	}	

	string SRR_list_sizeCategory_file = argv[1];
	string input_diff_results = argv[2];
	string results = argv[3];

	map<string, int> srr2sizeCategoryMap;
	ifstream srr_ifs(SRR_list_sizeCategory_file.c_str());
	while(!srr_ifs.eof())
	{
		string tmpStr;
		getline(srr_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpSRR = tmpStr.substr(0, tabLoc);
		string tmpSubStr = tmpStr.substr(tabLoc + 1);
		int tmpCategoryId = categoryStr2int(tmpSubStr);
		srr2sizeCategoryMap.insert(pair<string, int>(tmpSRR, tmpCategoryId));
	}
	srr_ifs.close();

	vector<int> category_numVec_0_1;
	category_numVec_0_1.push_back(0);
	category_numVec_0_1.push_back(0);
	category_numVec_0_1.push_back(0);
	category_numVec_0_1.push_back(0);
	category_numVec_0_1.push_back(0);
	vector<int> category_numVec_1_0;
	category_numVec_1_0.push_back(0);
	category_numVec_1_0.push_back(0);
	category_numVec_1_0.push_back(0);
	category_numVec_1_0.push_back(0);
	category_numVec_1_0.push_back(0);	

	int tmpQueryNum = 0;
	ifstream diff_ifs(input_diff_results.c_str());
	while(!diff_ifs.eof())
	{
		string tmpStr;
		getline(diff_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpQueryNum ++;
		if((tmpQueryNum / 1000000) * 1000000 == tmpQueryNum)
			cout << "tmpQueryNum: " << tmpQueryNum << endl;		
		if(tmpStr.substr(0,1) == ">")
			continue;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpSRR = tmpStr.substr(0, tabLoc_1);
		int tmpCategoryId = searchCategoryIdInMap(tmpSRR, srr2sizeCategoryMap);
		int tmpCategoryId_index = tmpCategoryId - 1;
		string tmpOccurrence_1st = tmpStr.substr(tabLoc_1 + 1, 1);
		string tmpOccurrence_2nd = tmpStr.substr(tabLoc_2 + 1, 1);
		if((tmpOccurrence_1st == "0")||(tmpOccurrence_2nd == "1"))
			category_numVec_0_1[tmpCategoryId_index] ++;
		else if((tmpOccurrence_1st == "1")||(tmpOccurrence_2nd == "0"))
			category_numVec_1_0[tmpCategoryId_index] ++;	
		else
		{
			cout << "incorrect tmpOccurrence_1st or tmpOccurrence_2nd!" << endl;
			cout << "tmpOccurrence_1st: " << tmpOccurrence_1st << endl;
			cout << "tmpOccurrence_2nd: " << tmpOccurrence_2nd << endl;
			exit(1);
		}
	}
	diff_ifs.close();

	ofstream results_ofs(results.c_str());
	results_ofs << "0-1:\t" << category_numVec_0_1[0]
		+ category_numVec_0_1[1] + category_numVec_0_1[2]
		+ category_numVec_0_1[3] + category_numVec_0_1[4] << endl;
	results_ofs << "\t0~300MB:\t" << category_numVec_0_1[0] << endl;
	results_ofs << "\t300~500MB:\t" << category_numVec_0_1[1] << endl;
	results_ofs << "\t500MB~1GB:\t" << category_numVec_0_1[2] << endl;
	results_ofs << "\t1~3GB:\t" << category_numVec_0_1[3] << endl;
	results_ofs << "\t>3GB:\t" << category_numVec_0_1[4] << endl;
	results_ofs << endl << "1-0:\t" << category_numVec_1_0[0]
		+ category_numVec_1_0[1] + category_numVec_1_0[2]
		+ category_numVec_1_0[3] + category_numVec_1_0[4] << endl;
	results_ofs << "\t0~300MB:\t" << category_numVec_1_0[0] << endl;
	results_ofs << "\t300~500MB:\t" << category_numVec_1_0[1] << endl;
	results_ofs << "\t500MB~1GB:\t" << category_numVec_1_0[2] << endl;
	results_ofs << "\t1~3GB:\t" << category_numVec_1_0[3] << endl;
	results_ofs << "\t>3GB:\t" << category_numVec_1_0[4] << endl;
	results_ofs.close();
	return 0;
}