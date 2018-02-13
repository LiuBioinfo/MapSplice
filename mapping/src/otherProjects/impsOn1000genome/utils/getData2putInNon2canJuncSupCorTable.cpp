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

using namespace std;

void get_non2can_sampleNum_supNum(
	vector<int>& tmpCanoHapNumVec, vector<int>& tmpReadSupNumVec, 
	int& tmp_non2can_sampleNum, int& tmp_non2can_supNum)
{
	tmp_non2can_sampleNum = 0;
	tmp_non2can_supNum = 0; 
	for(int tmp = 0; tmp < tmpCanoHapNumVec.size(); tmp++)
	{
		int tmpCanoHapNum = tmpCanoHapNumVec[tmp];
		int tmpReadSupNum = tmpReadSupNumVec[tmp];
		if(((tmpCanoHapNum == 1)||(tmpCanoHapNum == 2))&&(tmpReadSupNum > 0))
		{
			tmp_non2can_sampleNum ++;
			tmp_non2can_supNum += tmpReadSupNum;
		}
	}
}

int return_junc_eligible(vector< pair<int,int> >& non2can_sampleSupNumPairVec, 
	int sampleNum_thres, int supNum_thres)
{
	int tmpEligibleJuncNum = 0;
	for(int tmp = 0; tmp < non2can_sampleSupNumPairVec.size(); tmp++)
	{
		int tmpSampleNum = non2can_sampleSupNumPairVec[tmp].first;
		int tmpSupNum = non2can_sampleSupNumPairVec[tmp].second;
		if((tmpSampleNum >= sampleNum_thres)&&(tmpSupNum >= supNum_thres))
			tmpEligibleJuncNum ++;
	}
	return tmpEligibleJuncNum;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSJexpInEachSampleTable" << endl;
		cout << "#2 outputJuncCountSupCorrTable" << endl;
		exit(1);
	}
	string inputSJexpInEachSampleTable = argv[1];
	string outputSJreadCountSampleSupCorrTable = argv[2];
	ifstream expInEachSample_ifs(inputSJexpInEachSampleTable.c_str());
	string tmpHeader;
	getline(expInEachSample_ifs, tmpHeader);
	vector< pair<int,int> > non2can_sampleSupNumPairVec;
	while(!expInEachSample_ifs.eof())
	{
		string tmpStr;
		getline(expInEachSample_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> fieldVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{
				fieldVec.push_back(tmpStr.substr(startLoc));
				break;
			}
			string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
			fieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		string tmpChrName = fieldVec[0];
		string tmpStartPosStr = fieldVec[1];
		string tmpEndPosStr = fieldVec[2];
		string tmpRefCanSJhapNumStr = fieldVec[3];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		int tmpRefCanSJhapNum = atoi(tmpRefCanSJhapNumStr.c_str());
		int sampleNum = (fieldVec.size() - 4)/2;
		vector<int> tmpCanoHapNumVec;
		vector<int> tmpReadSupNumVec;
		for(int tmp = 0; tmp < sampleNum; tmp++)
		{
			string tmpSample_hapNumStr = fieldVec[tmp*2 + 4];
			string tmpSample_supNumStr = fieldVec[tmp*2 + 5];
			int tmpSample_hapNum = atoi(tmpSample_hapNumStr.c_str());
			int tmpSample_supNum = atoi(tmpSample_supNumStr.c_str());
			tmpCanoHapNumVec.push_back(tmpSample_hapNum);
			tmpReadSupNumVec.push_back(tmpSample_supNum);
		}
		if(tmpRefCanSJhapNum == 0)
		{	
			int tmp_non2can_sampleNum = 0, tmp_non2can_supNum = 0;
			get_non2can_sampleNum_supNum(tmpCanoHapNumVec, tmpReadSupNumVec, 
				tmp_non2can_sampleNum, tmp_non2can_supNum);
			if((tmp_non2can_sampleNum > 0)&&(tmp_non2can_supNum > 0))
				non2can_sampleSupNumPairVec.push_back(pair<int,int>(
					tmp_non2can_sampleNum, tmp_non2can_supNum));
		}
	}
	expInEachSample_ifs.close();

	ofstream readCountSampleSupCorr_ofs(outputSJreadCountSampleSupCorrTable.c_str());
	int junc_num_indiv_1_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 1, 2);
	int junc_num_indiv_1_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 1, 5);
	int junc_num_indiv_1_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 1, 10);
	int junc_num_indiv_1_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 1, 20);	
	int junc_num_indiv_1_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 1, 50);

	int junc_num_indiv_2_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 2, 2);
	int junc_num_indiv_2_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 2, 5);
	int junc_num_indiv_2_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 2, 10);
	int junc_num_indiv_2_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 2, 20);	
	int junc_num_indiv_2_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 2, 50);

	int junc_num_indiv_5_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 5, 2);
	int junc_num_indiv_5_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 5, 5);
	int junc_num_indiv_5_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 5, 10);
	int junc_num_indiv_5_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 5, 20);	
	int junc_num_indiv_5_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 5, 50);

	int junc_num_indiv_10_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 10, 2);
	int junc_num_indiv_10_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 10, 5);
	int junc_num_indiv_10_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 10, 10);
	int junc_num_indiv_10_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 10, 20);	
	int junc_num_indiv_10_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 10, 50);

	int junc_num_indiv_20_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 20, 2);
	int junc_num_indiv_20_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 20, 5);
	int junc_num_indiv_20_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 20, 10);
	int junc_num_indiv_20_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 20, 20);	
	int junc_num_indiv_20_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 20, 50);

	int junc_num_indiv_50_sup_2 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 50, 2);
	int junc_num_indiv_50_sup_5 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 50, 5);
	int junc_num_indiv_50_sup_10 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 50, 10);
	int junc_num_indiv_50_sup_20 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 50, 20);	
	int junc_num_indiv_50_sup_50 = return_junc_eligible(
		non2can_sampleSupNumPairVec, 50, 50);												

	readCountSampleSupCorr_ofs << "ind#\\sup#\t2\t5\t10\t20\t50" << endl;
	readCountSampleSupCorr_ofs << "1\t" << junc_num_indiv_1_sup_2 << "\t"
		<< junc_num_indiv_1_sup_5 << "\t" << junc_num_indiv_1_sup_10 << "\t"
		<< junc_num_indiv_1_sup_20 << "\t" << junc_num_indiv_1_sup_50 << endl;
	
	readCountSampleSupCorr_ofs << "2\t" << junc_num_indiv_2_sup_2 << "\t"
		<< junc_num_indiv_2_sup_5 << "\t" << junc_num_indiv_2_sup_10 << "\t"
		<< junc_num_indiv_2_sup_20 << "\t" << junc_num_indiv_2_sup_50 << endl;
	
	readCountSampleSupCorr_ofs << "5\t" << junc_num_indiv_5_sup_2 << "\t"
		<< junc_num_indiv_5_sup_5 << "\t" << junc_num_indiv_5_sup_10 << "\t"
		<< junc_num_indiv_5_sup_20 << "\t" << junc_num_indiv_5_sup_50 << endl;		
	
	readCountSampleSupCorr_ofs << "10\t" << junc_num_indiv_10_sup_2 << "\t"
		<< junc_num_indiv_10_sup_5 << "\t" << junc_num_indiv_10_sup_10 << "\t"
		<< junc_num_indiv_10_sup_20 << "\t" << junc_num_indiv_10_sup_50 << endl;

	readCountSampleSupCorr_ofs << "20\t" << junc_num_indiv_20_sup_2 << "\t"
		<< junc_num_indiv_20_sup_5 << "\t" << junc_num_indiv_20_sup_10 << "\t"
		<< junc_num_indiv_20_sup_20 << "\t" << junc_num_indiv_20_sup_50 << endl;		
	
	readCountSampleSupCorr_ofs << "50\t" << junc_num_indiv_50_sup_2 << "\t"
		<< junc_num_indiv_50_sup_5 << "\t" << junc_num_indiv_50_sup_10 << "\t"
		<< junc_num_indiv_50_sup_20 << "\t" << junc_num_indiv_50_sup_50 << endl;

	readCountSampleSupCorr_ofs.close();
	return 0;
}