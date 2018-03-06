/*
fusion list:
A,B,	C,
D,	E,F
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <vector>
#include <map>
#include <set>

using namespace std;

vector<string> parse_str_tab(string& tmpStr)
{
	vector<string> tmpFieldVec;
	int start_loc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmp_tab_loc = tmpStr.find("\t", start_loc);
		if(tmp_tab_loc == string::npos)
		{
			tmpFieldVec.push_back(tmpStr.substr(start_loc));
			break;
		}
		tmpFieldVec.push_back(tmpStr.substr(start_loc, tmp_tab_loc - start_loc));
		start_loc = tmp_tab_loc + 1;
		if(start_loc >= tmpStr.length())
			break;
	}
	return tmpFieldVec;
}

vector<string> parse_str_comma(string& tmpStr)
{
	vector<string> tmpFieldVec;
	int start_loc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmp_tab_loc = tmpStr.find(",", start_loc);
		if(tmp_tab_loc == string::npos)
		{
			tmpFieldVec.push_back(tmpStr.substr(start_loc));
			break;
		}
		tmpFieldVec.push_back(tmpStr.substr(start_loc, tmp_tab_loc - start_loc));
		start_loc = tmp_tab_loc + 1;
		if(start_loc >= tmpStr.length())
			break;
	}
	return tmpFieldVec;
}

vector< pair< vector<string>, vector<string> > > loadFusionList(string inputFile)
{
	vector< pair< vector<string>, vector<string> > > fusionVec;
	ifstream fus_ifs(inputFile.c_str());
	while(!fus_ifs.eof())
	{
		string tmpStr;
		getline(fus_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tabVec = parse_str_tab(tmpStr);
		vector<string> geneVec_1 = parse_str_comma(tabVec[0]);
		vector<string> geneVec_2 = parse_str_comma(tabVec[1]);
		fusionVec.push_back(pair< vector<string>, vector<string> >(geneVec_1, geneVec_2));
	}
	fus_ifs.close();
	return fusionVec;
}

bool atLeastOneInGeneListMatch(vector<string> geneList_1, vector<string> geneList_2)
{
	for(int tmp1 = 0; tmp1 < geneList_1.size(); tmp1++)
	{
		string tmpGene_1 = geneList_1[tmp1];
		for(int tmp2 = 0; tmp2 < geneList_2.size(); tmp2++)
		{
			string tmpGene_2 = geneList_2[tmp2];
			if(tmpGene_1 == tmpGene_2)
				return true;
		}
	}
	return false;
}

bool fusionPairMatchBool(pair< vector<string>, vector<string> > fusionResPair, 
	pair< vector<string>, vector<string> > fusionGtPair)
{
	if((atLeastOneInGeneListMatch(fusionResPair.first, fusionGtPair.first)
		&&atLeastOneInGeneListMatch(fusionResPair.second, fusionGtPair.second))
		||(atLeastOneInGeneListMatch(fusionResPair.first, fusionGtPair.second)
		&&atLeastOneInGeneListMatch(fusionResPair.second, fusionGtPair.first)))
		return true;
	else
		return false;
}

string print_fusionListPair(pair< vector<string>, vector<string> > fusionPair)
{
	string geneList_1 = "", geneList_2 = "";
	for(int tmp = 0; tmp < (fusionPair.first).size(); tmp++)
	{
		geneList_1 += (fusionPair.first)[tmp];
		geneList_1 += ",";
	}
	for(int tmp = 0; tmp < (fusionPair.second).size(); tmp++)
	{
		geneList_2 += (fusionPair.second)[tmp];
		geneList_2 += ",";
	}
	return geneList_1 + "\t" + geneList_2;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputFusionRes" << endl;
		cout << "#2 inputFusionGt" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string inputFusionRes = argv[1];
	string inputFusionGt = argv[2];
	string outputDir = argv[3];
	outputDir += "/";
	string cmd_mkdir =  "mkdir " + outputDir;
	system(cmd_mkdir.c_str());

	vector< pair< vector<string>, vector<string> > > fusionResVec = loadFusionList(inputFusionRes);
	vector< pair< vector<string>, vector<string> > > fusionGtVec = loadFusionList(inputFusionGt);
	int fusionResNum = fusionResVec.size();
	vector<bool> hitVec_res(fusionResNum, false);
	int fusionGtNum = fusionGtVec.size();
	vector<bool> hitVec_gt(fusionGtNum, false);

	for(int tmpResIndex = 0; tmpResIndex < fusionResNum; tmpResIndex++)
	{
		for(int tmpGtIndex = 0; tmpGtIndex < fusionGtNum; tmpGtIndex++)
		{
			if(fusionPairMatchBool(fusionResVec[tmpResIndex], fusionGtVec[tmpGtIndex]))
			{
				hitVec_res[tmpResIndex] = true;
				hitVec_gt[tmpGtIndex] = true;
			}
		}
	}

	string output_fusionRes = outputDir + "fusion.res.cmp";
	string output_fusionGt = outputDir + "fusion.gt.cmp";
	string output_fusionCmp = outputDir + "cmp.stats";
	
	int res_true_num = 0, res_false_num = 0,
		gt_hit_num = 0, gt_miss_num = 0;

	ofstream res_ofs(output_fusionRes.c_str());
	for(int tmp = 0; tmp < fusionResVec.size(); tmp++)
	{	
		string tmpFusionListPair = print_fusionListPair(fusionResVec[tmp]);
		if(hitVec_res[tmp])
		{
			res_ofs << tmpFusionListPair << "\tTrue" << endl;
			res_true_num++; 
		}
		else
		{
			res_ofs << tmpFusionListPair << "\tFalse" << endl;
			res_false_num++;			
		}
	}
	res_ofs.close();

	ofstream gt_ofs(output_fusionGt.c_str());
	for(int tmp = 0; tmp < fusionGtVec.size(); tmp++)
	{	
		string tmpFusionListPair = print_fusionListPair(fusionGtVec[tmp]);
		if(hitVec_gt[tmp])
		{
			gt_ofs << tmpFusionListPair << "\tHit" << endl;
			gt_hit_num++; 
		}
		else
		{
			gt_ofs << tmpFusionListPair << "\tMiss" << endl;
			gt_miss_num++;			
		}
	}
	gt_ofs.close();

	ofstream cmp_ofs(output_fusionCmp.c_str());
	cmp_ofs << "Detect_True_#:\t" << res_true_num << endl
		<< "Detect_False_#:\t" << res_false_num << endl
		<< "GroundTruth_Hit_#:\t" << gt_hit_num << endl
		<< "GroundTruth_Miss_#:\t" << gt_miss_num << endl;
	cmp_ofs.close();
	return 0;
}