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

bool twoGeneNameAtListSimilar(string& tmpGene_1, string& tmpGene_2)
{
	if((tmpGene_1.length() < 2)||(tmpGene_2.length() < 2))
		return false;
	else if(tmpGene_1.substr(0,2) == tmpGene_2.substr(0,2))
		return true;
	else
		return false;
}

void parseOutGenePair(string& tmpBreakPointStr, vector<string>& tmpGeneVec_1, vector<string>& tmpGeneVec_2)
{
	int lastTabLoc = tmpBreakPointStr.rfind("\t");
	int penultimateTabLoc = tmpBreakPointStr.rfind("\t", lastTabLoc - 1);
	//cout << "lastTabLoc: " << lastTabLoc << endl;
	//cout << "penultimateTabLoc: " << penultimateTabLoc << endl;
	string tmpGeneVecStr_1 = tmpBreakPointStr.substr(penultimateTabLoc + 1, lastTabLoc - 1 - penultimateTabLoc + 1 + 1);
	string tmpGeneVecStr_2 = tmpBreakPointStr.substr(lastTabLoc + 1);
	int startCommaLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabCommaLoc = tmpGeneVecStr_1.find(",", startCommaLoc);
		if(tabCommaLoc == string::npos)
			break;
		else
		{
			string tmpField = tmpGeneVecStr_1.substr(startCommaLoc, tabCommaLoc - startCommaLoc);
			//cout << "tmpGeneId_1: " << tmpField << endl;
			tmpGeneVec_1.push_back(tmpField);
		}
		startCommaLoc = tabCommaLoc + 1;
	}
	startCommaLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabCommaLoc = tmpGeneVecStr_2.find(",", startCommaLoc);
		if(tabCommaLoc == string::npos)
			break;
		else
		{
			string tmpField = tmpGeneVecStr_2.substr(startCommaLoc, tabCommaLoc - startCommaLoc);
			//cout << "tmpGeneId_2: " << tmpField << endl;
			tmpGeneVec_2.push_back(tmpField);
		}
		startCommaLoc = tabCommaLoc + 1;
	}	
}

bool geneVecPair_theSame_bool(vector<string>& gene1Vec_pair1, vector<string>& gene2Vec_pair1,
	vector<string>& gene1Vec_pair2, vector<string>& gene2Vec_pair2)
{
	vector< pair<string, string> > genePairVec_1;
	vector< pair<string, string> > genePairVec_2;
	for(int tmp1 = 0; tmp1 < gene1Vec_pair1.size(); tmp1++)
	{
		string tmpGene1 = gene1Vec_pair1[tmp1];
		for(int tmp2 = 0; tmp2 < gene2Vec_pair1.size(); tmp2++)
		{
			string tmpGene2 = gene2Vec_pair1[tmp2];
			genePairVec_1.push_back(pair<string,string>(tmpGene1, tmpGene2));
		}
	}
	for(int tmp1 = 0; tmp1 < gene1Vec_pair2.size(); tmp1++)
	{
		string tmpGene1 = gene1Vec_pair2[tmp1];
		for(int tmp2 = 0; tmp2 < gene2Vec_pair2.size(); tmp2++)
		{
			string tmpGene2 = gene2Vec_pair2[tmp2];
			genePairVec_2.push_back(pair<string,string>(tmpGene1, tmpGene2));
		}
	}
	for(int tmp1 = 0; tmp1 < genePairVec_1.size(); tmp1++)
	{
		string tmpGene1_pair1 = genePairVec_1[tmp1].first;
		string tmpGene2_pair1 = genePairVec_1[tmp1].second;
		for(int tmp2 = 0; tmp2 < genePairVec_2.size(); tmp2++)
		{
			string tmpGene1_pair2 = genePairVec_2[tmp2].first;
			string tmpGene2_pair2 = genePairVec_2[tmp2].second;			
			if((twoGeneNameAtListSimilar(tmpGene1_pair1, tmpGene1_pair2)&&twoGeneNameAtListSimilar(tmpGene2_pair1,tmpGene2_pair2))
				||(twoGeneNameAtListSimilar(tmpGene1_pair1,tmpGene2_pair2)&&twoGeneNameAtListSimilar(tmpGene2_pair1,tmpGene1_pair2)))
				return true;
		}
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputAnnotatedFusionBreakPointFile outputFusionGenePairFile" << endl;
		exit(1);
	}
	string inputAnnotatedFusionBreakPointFile = argv[1];
	string outputFusionGenePairFile = argv[2];
	//cout << "start to load break point file" << endl;
	vector< pair< vector<string>, vector<string> > > fusionGeneVecPairVec;
	ifstream bp_ifs(inputAnnotatedFusionBreakPointFile.c_str());
	while(!bp_ifs.eof())
	{
		string tmpStr;
		getline(bp_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpGeneIdVec_1, tmpGeneIdVec_2;
		parseOutGenePair(tmpStr, tmpGeneIdVec_1, tmpGeneIdVec_2);
		bool already_exists_bool = false;
		int existing_geneVecPair_size = fusionGeneVecPairVec.size();
		for(int tmp = 0; tmp < existing_geneVecPair_size; tmp++)
		{
			bool tmp_theSameOrNot_bool = geneVecPair_theSame_bool(tmpGeneIdVec_1, tmpGeneIdVec_2,
				fusionGeneVecPairVec[tmp].first, fusionGeneVecPairVec[tmp].second);
			if(tmp_theSameOrNot_bool)
			{
				already_exists_bool = true;
				break;				
			}
		}
		if(!already_exists_bool)
			fusionGeneVecPairVec.push_back(pair< vector<string>, vector<string> >
				(tmpGeneIdVec_1, tmpGeneIdVec_2));
	}
	bp_ifs.close();
	//cout << "fusionGeneVecPairVec.size(): " << fusionGeneVecPairVec.size() << endl;

	ofstream fusionGenePair_ofs(outputFusionGenePairFile.c_str());
	for(int tmp1 = 0; tmp1 < fusionGeneVecPairVec.size(); tmp1++)
	{
		string tmpGene1VecStr = "";
		int tmp_gene1_vec_size = (fusionGeneVecPairVec[tmp1].first).size();
		for(int tmp2 = 0; tmp2 < tmp_gene1_vec_size; tmp2++)
		{	
			tmpGene1VecStr += (fusionGeneVecPairVec[tmp1].first)[tmp2];
			tmpGene1VecStr += ",";
		}
		string tmpGene2VecStr = "";
		int tmp_gene2_vec_size = (fusionGeneVecPairVec[tmp1].second).size();
		for(int tmp2 = 0; tmp2 < tmp_gene2_vec_size; tmp2++)
		{
			tmpGene2VecStr += (fusionGeneVecPairVec[tmp1].second)[tmp2];
			tmpGene2VecStr += ",";
		}
		fusionGenePair_ofs << tmpGene1VecStr << "\t" << tmpGene2VecStr << endl;
	}
	fusionGenePair_ofs.close();
	return 0;
}