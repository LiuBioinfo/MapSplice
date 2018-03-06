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
#include "../general/paralogGene.h"
using namespace std;

void parseGeneInfoStr2geneNameVec(vector<string>& tmpGeneNameVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find(",", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpGeneNameVec.push_back(tmpField);
		//cout << "tmpField: " << tmpField << endl;
		startLoc = tabLoc + 1;
	}
	//tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

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
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 supNumMin" << endl;
		cout << "#2 paralogGeneFile" << endl; 
		cout << "#3 inputRawMPS3fusionResults" << endl;
		cout << "#4 outputFilteredFusionResultsPrefix" << endl;
		exit(1);
	}

	string supNumMinStr = argv[1];
	int supNumMin = atoi(supNumMinStr.c_str());
	//string paralogGeneFile = argv[2];
	string inputRawMPS3fusionResults = argv[3];
	string inputReformattedParalogGeneGroupFile = argv[2];
	ParalogGeneGroupVec tmpParalogGeneGroupVecInfo;
	tmpParalogGeneGroupVecInfo.initiate_reformattedParalogGeneGroupFile(inputReformattedParalogGeneGroupFile);

	string outputFilteredFusionResultsPrefix = argv[4];
	string outputFilteredFusionResultsPrefix_kept = outputFilteredFusionResultsPrefix + ".kept";
	string outputFilteredFusionResultsPrefix_filteredOut = outputFilteredFusionResultsPrefix + ".filteredOut";
	
	// filtering starts
	ifstream rawMPS3fusionResults_ifs(inputRawMPS3fusionResults.c_str());
	ofstream kept_ofs(outputFilteredFusionResultsPrefix_kept.c_str());
	ofstream filteredOut_ofs(outputFilteredFusionResultsPrefix_filteredOut.c_str());
	while(!rawMPS3fusionResults_ifs.eof())
	{
		string tmpStr;
		getline(rawMPS3fusionResults_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpSupNumStr = tmpFieldVec[9];
		string tmpGeneNameVec_1_infoStr = tmpFieldVec[tmpFieldVec.size()-2];
		string tmpGeneNameVec_2_infoStr = tmpFieldVec[tmpFieldVec.size()-1];
		vector<string> tmpGeneName_vec_1;
		parseGeneInfoStr2geneNameVec(tmpGeneName_vec_1, tmpGeneNameVec_1_infoStr);
		vector<string> tmpGeneName_vec_2;
		parseGeneInfoStr2geneNameVec(tmpGeneName_vec_2, tmpGeneNameVec_2_infoStr);
		bool tmpGenePairVec_paralogOrNot_or_similarGeneName_bool = false; 
		for(int tmp1 = 0; tmp1 < tmpGeneName_vec_1.size(); tmp1 ++)
		{
			string tmpGeneName_1 = tmpGeneName_vec_1[tmp1];
			for(int tmp2 = 0; tmp2 < tmpGeneName_vec_2.size(); tmp2 ++)
			{
				string tmpGeneName_2 = tmpGeneName_vec_2[tmp2];
				bool tmpGenePair_paralogOrNot_bool = tmpParalogGeneGroupVecInfo.searchGenePairParalogOrNot(tmpGeneName_1, tmpGeneName_2);
				bool tmpGenePair_similarGeneName_bool = twoGeneNameAtListSimilar(tmpGeneName_1, tmpGeneName_2);
				if(tmpGenePair_paralogOrNot_bool || tmpGenePair_similarGeneName_bool)
				{
					tmpGenePairVec_paralogOrNot_or_similarGeneName_bool = true;
					break;
				}
			}
			if(tmpGenePairVec_paralogOrNot_or_similarGeneName_bool)
				break;
		}
		int tmpSupNum = atoi(tmpSupNumStr.c_str());
		if((tmpSupNum >= supNumMin)&&(!tmpGenePairVec_paralogOrNot_or_similarGeneName_bool))
			kept_ofs << tmpStr << endl;
		else
			filteredOut_ofs << tmpStr << endl;
	}
	rawMPS3fusionResults_ifs.close();
	filteredOut_ofs.close();
	kept_ofs.close();

	//////////////////////////////
	// start to get genePair file 
	//////////////////////////////
	string kept_breakPoint_file = outputFilteredFusionResultsPrefix_kept;
	string kept_genePair_file = outputFilteredFusionResultsPrefix_kept + ".genePair";
	vector< pair< vector<string>, vector<string> > > fusionGeneVecPairVec;
	ifstream bp_ifs(kept_breakPoint_file.c_str());
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

	ofstream kept_fusionGenePair_ofs(kept_genePair_file.c_str());
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
		kept_fusionGenePair_ofs << tmpGene1VecStr << "\t" << tmpGene2VecStr << endl;
	}
	kept_fusionGenePair_ofs.close();
	return 0;
}