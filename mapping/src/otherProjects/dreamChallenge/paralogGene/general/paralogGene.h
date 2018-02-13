// paralogGeneStr: group_id NB_Genes Gene_Name

#ifndef PARALOGGENE_H
#define PARALOGGENE_H

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

class ParalogGeneGroup
{
private:
	int group_id;
	int gene_num;
	vector<string> geneNameVec;
public:
	ParalogGeneGroup()
	{}

	int return_groupId()
	{
		return group_id;
	}

	int return_geneNum()
	{
		return gene_num;
	}

	string return_geneName(int index)
	{
		return geneNameVec[index];
	}

	void initiate(int tmp_group_id, vector<string>& tmpGeneNameVec)
	{
		group_id = tmp_group_id;
		gene_num = tmpGeneNameVec.size();
		for(int tmp = 0; tmp < gene_num; tmp++)
			geneNameVec.push_back(tmpGeneNameVec[tmp]);
	}
};

class ParalogGeneGroupVec
{
private:
	vector<ParalogGeneGroup> paralogGeneGroupVec;
public:
	ParalogGeneGroupVec()
	{}

	int searchSingleGene_returnGroupIndex(string& tmpGeneName)
	{
		for(int tmp1 = 0; tmp1 < paralogGeneGroupVec.size(); tmp1++)
		{
			int tmpTargetGeneNum = paralogGeneGroupVec[tmp1].return_geneNum();
			for(int tmp2 = 0; tmp2 < tmpTargetGeneNum; tmp2++)
			{
				string tmpTargetGeneName = (paralogGeneGroupVec[tmp1]).return_geneName(tmp2);
				if(tmpGeneName == tmpTargetGeneName)
					return tmp1;
			}
		}
		return -1;
	}

	bool searchGenePairParalogOrNot(string& geneName_1, string& geneName_2)
	{
		int searchSingleGene_1_groupIndex = searchSingleGene_returnGroupIndex(geneName_1);
		int searchSingleGene_2_groupIndex = searchSingleGene_returnGroupIndex(geneName_2);
		if(searchSingleGene_1_groupIndex != searchSingleGene_2_groupIndex)
			return false;
		else
		{
			if(searchSingleGene_1_groupIndex == -1)
				return false;
			else
				return true;
		}
	}

	void initiate_duplicatedGeneDatabaseFile(string& tmpParalogGeneFile)
	{
		ifstream paralogGene_ifs(tmpParalogGeneFile.c_str());
		string tmp1stLine;
		getline(paralogGene_ifs, tmp1stLine);
		//int currentGroupId = -1;
		while(!paralogGene_ifs.eof())
		{
			string tmpStr;
			getline(paralogGene_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpGroupGeneNameVec;
			int tmpGroupId, tmpGroupGeneNum;

			string tmpGeneName_1st;
			parseParalogGeneStr(tmpStr, tmpGroupId, tmpGroupGeneNum, tmpGeneName_1st);
			tmpGroupGeneNameVec.push_back(tmpGeneName_1st);
			
			int tmpLeftGeneNum = tmpGroupGeneNum - 1;
			for(int tmp = 0; tmp < tmpLeftGeneNum; tmp++)
			{
				string tmpLeftParalogGeneStr;
				getline(paralogGene_ifs, tmpLeftParalogGeneStr);
				int tmpLeftGene_groupId, tmpLeftGene_groupGeneNum;
				string tmpLeftGene_geneName;
				parseParalogGeneStr(tmpLeftParalogGeneStr, tmpLeftGene_groupId, 
					tmpLeftGene_groupGeneNum, tmpLeftGene_geneName);
				tmpGroupGeneNameVec.push_back(tmpLeftGene_geneName);
			}
			ParalogGeneGroup tmpParalogGeneGroup;
			tmpParalogGeneGroup.initiate(tmpGroupId, tmpGroupGeneNameVec);
			paralogGeneGroupVec.push_back(tmpParalogGeneGroup);
		}
		paralogGene_ifs.close();
	}

	void initiate_reformattedParalogGeneGroupFile(string& tmpParalogGeneFile)
	{
		ifstream paralogGene_ifs(tmpParalogGeneFile.c_str());
		while(!paralogGene_ifs.eof())
		{		
			string tmpStr;
			getline(paralogGene_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			string tmpGroupIdStr = tmpStr.substr(0, tabLoc_1);
			int tmpGroupId = atoi(tmpGroupIdStr.c_str());
			string tmpGroupGeneNumStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			int tmpGroupGeneNum = atoi(tmpGroupGeneNumStr.c_str());
			string tmpGeneNameVecStr = tmpStr.substr(tabLoc_2 + 1);
			vector<string> tmpGroupGeneNameVec;
			int tmpStartLoc = 0;
			for(int tmp = 0; tmp < tmpGroupGeneNum; tmp++)
			{
				int tmpCommaLoc = tmpGeneNameVecStr.find(",", tmpStartLoc);
				string tmpGeneNameStr = tmpGeneNameVecStr.substr(tmpStartLoc, tmpCommaLoc - tmpStartLoc);
				tmpGroupGeneNameVec.push_back(tmpGeneNameStr);
				tmpStartLoc = tmpCommaLoc + 1;
			}
			ParalogGeneGroup tmpParalogGeneGroup;
			tmpParalogGeneGroup.initiate(tmpGroupId, tmpGroupGeneNameVec);
			paralogGeneGroupVec.push_back(tmpParalogGeneGroup);			
		}
		paralogGene_ifs.close();
	}

	void outputDuplicatedGeneDatabaseFile(string& tmpFile)
	{
		ofstream paralogGeneGroup_ofs(tmpFile.c_str());
		paralogGeneGroup_ofs << "group_id\tNB_Genes\tName" << endl;
		int groupNum = paralogGeneGroupVec.size();
		for(int tmp1 = 0; tmp1 < groupNum; tmp1++)
		{
			int groupId = paralogGeneGroupVec[tmp1].return_groupId();
			int groupGeneNum = paralogGeneGroupVec[tmp1].return_geneNum();
			for(int tmp2 = 0; tmp2 < groupGeneNum; tmp2++)
			{
				string tmpGeneName = paralogGeneGroupVec[tmp1].return_geneName(tmp2);
				paralogGeneGroup_ofs << groupId << "\t" << groupGeneNum << "\t" << tmpGeneName << endl;	
			}
		}
		paralogGeneGroup_ofs.close();
	}

	void outputReformattedParalogGeneGroupFile(string& tmpFile)
	{
		ofstream paralogGeneGroup_ofs(tmpFile.c_str());
		int groupNum = paralogGeneGroupVec.size();
		for(int tmp1 = 0; tmp1 < groupNum; tmp1++)
		{
			int groupId = paralogGeneGroupVec[tmp1].return_groupId();
			int groupGeneNum = paralogGeneGroupVec[tmp1].return_geneNum();
			string tmpGeneNameVecStr;
			for(int tmp2 = 0; tmp2 < groupGeneNum; tmp2++)
			{
				tmpGeneNameVecStr += paralogGeneGroupVec[tmp1].return_geneName(tmp2);
				tmpGeneNameVecStr += ",";
			}
			paralogGeneGroup_ofs << groupId << "\t" << groupGeneNum << "\t" << tmpGeneNameVecStr << endl;
		}
		paralogGeneGroup_ofs.close();
	}

	void parseParalogGeneStr(string& tmpStr, int& tmpGroupId, int& tmpGroupGeneNum, string& tmpGeneName)
	{
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpGroupIdStr = tmpStr.substr(0, tabLoc_1);
		string tmpGroupGeneNumStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		tmpGeneName = tmpStr.substr(tabLoc_2 + 1);
		tmpGroupId = atoi(tmpGroupIdStr.c_str());
		tmpGroupGeneNum = atoi(tmpGroupGeneNumStr.c_str());
	}
};

#endif