// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALU_ENTRY_HASH_H
#define ALU_ENTRY_HASH_H

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

#include "alu_entry.h"

using namespace std;

typedef map<int, set<int> > AluEntryArea2infoIndexMap;

class Alu_Entry_Hash
{
private:
	vector<AluEntryArea2infoIndexMap> aluEntryArea2infoIndexMapVec;
	int areaSize;
public:
	vector<Alu_Entry> aluAnnEntryInfoVec;

	Alu_Entry_Hash()
	{
		areaSize = 1000;
	}

	void initiate_aluAnnEntryArea2infoIndexMapVec(int tmpChromTotalNum)
	{
		for(int tmp = 0; tmp < tmpChromTotalNum; tmp ++)
		{
			AluEntryArea2infoIndexMap tmpAluEntryArea2infoIndexMap;
			aluEntryArea2infoIndexMapVec.push_back(tmpAluEntryArea2infoIndexMap);
		}		
	}

	void loadAluAnn(string& tmpSimplifiedAluAnnFile, Index_Info* indexInfo, 
		string& alu_valid_file, string& alu_invalid_file)
	{
		ofstream alu_valid_ofs(alu_valid_file.c_str());
		ofstream alu_invalid_ofs(alu_invalid_file.c_str());
		ifstream tmpSimplifiedAluAnn_ifs(tmpSimplifiedAluAnnFile.c_str());
		while(!tmpSimplifiedAluAnn_ifs.eof())
		{
			string tmpStr;
			getline(tmpSimplifiedAluAnn_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			string tmpAlu_chrName = tmpStr.substr(0, tabLoc_1);
			int tmpAlu_chrNameInt = indexInfo->convertStringToInt(tmpAlu_chrName);
			if(tmpAlu_chrNameInt < 0)
				alu_invalid_ofs << tmpStr << endl;
			else
			{	
				alu_valid_ofs << tmpStr << endl;
				this->insert_newAluAnnEntry(tmpStr, indexInfo);
			}
		}
		tmpSimplifiedAluAnn_ifs.close();
		alu_invalid_ofs.close();
		alu_valid_ofs.close();
	}

	void insert_newAluAnnEntry(string& tmpNewAluAnnEntry, Index_Info* indexInfo)
	{
		//string tmpChrName, tmpStrand, tmpGeneId, tmpTranscriptId,
		//int tmpStartPos, tmpEndPos;
		vector<string> tmpFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 8; tmp++)
		{
			int tabLoc = tmpNewAluAnnEntry.find("\t", startLoc);
			string tmpField = tmpNewAluAnnEntry.substr(startLoc, tabLoc-startLoc);
			tmpFieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		tmpFieldVec.push_back(tmpNewAluAnnEntry.substr(startLoc));
		string tmpChrName = tmpFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			//cout << "invalid chr name, exiting ......: " << tmpChrName << endl;
			return;
		}
		
		string tmpStartPosStr = tmpFieldVec[3];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		string tmpEndPosStr = tmpFieldVec[4];
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpStrand = tmpFieldVec[6];
		string tmpOtherStr = tmpFieldVec[8];
		int quoteLoc_1 = tmpOtherStr.find("\"");
		int quoteLoc_2 = tmpOtherStr.find("\"", quoteLoc_1 + 1);
		int quoteLoc_3 = tmpOtherStr.find("\"", quoteLoc_2 + 1);
		int quoteLoc_4 = tmpOtherStr.find("\"", quoteLoc_3 + 1);
		string tmpGeneId = tmpOtherStr.substr(quoteLoc_1 + 1, quoteLoc_2 - quoteLoc_1 - 1);
		string tmpTranscriptId = tmpOtherStr.substr(quoteLoc_3 + 1, quoteLoc_4 - quoteLoc_3 - 1);
		//string tmpGeneId = tmpFieldVec[8];
		//string tmpTranscriptId = tmpFieldVec[9];
		//cout << "tmpGeneId: " << tmpGeneId << endl;
		//cout << "tmpTranscriptId: " << tmpTranscriptId << endl;
		//cout << "inserting new Entry: " << endl << tmpNewGeneAnnEntry << endl;
		this->insert_newAlu(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpGeneId, tmpTranscriptId);
	}

	void insert_newAlu(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, 
		string& tmpStrand, string& tmpGeneId, string& tmpTranscriptId)
	{
		Alu_Entry tmpAluAnnEntryInfo;
		tmpAluAnnEntryInfo.initiate(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpGeneId, tmpTranscriptId);
		aluAnnEntryInfoVec.push_back(tmpAluAnnEntryInfo);
		int tmpAluAnnEntryInfoIndex = this->returnAluAnnEntryNum() - 1;
		//cout << "current geneAnnEntryNum: " << this->returnGeneAnnEntryNum() << endl;
		int areaNO_startPos = (int)(tmpStartPos/areaSize);
		int areaNO_endPos = (int)(tmpEndPos/areaSize);
		//cout << "areaNO_startPos: " << areaNO_startPos << endl;
		//cout << "areaNO_endPos: " << areaNO_endPos << endl;
		for(int tmpAreaNO = areaNO_startPos; tmpAreaNO <= areaNO_endPos; tmpAreaNO ++)
		{
			AluEntryArea2infoIndexMap::iterator tmpAreaIter;
			tmpAreaIter	= aluEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
			if(tmpAreaIter != aluEntryArea2infoIndexMapVec[tmpChrNameInt].end())
				(tmpAreaIter->second).insert(tmpAluAnnEntryInfoIndex);
			else
			{
				//cout << "inserting area code " << endl;
				set<int> tmpSet;
				tmpSet.insert(tmpAluAnnEntryInfoIndex);
				aluEntryArea2infoIndexMapVec[tmpChrNameInt].insert(pair<int, set<int> >(tmpAreaNO, tmpSet));
			}	
		}
	}

	void searchAndReturnAluAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		string& tmpChrName, int tmpPos, Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
			return;
		this->searchAndReturnAluAnnEntryStrVec(tmpAnnEntryInfoStrVec, tmpChrNameInt, tmpPos, indexInfo);
	}

	void searchAndReturnAluAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		if(tmpChrNameInt < 0)
			return;
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnAluAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo);
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpEntryInfoStr = this->returnAluAnnEntry_str(tmpIndex, indexInfo);
			tmpAnnEntryInfoStrVec.push_back(tmpEntryInfoStr);
		}
	}

	void searchAndReturnAluTranscriptIdVec(vector<string>& tmpAnnTranscriptIdVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		if(tmpChrNameInt < 0)
			return;
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnAluAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo);
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpAluTranscriptId = this->returnAluAnnEntry_transcriptId(tmpIndex);
			tmpAnnTranscriptIdVec.push_back(tmpAluTranscriptId);
		}
	}

	void searchAndReturnAluGeneIdVec(vector<string>& tmpAnnGeneIdVec, 
		string& tmpChrName, int tmpPos, Index_Info* indexInfo)	
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		this->searchAndReturnAluGeneIdVec(tmpAnnGeneIdVec, tmpChrNameInt, tmpPos, indexInfo);
	}

	void searchAndReturnAluGeneIdVec(vector<string>& tmpAnnGeneIdVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		if(tmpChrNameInt < 0)
			return;
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnAluAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo);
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpAluGeneId = this->returnAluAnnEntry_geneId(tmpIndex);
			tmpAnnGeneIdVec.push_back(tmpAluGeneId);
		}
	}

	void searchAndReturnAluAnnEntryInfoIndexVec(vector<int>& tmpAnnEntryInfoIndexVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		if(tmpChrNameInt < 0)
			return;
		int tmpAreaNO = (int)(tmpPos/areaSize);
		//cout << "tmpAreaNO: " << tmpAreaNO << endl;
		AluEntryArea2infoIndexMap::iterator tmpAreaIter 
			= aluEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
		if(tmpAreaIter == aluEntryArea2infoIndexMapVec[tmpChrNameInt].end())
			return;
		else
		{
			//cout << "area found !" << endl;
			for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); 
				tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
			{
				int tmpEntryInfoIndex = (*tmpSetIter);
				//cout << "tmpEntryInfoIndex: " << tmpEntryInfoIndex << endl;
				int tmpEntry_startPos = this->returnAluAnnEntry_startPos(tmpEntryInfoIndex);
				int tmpEntry_endPos = this->returnAluAnnEntry_endPos(tmpEntryInfoIndex);
				//cout << "tmpEntry_startPos: " << tmpEntry_startPos << endl;
				//cout << "tmpEntry_endPos: " << tmpEntry_endPos << endl;
				if((tmpEntry_startPos <= tmpPos)&&(tmpPos <= tmpEntry_endPos))
					tmpAnnEntryInfoIndexVec.push_back(tmpEntryInfoIndex);
			}
		}		
	}

	string returnAluAnnEntry_str(int index, Index_Info* indexInfo)
	{
		return aluAnnEntryInfoVec[index].returnAnnEntry(indexInfo);
	}

	int returnAluAnnEntry_chrNameInt(int index)
	{
		return aluAnnEntryInfoVec[index].returnChrNameInt();
	}

	int returnAluAnnEntry_startPos(int index)
	{
		return aluAnnEntryInfoVec[index].returnStartPos();
	}

	int returnAluAnnEntry_endPos(int index)
	{
		return aluAnnEntryInfoVec[index].returnEndPos();
	}

	string returnAluAnnEntry_geneId(int index)
	{
		return aluAnnEntryInfoVec[index].returnGeneId();
	}

	string returnAluAnnEntry_transcriptId(int index)
	{
		return aluAnnEntryInfoVec[index].returnTranscriptId();
	}

	int returnAluAnnEntryNum()
	{
		return aluAnnEntryInfoVec.size();
	}	
};

#endif