// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef GENEANNENTRYHASH_INFO_H
#define GENEANNENTRYHASH_INFO_H

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

#include "geneAnnEntry.h"

using namespace std;

//typedef map<int, int> GeneAnnEntryInfoIndexMap;
typedef map<int, set<int> > GeneAnnEntryArea2infoIndexMap; //( areaNO = pos/1000 ) intermediate hash to directly get all possible SNPs near a certain position
//typedef SNPareaHash::iterator SNPareaMapIter

class GeneAnnEntry_Hash_Info
{
private:
	vector<GeneAnnEntryArea2infoIndexMap> geneAnnEntryArea2infoIndexMapVec;
	//vector<GeneAnnEntryAreaMap> geneAnnEntryAreaMapVec;
	int areaSize;
public:
	vector<GeneAnnEntry_Info> geneAnnEntryInfoVec;

	GeneAnnEntry_Hash_Info()
	{
		areaSize = 1000;
	}

	string searchAndReturnGeneName_withGeneId(string& tmpGeneId)
	{
		for(int tmp = 0; tmp < geneAnnEntryInfoVec.size(); tmp++)
		{
			string tmpGeneId_inGeneAnnEntry = geneAnnEntryInfoVec[tmp].returnGeneId();
			if(tmpGeneId_inGeneAnnEntry == tmpGeneId)
				return geneAnnEntryInfoVec[tmp].returnGeneName();
		}
		cout << "invalid tmpGeneId ..." << endl;
		exit(1);
	}

	void initiateGeneAnnEntryInfoVec(string& tmpSimplifiedGeneAnnFile, Index_Info* indexInfo)
	{
		ifstream tmpSimplifiedGeneAnn_ifs(tmpSimplifiedGeneAnnFile.c_str());
		while(!tmpSimplifiedGeneAnn_ifs.eof())
		{
			string tmpStr;
			getline(tmpSimplifiedGeneAnn_ifs, tmpStr);
			if(tmpStr == "")
				break;
			this->addNewGeneAnnEntry2infoVec(tmpStr, indexInfo);
		}
		tmpSimplifiedGeneAnn_ifs.close();
	}

	void initiateGeneAnnEntryArea2infoIndexMapVec(string& tmpSimplifiedGeneAnnIndexFile, Index_Info* indexInfo)
	{
		int tmpChromTotalNum = indexInfo->returnChromNum();
		for(int tmp = 0; tmp < tmpChromTotalNum; tmp ++)
		{
			GeneAnnEntryArea2infoIndexMap tmpGeneAnnEntryArea2infoIndexMap;
			geneAnnEntryArea2infoIndexMapVec.push_back(tmpGeneAnnEntryArea2infoIndexMap);
		}	

		ifstream tmpSimplifiedGeneAnnIndex_ifs(tmpSimplifiedGeneAnnIndexFile.c_str());
		while(!tmpSimplifiedGeneAnnIndex_ifs.eof())
		{
			string tmpStr;
			getline(tmpSimplifiedGeneAnnIndex_ifs, tmpStr);
			if(tmpStr == "")
				break;
			this->addNewGeneAnnEntry2indexMap(tmpStr, indexInfo);
		}
		tmpSimplifiedGeneAnnIndex_ifs.close();		
	}	

	void addNewGeneAnnEntry2indexMap(string& tmpStr, Index_Info* indexInfo)
	{
		int tabLoc_1 = tmpStr.find(":");
		int tabLoc_2 = tmpStr.find(":", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find(":", tabLoc_2 + 1);
		string tmpChrName = tmpStr.substr(0, tabLoc_1);
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		string tmpAreaNOstr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tmpAreaNO = atoi(tmpAreaNOstr.c_str());
		string tmpTotalEntryNumStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpTotalEntryNum = atoi(tmpTotalEntryNumStr.c_str());
		string tmpEntryStr = tmpStr.substr(tabLoc_3 + 1);
		set<int> tmpSet;
		int startLoc = 0;
		for(int tmpEntry = 0; tmpEntry < tmpTotalEntryNum; tmpEntry ++)
		{
			int tmpSemiCommaLoc = tmpEntryStr.find(",", startLoc);
			string tmpEntryIndexStr = tmpEntryStr.substr(startLoc, tmpSemiCommaLoc - startLoc);
			int tmpEntryIndexInt = atoi(tmpEntryIndexStr.c_str());
			tmpSet.insert(tmpEntryIndexInt);
			startLoc = tmpSemiCommaLoc + 1;
		}
		geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].insert(pair<int, set<int> >(tmpAreaNO, tmpSet));
	}

	void outputIndex(string& indexFile, Index_Info* indexInfo) // each line: chrNameInt:tmpAreaNO:tmpAreaEntryNum:tmpAreaEntryStr
	{
		ofstream index_ofs(indexFile.c_str());
		for(int tmp = 0; tmp < geneAnnEntryArea2infoIndexMapVec.size(); tmp++)
		{
			string tmpChrName = indexInfo->returnChrNameStr(tmp);
			GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter;
			for(tmpAreaIter = geneAnnEntryArea2infoIndexMapVec[tmp].begin();
				tmpAreaIter != geneAnnEntryArea2infoIndexMapVec[tmp].end(); tmpAreaIter ++)
			{
				int tmpAreaNO = tmpAreaIter->first;
				index_ofs << tmpChrName << ":" << tmpAreaNO << ":";
				int tmpAreaEntryNum = 0;
				string tmpAreaEntryStr = "";
				for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin();
					tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
				{
					tmpAreaEntryNum ++;
					int tmpEntryIndex = (*tmpSetIter);
					tmpAreaEntryStr += int_to_str(tmpEntryIndex);
					tmpAreaEntryStr += ",";
				}
				index_ofs << tmpAreaEntryNum << ":" << tmpAreaEntryStr << endl;
			}
		}
		index_ofs.close();
	}

	void returnLeftBoundaryPosVec(vector<int>& leftBoundaryPosVec, int tmpChrNameInt, int tmpStartPosInChr, int tmpEndPosInChr) // originally developed for fusion detection	
	{
		vector<int> leftBoundaryPosVec_raw;
		int areaNO_startPos = (int)(tmpStartPosInChr/areaSize);
		int areaNO_endPos = (int)(tmpEndPosInChr/areaSize);
		for(int tmpAreaNO = areaNO_startPos; tmpAreaNO <= areaNO_endPos; tmpAreaNO ++)
		{
			GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter 
				= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
			if(tmpAreaIter == geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
				continue;
			else
			{
				for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
				{
					int tmpEntryInfoIndex = (*tmpSetIter);
					int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpEntryInfoIndex);
					//int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpEntryInfoIndex);
					if((tmpEntry_startPos >= tmpStartPosInChr)&&(tmpEntry_startPos <= tmpEndPosInChr))
						leftBoundaryPosVec_raw.push_back(tmpEntry_startPos);
				}
			}			
		}
		for(int tmp1 = 0; tmp1 < leftBoundaryPosVec_raw.size(); tmp1 ++)
		{
			int tmp1_pos = leftBoundaryPosVec_raw[tmp1];
			int currentBoundaryPosVecSize = leftBoundaryPosVec.size();
			bool alreadyExists_bool = false;
			for(int tmp2 = 0; tmp2 < currentBoundaryPosVecSize; tmp2++)
			{
				int tmp2_pos = leftBoundaryPosVec[tmp2];
				if(tmp1_pos == tmp2_pos)
				{
					alreadyExists_bool = true;
					break;
				}
			}
			if(!alreadyExists_bool)
				leftBoundaryPosVec.push_back(tmp1_pos);
		}
	}

	void returnRightBoundaryPosVec(vector<int>& rightBoundaryPosVec, int tmpChrNameInt, int tmpStartPosInChr, int tmpEndPosInChr) // originally developed for fusion detection	
	{
		vector<int> rightBoundaryPosVec_raw;
		int areaNO_startPos = (int)(tmpStartPosInChr/areaSize);
		int areaNO_endPos = (int)(tmpEndPosInChr/areaSize);
		for(int tmpAreaNO = areaNO_startPos; tmpAreaNO <= areaNO_endPos; tmpAreaNO ++)
		{
			GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter 
				= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
			if(tmpAreaIter == geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
				continue;
			else
			{
				for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
				{
					int tmpEntryInfoIndex = (*tmpSetIter);
					//int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpEntryInfoIndex);
					int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpEntryInfoIndex);
					if((tmpEntry_endPos >= tmpStartPosInChr)&&(tmpEntry_endPos <= tmpEndPosInChr))
						rightBoundaryPosVec_raw.push_back(tmpEntry_endPos);
				}
			}			
		}
		for(int tmp1 = 0; tmp1 < rightBoundaryPosVec_raw.size(); tmp1 ++)
		{
			int tmp1_pos = rightBoundaryPosVec_raw[tmp1];
			int currentBoundaryPosVecSize = rightBoundaryPosVec.size();
			bool alreadyExists_bool = false;
			for(int tmp2 = 0; tmp2 < currentBoundaryPosVecSize; tmp2++)
			{
				int tmp2_pos = rightBoundaryPosVec[tmp2];
				if(tmp1_pos == tmp2_pos)
				{
					alreadyExists_bool = true;
					break;
				}
			}
			if(!alreadyExists_bool)
				rightBoundaryPosVec.push_back(tmp1_pos);
		}
	}

	void searchAndReturnGeneAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		string& tmpChrNameStr, int tmpPos, Index_Info* indexInfo)
	{
		this->searchAndReturnGeneAnnEntryStrVec(tmpAnnEntryInfoStrVec, tmpChrNameStr, tmpPos, indexInfo, 0);
	}

	void searchAndReturnGeneAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		string& tmpChrNameStr, int tmpPos, Index_Info* indexInfo, int mode)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		this->searchAndReturnGeneAnnEntryStrVec(tmpAnnEntryInfoStrVec, tmpChrNameInt, tmpPos, indexInfo, mode);
	}

	void searchAndReturnGeneAnnEntryGeneIdVec(vector<string>& tmpEntryGeneIdVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		this->searchAndReturnGeneAnnEntryGeneIdVec(tmpEntryGeneIdVec, tmpChrNameInt, tmpPos, indexInfo, 0);
	}	

	void searchAndReturnGeneAnnEntryGeneIdVec(vector<string>& tmpEntryGeneIdVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo, int mode)	
	{
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnGeneAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo, mode);
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpEntryGeneIdStr = this->returnGeneAnnEntry_geneId(tmpIndex);
			tmpEntryGeneIdVec.push_back(tmpEntryGeneIdStr);
		}		
	}

	string searchAndReturnGeneAnnEntryStrand_within(string& tmpChrName, int tmpChrPos, Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		return this->searchAndReturnGeneAnnEntryStrand_within(tmpChrNameInt, tmpChrPos, indexInfo);
	}

	string searchAndReturnGeneAnnEntryStrand_within(int tmpChrNameInt, int tmpChrPos, Index_Info* indexInfo)
	{
		return this->searchAndReturnGeneAnnEntryStrand(tmpChrNameInt, tmpChrPos, indexInfo, 0);
	}

	string searchAndReturnGeneAnnEntryStrand_bondaryOnly(string& tmpChrName, int tmpChrPos, Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		return this->searchAndReturnGeneAnnEntryStrand_bondaryOnly(tmpChrNameInt, tmpChrPos, indexInfo);
	}

	string searchAndReturnGeneAnnEntryStrand_bondaryOnly(int tmpChrNameInt, int tmpChrPos, Index_Info* indexInfo)
	{
		return this->searchAndReturnGeneAnnEntryStrand(tmpChrNameInt, tmpChrPos, indexInfo, 1);
	}

	string searchAndReturnGeneAnnEntryStrand(int tmpChrNameInt, int tmpChrPos, Index_Info* indexInfo, int mode)
	{
		vector<int> tmpAnnEntryInfoIndexVec;
		this->searchAndReturnGeneAnnEntryInfoIndexVec(tmpAnnEntryInfoIndexVec, tmpChrNameInt, tmpChrPos, indexInfo, mode);
		if(tmpAnnEntryInfoIndexVec.size() == 0)
			return "O";
		else
		{
			bool sense_strand_geneAnnEntry_exists_bool = false;
			bool antisense_strand_geneAnnEntry_exists_bool = false;
			for(int tmp = 0; tmp < tmpAnnEntryInfoIndexVec.size(); tmp++)
			{
				int tmpIndex = tmpAnnEntryInfoIndexVec[tmp];
				string tmpStrand = this->returnGeneAnnEntry_strand(tmpIndex);
				if(tmpStrand == "+")
					sense_strand_geneAnnEntry_exists_bool = true;
				else if(tmpStrand == "-")
					antisense_strand_geneAnnEntry_exists_bool = true;
				else
				{
					cout << "error ! invalid strand: " << tmpStrand << endl;
					exit(1);
				}
			}
			if(sense_strand_geneAnnEntry_exists_bool && (!antisense_strand_geneAnnEntry_exists_bool))
				return "+";
			else if((!sense_strand_geneAnnEntry_exists_bool) && antisense_strand_geneAnnEntry_exists_bool)
				return "-";
			else if(sense_strand_geneAnnEntry_exists_bool && antisense_strand_geneAnnEntry_exists_bool)
				return "X";
		}
	}

	void searchAndReturnGeneAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		this->searchAndReturnGeneAnnEntryStrVec(tmpAnnEntryInfoStrVec, tmpChrNameInt, tmpPos, indexInfo, 0);
	}

	void searchAndReturnGeneAnnEntryStrVec(vector<string>& tmpAnnEntryInfoStrVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo, int mode)
	{
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnGeneAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo, mode);
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpEntryInfoStr = this->returnGeneAnnEntry_str(tmpIndex, indexInfo);
			tmpAnnEntryInfoStrVec.push_back(tmpEntryInfoStr);
		}
	}

	void searchAndReturnGeneAnnEntryInfoIndexVec(vector<int>& tmpAnnEntryInfoIndexVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo)
	{
		this->searchAndReturnGeneAnnEntryInfoIndexVec(tmpAnnEntryInfoIndexVec, tmpChrNameInt, tmpPos, indexInfo, 0);
	}

	string searchStartPosAndReturn1stGeneName(int tmpChrNameInt, int tmpPos)
	{
		int tmpAreaNO = (int)(tmpPos/areaSize);
		GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter 
			= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
		if(tmpAreaIter == geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
			return "NULL";
		else
		{
			//cout << "area found !" << endl;
			for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); 
				tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
			{
				int tmpEntryInfoIndex = (*tmpSetIter);
				//cout << "tmpEntryInfoIndex: " << tmpEntryInfoIndex << endl;
				int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpEntryInfoIndex);
				//int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpEntryInfoIndex);
				//cout << "tmpEntry_startPos: " << tmpEntry_startPos << endl;
				//cout << "tmpEntry_endPos: " << tmpEntry_endPos << endl;
				if(tmpEntry_startPos == tmpPos)
					return returnGeneAnnEntry_geneName(tmpEntryInfoIndex);
			}
		}
		return "NULL";	
	}

	string searchEndPosAndReturn1stGeneName(int tmpChrNameInt, int tmpPos)
	{
		int tmpAreaNO = (int)(tmpPos/areaSize);
		GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter 
			= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
		if(tmpAreaIter == geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
			return "NULL";
		else
		{
			//cout << "area found !" << endl;
			for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); 
				tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
			{
				int tmpEntryInfoIndex = (*tmpSetIter);
				//cout << "tmpEntryInfoIndex: " << tmpEntryInfoIndex << endl;
				//int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpEntryInfoIndex);
				int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpEntryInfoIndex);
				//cout << "tmpEntry_startPos: " << tmpEntry_startPos << endl;
				//cout << "tmpEntry_endPos: " << tmpEntry_endPos << endl;
				if(tmpEntry_endPos == tmpPos)
					return returnGeneAnnEntry_geneName(tmpEntryInfoIndex);
			}
		}
		return "NULL";	
	}

	string searchAndReturnSingleGeneId_geneNameVec_strand(string& tmpGeneNameVecStr, string& tmpStrand, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo) // mode = 0, within; mode = 1, boundary only
	{
		vector<int> tmpEntryIndexVec;
		this->searchAndReturnGeneAnnEntryInfoIndexVec(tmpEntryIndexVec, tmpChrNameInt, tmpPos, indexInfo, 0);
		
		vector<string> tmpEntryGeneIdVec;
		vector< pair<int,int> > tmpEntryPosPairVec;
		for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
		{
			int tmpIndex = tmpEntryIndexVec[tmp];
			string tmpEntryStrand = this->returnGeneAnnEntry_strand(tmpIndex);
			if(tmpEntryStrand == tmpStrand)
			{
				string tmpEntry_geneId = this->returnGeneAnnEntry_geneId(tmpIndex);
				int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpIndex);
				int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpIndex);
				tmpEntryGeneIdVec.push_back(tmpEntry_geneId);
				tmpEntryPosPairVec.push_back(pair<int,int>(tmpEntry_startPos, tmpEntry_endPos));
			}
		}
		if(tmpEntryGeneIdVec.size() == 0)
		{
			for(int tmp = 0; tmp < tmpEntryIndexVec.size(); tmp++)
			{
				int tmpIndex = tmpEntryIndexVec[tmp];
				//string tmpEntryStrand = this->returnGeneAnnEntry_strand(tmpIndex);
				//if(tmpEntryStrand == tmpStrand)
				//{
					string tmpEntry_geneId = this->returnGeneAnnEntry_geneId(tmpIndex);
					int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpIndex);
					int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpIndex);
					tmpEntryGeneIdVec.push_back(tmpEntry_geneId);
					tmpEntryPosPairVec.push_back(pair<int,int>(tmpEntry_startPos, tmpEntry_endPos));
				//}
			}			
		}
		vector<string> tmpEntryGeneIdVec_boundaryOnly;
		vector< pair<int,int> > tmpEntryPosPairVec_boundaryOnly;
		for(int tmp = 0; tmp < tmpEntryGeneIdVec.size(); tmp++)
		{
			string tmpEntry_geneId = tmpEntryGeneIdVec[tmp];
			int tmpEntry_startPos = tmpEntryPosPairVec[tmp].first;
			int tmpEntry_endPos = tmpEntryPosPairVec[tmp].second;
			if((tmpEntry_startPos == tmpPos)||(tmpEntry_endPos == tmpPos))
			{
				tmpEntryGeneIdVec_boundaryOnly.push_back(tmpEntry_geneId);
				tmpEntryPosPairVec_boundaryOnly.push_back(pair<int,int>(tmpEntry_startPos, tmpEntry_endPos));
			}
		}
		if(tmpEntryGeneIdVec_boundaryOnly.size() > 0)
			return tmpEntryGeneIdVec_boundaryOnly[0];
		else
			return tmpEntryGeneIdVec[0];
	}

	void searchAndReturnGeneAnnEntryInfoIndexVec(vector<int>& tmpAnnEntryInfoIndexVec, 
		int tmpChrNameInt, int tmpPos, Index_Info* indexInfo, int mode) // mode = 0, within; mode = 1, boundary only
	{
		int tmpAreaNO = (int)(tmpPos/areaSize);
		//cout << "tmpAreaNO: " << tmpAreaNO << endl;
		GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter 
			= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
		if(tmpAreaIter == geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
			return;
		else
		{
			//cout << "area found !" << endl;
			for(set<int>::iterator tmpSetIter = (tmpAreaIter->second).begin(); 
				tmpSetIter != (tmpAreaIter->second).end(); tmpSetIter ++)
			{
				int tmpEntryInfoIndex = (*tmpSetIter);
				//cout << "tmpEntryInfoIndex: " << tmpEntryInfoIndex << endl;
				int tmpEntry_startPos = this->returnGeneAnnEntry_startPos(tmpEntryInfoIndex);
				int tmpEntry_endPos = this->returnGeneAnnEntry_endPos(tmpEntryInfoIndex);
				//cout << "tmpEntry_startPos: " << tmpEntry_startPos << endl;
				//cout << "tmpEntry_endPos: " << tmpEntry_endPos << endl;
				if(mode == 0)
				{	
					if((tmpEntry_startPos <= tmpPos)&&(tmpPos <= tmpEntry_endPos))
						tmpAnnEntryInfoIndexVec.push_back(tmpEntryInfoIndex);
				}
				else if(mode == 1)
				{
					if((tmpEntry_startPos == tmpPos)||(tmpPos == tmpEntry_endPos))
						tmpAnnEntryInfoIndexVec.push_back(tmpEntryInfoIndex);
				}
				else
				{
					cout << "invalid mode " << endl;
					exit(1);
				}
			}
		}
	}

	void loadGeneAnn(string& tmpSimplifiedGeneAnnFile, Index_Info* indexInfo)
	{
		ifstream tmpSimplifiedGeneAnn_ifs(tmpSimplifiedGeneAnnFile.c_str());
		while(!tmpSimplifiedGeneAnn_ifs.eof())
		{
			string tmpStr;
			getline(tmpSimplifiedGeneAnn_ifs, tmpStr);
			if(tmpStr == "")
				break;
			this->insert_newGeneAnnEntry(tmpStr, indexInfo);
		}
		tmpSimplifiedGeneAnn_ifs.close();
	}

	void addNewGeneAnnEntry2infoVec(string& tmpNewGeneAnnEntry, Index_Info* indexInfo)
	{
		vector<string> tmpFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 9; tmp++)
		{
			int tabLoc = tmpNewGeneAnnEntry.find("\t", startLoc);
			string tmpField = tmpNewGeneAnnEntry.substr(startLoc, tabLoc-startLoc);
			tmpFieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		tmpFieldVec.push_back(tmpNewGeneAnnEntry.substr(startLoc));
		string tmpChrName = tmpFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			cout << "invalid chr name, exiting ......" << endl;
			exit(1);
		}		
		string tmpStartPosStr = tmpFieldVec[1];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		string tmpEndPosStr = tmpFieldVec[2];
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpStrand = tmpFieldVec[3];
		string tmpFeatureType = tmpFieldVec[4];
		string tmpGeneId = tmpFieldVec[5];
		string tmpGeneName = tmpFieldVec[6];
		string tmpTranscriptId = tmpFieldVec[7];
		string tmpTranscriptName = tmpFieldVec[8];
		string tmpSource = tmpFieldVec[9];		
		GeneAnnEntry_Info tmpGeneAnnEntryInfo;
		tmpGeneAnnEntryInfo.initiate(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpSource, 
			tmpFeatureType, tmpGeneName, tmpGeneId, tmpTranscriptName, tmpTranscriptId);
		geneAnnEntryInfoVec.push_back(tmpGeneAnnEntryInfo);
	}

	void insert_newGeneAnnEntry(string& tmpNewGeneAnnEntry, Index_Info* indexInfo)
	{
		//string tmpChrName, tmpFeatureType, tmpStrand, tmpGeneId, tmpGeneName, tmpTranscriptId, tmpTranscriptName;
		//int tmpStartPos, tmpEndPos;
		vector<string> tmpFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 9; tmp++)
		{
			int tabLoc = tmpNewGeneAnnEntry.find("\t", startLoc);
			string tmpField = tmpNewGeneAnnEntry.substr(startLoc, tabLoc-startLoc);
			tmpFieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		tmpFieldVec.push_back(tmpNewGeneAnnEntry.substr(startLoc));
		string tmpChrName = tmpFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			cout << "invalid chr name, exiting ......" << endl;
			exit(1);
		}
		
		string tmpStartPosStr = tmpFieldVec[1];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		string tmpEndPosStr = tmpFieldVec[2];
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpStrand = tmpFieldVec[3];
		string tmpFeatureType = tmpFieldVec[4];
		string tmpGeneId = tmpFieldVec[5];
		string tmpGeneName = tmpFieldVec[6];
		string tmpTranscriptId = tmpFieldVec[7];
		string tmpTranscriptName = tmpFieldVec[8];
		string tmpSource = tmpFieldVec[9];
		//cout << "inserting new Entry: " << endl << tmpNewGeneAnnEntry << endl;
		this->insert_newGeneAnnEntry(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpSource,
			tmpFeatureType,	tmpGeneName, tmpGeneId, tmpTranscriptName, tmpTranscriptId);
	}

	void insert_newGeneAnnEntry(string& tmpChrNameStr, int tmpStartPos, int tmpEndPos, 
		string& tmpStrand, string& tmpSource, string& tmpFeatureType, string& tmpGeneName,  
		string& tmpGeneId, string& tmpTranscriptName, string& tmpTranscriptId, Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		this->insert_newGeneAnnEntry(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpSource, 
			tmpFeatureType, tmpGeneName, tmpGeneId, tmpTranscriptName, tmpTranscriptId);
	}

	void insert_newGeneAnnEntry(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, 
		string& tmpStrand, string& tmpSource, string& tmpFeatureType, string& tmpGeneName, 
		string& tmpGeneId, string& tmpTranscriptName, string& tmpTranscriptId)
	{
		GeneAnnEntry_Info tmpGeneAnnEntryInfo;
		tmpGeneAnnEntryInfo.initiate(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpSource, 
			tmpFeatureType, tmpGeneName, tmpGeneId, tmpTranscriptName, tmpTranscriptId);
		geneAnnEntryInfoVec.push_back(tmpGeneAnnEntryInfo);
		int tmpGeneAnnEntryInfoIndex = this->returnGeneAnnEntryNum() - 1;
		//cout << "current geneAnnEntryNum: " << this->returnGeneAnnEntryNum() << endl;
		int areaNO_startPos = (int)(tmpStartPos/areaSize);
		int areaNO_endPos = (int)(tmpEndPos/areaSize);
		//cout << "areaNO_startPos: " << areaNO_startPos << endl;
		//cout << "areaNO_endPos: " << areaNO_endPos << endl;
		for(int tmpAreaNO = areaNO_startPos; tmpAreaNO <= areaNO_endPos; tmpAreaNO ++)
		{
			GeneAnnEntryArea2infoIndexMap::iterator tmpAreaIter;
			tmpAreaIter	= geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].find(tmpAreaNO);
			if(tmpAreaIter != geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].end())
				(tmpAreaIter->second).insert(tmpGeneAnnEntryInfoIndex);
			else
			{
				//cout << "inserting area code " << endl;
				set<int> tmpSet;
				tmpSet.insert(tmpGeneAnnEntryInfoIndex);
				geneAnnEntryArea2infoIndexMapVec[tmpChrNameInt].insert(pair<int, set<int> >(tmpAreaNO, tmpSet));
			}	
		}
	}

	void initiate_geneAnnEntryArea2infoIndexMapVec(int tmpChromTotalNum)
	{
		for(int tmp = 0; tmp < tmpChromTotalNum; tmp ++)
		{
			GeneAnnEntryArea2infoIndexMap tmpGeneAnnEntryArea2infoIndexMap;
			geneAnnEntryArea2infoIndexMapVec.push_back(tmpGeneAnnEntryArea2infoIndexMap);
		}		
	}

	string returnGeneAnnEntry_str(int index, Index_Info* indexInfo)
	{
		return geneAnnEntryInfoVec[index].returnAnnEntry(indexInfo);
	}

	int returnGeneAnnEntry_chrNameInt(int index)
	{
		return geneAnnEntryInfoVec[index].returnChrNameInt();
	}

	int returnGeneAnnEntry_startPos(int index)
	{
		return geneAnnEntryInfoVec[index].returnStartPos();
	}

	int returnGeneAnnEntry_endPos(int index)
	{
		return geneAnnEntryInfoVec[index].returnEndPos();
	}

	string returnGeneAnnEntry_strand(int index)
	{
		return geneAnnEntryInfoVec[index].returnStrand();
	}

	string returnGeneAnnEntry_source(int index)
	{
		return geneAnnEntryInfoVec[index].returnSource();
	}

	string returnGeneAnnEntry_featureType(int index)
	{
		return geneAnnEntryInfoVec[index].returnFeatureType();
	}

	string returnGeneAnnEntry_geneId(int index)
	{
		return geneAnnEntryInfoVec[index].returnGeneId();
	}

	string returnGeneAnnEntry_geneName(int index)
	{
		return geneAnnEntryInfoVec[index].returnGeneName();
	}

	string returnGeneAnnEntry_transcriptId(int index)
	{
		return geneAnnEntryInfoVec[index].returnTranscriptId();
	}

	string returnGeneAnnEntry_transcriptName(int index)
	{
		return geneAnnEntryInfoVec[index].returnTranscriptName();
	}

	int returnGeneAnnEntryNum()
	{
		return geneAnnEntryInfoVec.size();
	}
};

#endif