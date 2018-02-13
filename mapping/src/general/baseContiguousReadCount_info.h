// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef BASECONTIGUOUSREADCOUNT_INFO_H
#define BASECONTIGUOUSREADCOUNT_INFO_H

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

typedef map<int,int> Pos2readCountMap;

typedef map<int, set<int> > PosAreaMap;

class BaseContiguousReadCount_Info
{
private:
	vector<PosAreaMap> posAreaMapVec;
	int areaSize;
public:
	vector<Pos2readCountMap> pos2readCountMapVec;

	BaseContiguousReadCount_Info()
	{
		areaSize = 1000;
	}

	bool searchBaseWithinRegion_returnBasePosVec(int tmpSam_chrNameInt, int tmpSam_tmpStartPosInChr, 
		int tmpSam_tmpEndPosInChr, vector<int>& tmpSam_tmpBasePosVec)
	{
		int areaNOmin = (int)(tmpSam_tmpStartPosInChr/areaSize);
		int areaNOmax = (int)(tmpSam_tmpEndPosInChr/areaSize);
		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			PosAreaMap::iterator tmpPosAreaIter = posAreaMapVec[tmpSam_chrNameInt].find(tmpArea);
			if(tmpPosAreaIter != posAreaMapVec[tmpSam_chrNameInt].end()) // found
			{
				for(set<int>::iterator tmpIntSetIter = (tmpPosAreaIter->second).begin();
					tmpIntSetIter != (tmpPosAreaIter->second).end(); tmpIntSetIter ++)
				{
					int tmpCandiPos = (*tmpIntSetIter);
					if((tmpCandiPos >= tmpSam_tmpStartPosInChr)&&(tmpCandiPos + 1 <= tmpSam_tmpEndPosInChr))
						tmpSam_tmpBasePosVec.push_back(tmpCandiPos);
				}
			}
		}
	}	

	int returnPos2readCountMapVecSize()
	{
		return pos2readCountMapVec.size();
	}

	void initaite_chromNum(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			Pos2readCountMap tmpPos2readCountMap;
			pos2readCountMapVec.push_back(tmpPos2readCountMap);
			PosAreaMap tmpPosAreaMap;
			posAreaMapVec.push_back(tmpPosAreaMap);
		}
	}

	void generatePosAreaMap_withStoredTotalBase()
	{
		for(int tmpChr = 0; tmpChr < pos2readCountMapVec.size(); tmpChr ++)
		{
			for(Pos2readCountMap::iterator tmpIter = pos2readCountMapVec[tmpChr].begin();
				tmpIter != pos2readCountMapVec[tmpChr].end(); tmpIter ++)
			{
				int tmpPos = tmpIter->first;
				int tmpPosAreaNO = (int)(tmpPos/areaSize);
				PosAreaMap::iterator tmpAreaMapIter;
				tmpAreaMapIter = posAreaMapVec[tmpChr].find(tmpPosAreaNO);
				if(tmpAreaMapIter == posAreaMapVec[tmpChr].end())
				{
					set<int> newPosSet;
					newPosSet.insert(tmpPos);
					posAreaMapVec[tmpChr].insert(pair<int, set<int> >(tmpPosAreaNO, newPosSet));
				}
				else
				{
					if((tmpAreaMapIter->second).find(tmpPos) == (tmpAreaMapIter->second).end())
						(tmpAreaMapIter->second).insert(tmpPos);
					else
					{}
				}
			}
		}		
	}


	int getEndLocInReadOfSpecificJumpCode(
		vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		if(jumpCodeIndex < 0)
			return 0;
		int tmpEndLocInRead = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "D")
			{
			}
			else if(tmpJumpCodeType == "N")
			{
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return tmpEndLocInRead;
	}	

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{}
			else if(tmpJumpCodeType == "M")
				tmpEndPos += tmpJumpCodeLength;
			else if(tmpJumpCodeType == "I")
			{}
			else if(tmpJumpCodeType == "D")
				tmpEndPos += tmpJumpCodeLength;
			else if(tmpJumpCodeType == "N")
				tmpEndPos += tmpJumpCodeLength;
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	void generateExonLocInReadPosInChr(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector<int>& endLocVecInRead, vector<int>& endPosVecInChr, vector<int>& lenVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp ++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmp].len;
			int tmpJumpCodeIndex = tmp;
			if(tmpJumpCodeType == "M")
			{
				int tmpEndLocInRead = getEndLocInReadOfSpecificJumpCode(cigarStringJumpCodeVec, tmpJumpCodeIndex);
				int tmpEndPosInChr = getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec, tmpJumpCodeIndex);
				endLocVecInRead.push_back(tmpEndLocInRead);
				endPosVecInChr.push_back(tmpEndPosInChr);
				lenVec.push_back(tmpJumpCodeLength);
			}
		}
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
			
		string candidateJumpCodeType = "SMNIDXH";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				if(tmpJumpCodeType != "H")
					cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	bool parseSam2chrNamePosCigarString(string& tmpSamStr, int& tmpSam_chrNameInt, 
		int& tmpSam_chrMapPos, string& tmpSam_cigarString, Index_Info* indexInfo)
	{
		int tabLoc_1 = tmpSamStr.find("\t");
		int tabLoc_2 = tmpSamStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpSamStr.find("\t", tabLoc_2 + 1);		
		int tabLoc_4 = tmpSamStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpSamStr.find("\t", tabLoc_4 + 1);
		int tabLoc_6 = tmpSamStr.find("\t", tabLoc_5 + 1);
		string tmpSam_tmpChrName = tmpSamStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		tmpSam_chrNameInt = indexInfo->convertStringToInt(tmpSam_tmpChrName);
		if(tmpSam_chrNameInt < 0)
			return false;
		string tmpSam_tmpChrPosStr = tmpSamStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		tmpSam_chrMapPos = atoi(tmpSam_tmpChrPosStr.c_str());
		if(tmpSam_chrMapPos < 0)
			return false;
		tmpSam_cigarString = tmpSamStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		if(tmpSam_cigarString == "*")
			return false;
		return true;
	}

	void update_contiguousReadCount_withSamStr(string& tmpSamStr, Index_Info* indexInfo)
	{
		int tmpSam_chrNameInt, tmpSam_chrMapPos;
		string tmpSam_cigarString, tmpSam_readSeq;
		//cout << "start to parse " << endl;
		bool tmpSam_parseBool = this->parseSam2chrNamePosCigarString(tmpSamStr, tmpSam_chrNameInt, 
			tmpSam_chrMapPos, tmpSam_cigarString, indexInfo);
		//cout << "parseBool: " << parseBool << endl;
		if(!tmpSam_parseBool)
			return;
		vector<Jump_Code> tmpSam_cigarStringJumpCodeVec;
		this->cigarString2jumpCodeVec(tmpSam_cigarString, tmpSam_cigarStringJumpCodeVec);
		int tmpSam_cigarStringJumpCodeVecSize = tmpSam_cigarStringJumpCodeVec.size();
		int tmpSam_chrMapPos_end = this->getEndPosOfSpecificJumpCode(tmpSam_chrMapPos, 
				tmpSam_cigarStringJumpCodeVec, tmpSam_cigarStringJumpCodeVecSize - 1);
		vector<int> tmpSam_endLocVecInRead;
		vector<int> tmpSam_endPosVecInChr;
		vector<int> tmpSam_lenVec;
		//cout << "start to do generateExonLocInReadPosInChr " << endl;
		this->generateExonLocInReadPosInChr(tmpSam_chrMapPos, tmpSam_cigarStringJumpCodeVec, 
			tmpSam_endLocVecInRead, tmpSam_endPosVecInChr, tmpSam_lenVec);		
		for(int tmp = 0; tmp < tmpSam_endLocVecInRead.size(); tmp++)
		{
			int tmpSam_tmpEndLocInRead = tmpSam_endLocVecInRead[tmp];
			int tmpSam_tmpEndPosInChr = tmpSam_endPosVecInChr[tmp];
			int tmpSam_tmpLen = tmpSam_lenVec[tmp];
			int tmpSam_tmpStartPosInChr = tmpSam_tmpEndPosInChr - tmpSam_tmpLen + 1;
			vector<int> tmpSam_tmpBasePosVec;
			bool tmp_baseExistsWithinRegion_bool = this->searchBaseWithinRegion_returnBasePosVec(tmpSam_chrNameInt,
				tmpSam_tmpStartPosInChr, tmpSam_tmpEndPosInChr, tmpSam_tmpBasePosVec);
			if(tmp_baseExistsWithinRegion_bool)
			{
				int tmpSam_tmpBasePosVecSize = tmpSam_tmpBasePosVec.size();
				for(int tmpBaseIndex = 0; tmpBaseIndex < tmpSam_tmpBasePosVecSize; tmpBaseIndex ++)
				{
					int tmpBasePos = tmpSam_tmpBasePosVec[tmpBaseIndex];
					Pos2readCountMap::iterator tmpIter = pos2readCountMapVec[tmpSam_chrNameInt].find(tmpBasePos);
					if(tmpIter != pos2readCountMapVec[tmpSam_chrNameInt].end())
						(tmpIter->second)++;
					else
					{
						cout << "error ! base should be in hash" << endl;
						exit(1);
					}
				}
			}
		}	
	}

	int return_contiguousReadCount(int tmpChrNameInt, int tmpPos)
	{
		Pos2readCountMap::iterator tmpIter = pos2readCountMapVec[tmpChrNameInt].find(tmpPos);
		if(tmpIter == pos2readCountMapVec[tmpChrNameInt].end())
			return 0;
		else
			return (tmpIter->second);
	}

	void insert_initiate_base(int tmpChrNameInt, int tmpPos)
	{
		pos2readCountMapVec[tmpChrNameInt].insert(pair<int,int>(tmpPos, 0));
	}

	void initiateBase_spliceSiteFromAlignInferJuncHash(AlignInferJunctionHash_Info* juncHash)
	{
		int juncNum = juncHash->returnAlignInferInfoVecSize();
		for(int tmpJunc = 0; tmpJunc < juncNum; tmpJunc ++)
		{
			int tmpJunc_chrNameInt = juncHash->returnAlignInferInfo_chrNameInt(tmpJunc);
			int tmpJunc_donerEndPos = juncHash->returnAlignInferInfo_donerEndPos(tmpJunc);
			int tmpJunc_acceptorStartPos = juncHash->returnAlignInferInfo_acceptorStartPos(tmpJunc);
			int toInsertPos_1 = tmpJunc_donerEndPos;
			int toInsertPos_2 = tmpJunc_acceptorStartPos - 1;
			// insert two bases
			this->insert_initiate_base(tmpJunc_chrNameInt, toInsertPos_1);
			this->insert_initiate_base(tmpJunc_chrNameInt, toInsertPos_2);
		}
	}
};
#endif