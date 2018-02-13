#ifndef LEARNEDCANDISNPHASH_INFO_H
#define LEARNEDCANDISNPHASH_INFO_H

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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../general/learnedCandiSNP_info.h"

using namespace std;

typedef map<int, int> LearnedCandiSNPinfoIndexMap;
typedef map<int, set<int> > LearnedCandiSNPareaMap; //( areaNO = pos/100 ) intermediate hash to directly get all possible SNPs near a certain position

class LearnedCandiSNPhash_Info
{
private:
	vector <LearnedCandiSNPareaMap> learnedCandiSNPareaMapVec;
	vector <LearnedCandiSNPinfoIndexMap> learnedCandiSNPinfoIndexMapVec;

	int areaSize;
public:
	vector<LearnedCandiSNP_Info> learnedCandiSNPinfoVec;

	LearnedCandiSNPhash_Info()
	{
		areaSize = 1000;
	}

	void initiate(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			LearnedCandiSNPareaMap tmpLearnedCandiSNPareaMap;
			learnedCandiSNPareaMapVec.push_back(tmpLearnedCandiSNPareaMap);
			LearnedCandiSNPinfoIndexMap tmpLearnedCandiSNPinfoIndexMap;
			learnedCandiSNPinfoIndexMapVec.push_back(tmpLearnedCandiSNPinfoIndexMap);
		}
	}

	int returnLearnedCandiSNPinfoVecSize()
	{
		return learnedCandiSNPinfoVec.size();
	}

	void mergeWithAnotherLearnedCandiSNPhash(LearnedCandiSNPhash_Info* theOtherLearnedCandiSNPhashInfo, 
		Index_Info* indexInfo)
	{
		int tmpSNPinfoVecSize = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNPinfoVecSize();
		for(int tmp = 0; tmp < tmpSNPinfoVecSize; tmp++)
		{
			int tmpChrNameInt = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_chrNameInt(tmp);
			int tmpPos = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_pos(tmp);
			string tmpReferBase = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_referBase(tmp);
			int tmpTotalMismatchNum = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_totalMismatchNum(tmp);
			int tmpMismatchReadCount_A = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(tmp, 0);
			int tmpMismatchReadCount_C = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(tmp, 1);
			int tmpMismatchReadCount_G = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(tmp, 2);
			int tmpMismatchReadCount_T = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(tmp, 3);
			int tmpMismatchReadCount_N = theOtherLearnedCandiSNPhashInfo->returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(tmp, 4);
			LearnedCandiSNPinfoIndexMap::iterator tmpIter = learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].find(tmpPos);
			if(tmpIter != learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].end()) // found
			{
				int foundIndex_in_learnedCandiSNPinfoVec = tmpIter->second;
				learnedCandiSNPinfoVec[foundIndex_in_learnedCandiSNPinfoVec].mergedWithLearnedCandiSNPinfoFromAnotherHash(tmpTotalMismatchNum,
					tmpMismatchReadCount_A, tmpMismatchReadCount_C, tmpMismatchReadCount_G, tmpMismatchReadCount_T, tmpMismatchReadCount_N);
			}
			else // not found
			{
				this->insertLearnedCandiSNPpos2areaHash(tmpChrNameInt, tmpPos);
				LearnedCandiSNP_Info tmpLearnedCandiSNPinfo;
				tmpLearnedCandiSNPinfo.initiate_withLearnedCandiSNPinfoFromAnotherHash(tmpChrNameInt, 
					tmpPos, tmpTotalMismatchNum, tmpMismatchReadCount_A, tmpMismatchReadCount_C, 
					tmpMismatchReadCount_G, tmpMismatchReadCount_T, tmpMismatchReadCount_N, indexInfo);
				learnedCandiSNPinfoVec.push_back(tmpLearnedCandiSNPinfo);
				int tmpNewAddedLearnedCandiSNPinfoIndex = learnedCandiSNPinfoVec.size() - 1;
				learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].insert(pair<int,int>(tmpPos, tmpNewAddedLearnedCandiSNPinfoIndex));
			}
		}
	}

	int returnLearnedCandiSNP_chrNameInt(int tmp)
	{
		return learnedCandiSNPinfoVec[tmp].returnChrNameInt();
	}

	int returnLearnedCandiSNP_pos(int tmp)
	{
		return learnedCandiSNPinfoVec[tmp].returnPos();
	}

	string returnLearnedCandiSNP_referBase(int tmp)
	{
		return learnedCandiSNPinfoVec[tmp].returnReferBase();
	}

	int returnLearnedCandiSNP_totalMismatchNum(int tmp)
	{
		return learnedCandiSNPinfoVec[tmp].returnTotalMismatchNum();
	}

	int returnLearnedCandiSNP_mismatchReadCountForSpecificBaseIndex(int tmp, int baseIndex)
	{
		return learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(baseIndex);
	}

	void insertLearnedCandiSNPpos2areaHash(int tmpChrNameInt, int tmpPos)
	{
		int tmpAreaNO = (int)(tmpPos/areaSize);		
		LearnedCandiSNPareaMap::iterator tmpAreaMapIter;
		tmpAreaMapIter = learnedCandiSNPareaMapVec[tmpChrNameInt].find(tmpAreaNO);
		if(tmpAreaMapIter == learnedCandiSNPareaMapVec[tmpChrNameInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(tmpPos);
			learnedCandiSNPareaMapVec[tmpChrNameInt].insert(pair<int, set<int> >(tmpAreaNO, newPosSet));
		}
		else
		{
			if((tmpAreaMapIter->second).find(tmpPos) == (tmpAreaMapIter->second).end())
				(tmpAreaMapIter->second).insert(tmpPos);
			else
			{}
		}
	}

	void addMismatchInSam(int tmpChrNameInt, int tmpPos, string& tmpMismatchBase, Index_Info* indexInfo)
	{
		LearnedCandiSNPinfoIndexMap::iterator tmpIter = learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].find(tmpPos);
		if(tmpIter != learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].end()) // found
		{
			int foundIndex_in_learnedCandiSNPinfoVec = tmpIter->second;
			learnedCandiSNPinfoVec[foundIndex_in_learnedCandiSNPinfoVec].addMismatchReadCount(tmpMismatchBase);
		}
		else // not found
		{
			this->insertLearnedCandiSNPpos2areaHash(tmpChrNameInt, tmpPos);
			LearnedCandiSNP_Info tmpLearnedCandiSNPinfo;
			tmpLearnedCandiSNPinfo.initiate_1stMismatch(tmpChrNameInt, tmpPos, tmpMismatchBase, indexInfo);
			learnedCandiSNPinfoVec.push_back(tmpLearnedCandiSNPinfo);
			int tmpNewAddedLearnedCandiSNPinfoIndex = learnedCandiSNPinfoVec.size() - 1;
			learnedCandiSNPinfoIndexMapVec[tmpChrNameInt].insert(pair<int,int>(tmpPos, tmpNewAddedLearnedCandiSNPinfoIndex));
		}
	}

	int returnLearnedCandiSNP_chrPos(int index)
	{
		return learnedCandiSNPinfoVec[index].returnPos();
	}

	void outputCandiSNP(string& outputFile, Index_Info* indexInfo)
	{
		ofstream candiSNP_ofs(outputFile.c_str());
		int candiSNPnum = learnedCandiSNPinfoVec.size();
		for(int tmp = 0; tmp < candiSNPnum; tmp++)
		{
			string tmpCandiSNPstr = learnedCandiSNPinfoVec[tmp].learnedCandiSNPinfoStr(indexInfo);
			candiSNP_ofs << tmpCandiSNPstr << endl;
		}
		candiSNP_ofs.close();
	}

	void output_filteredLearnedSNP(string& filteredLearnedSNP_file, int& called_snp_num, int filterCandiSNP_supNumMin, 
		double filterCandiSNP_ratioMin, Index_Info* indexInfo)
	{
		int tmp_called_SNP_num = 0;
		ofstream SNP_ofs(filteredLearnedSNP_file.c_str());
		int candiSNPnum = learnedCandiSNPinfoVec.size();
		for(int tmp = 0; tmp < candiSNPnum; tmp++)
		{
			int tmpChrNameInt = learnedCandiSNPinfoVec[tmp].returnChrNameInt();
			string tmpChrName = indexInfo->returnChrNameStr(tmpChrNameInt);
			int tmpChrPos = learnedCandiSNPinfoVec[tmp].returnPos();
			string tmpRefBase = learnedCandiSNPinfoVec[tmp].returnReferBase();
			int tmpTotalMismatchNum = learnedCandiSNPinfoVec[tmp].returnTotalMismatchNum();
			int tmpMismatchNum_A = learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(0);
			int tmpMismatchNum_C = learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(1);
			int tmpMismatchNum_G = learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(2);
			int tmpMismatchNum_T = learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(3);
			int tmpMismatchNum_N = learnedCandiSNPinfoVec[tmp].returnMismatchReadCountForSpecificBaseIndex(4);

			double tmpMismatchRatio_A = (double)tmpMismatchNum_A / (double)tmpTotalMismatchNum;
			double tmpMismatchRatio_C = (double)tmpMismatchNum_C / (double)tmpTotalMismatchNum;
			double tmpMismatchRatio_G = (double)tmpMismatchNum_G / (double)tmpTotalMismatchNum;
			double tmpMismatchRatio_T = (double)tmpMismatchNum_T / (double)tmpTotalMismatchNum;			
			
			if((tmpMismatchNum_A >= filterCandiSNP_supNumMin)&&(tmpMismatchRatio_A >= filterCandiSNP_ratioMin))
			{
				SNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tA" << endl;
				tmp_called_SNP_num ++;
			}
			else if((tmpMismatchNum_C >= filterCandiSNP_supNumMin)&&(tmpMismatchRatio_C >= filterCandiSNP_ratioMin))
			{
				SNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tC" << endl;
				tmp_called_SNP_num ++;
			}
			else if((tmpMismatchNum_G >= filterCandiSNP_supNumMin)&&(tmpMismatchRatio_G >= filterCandiSNP_ratioMin))
			{
				SNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tG" << endl;
				tmp_called_SNP_num ++;
			}
			else if((tmpMismatchNum_T >= filterCandiSNP_supNumMin)&&(tmpMismatchRatio_T >= filterCandiSNP_ratioMin))
			{
				SNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tT" << endl;
				tmp_called_SNP_num ++;
			}			
			else
			{}
		}
		called_snp_num = tmp_called_SNP_num;
		SNP_ofs.close();
	}
};

#endif