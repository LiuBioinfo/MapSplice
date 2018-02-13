// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef LEARNEDCANDISNPHASH_INFO_VEC_H
#define LEARNEDCANDISNPHASH_INFO_VEC_H

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

#include "../general/learnedCandiSNPhash_info.h"

using namespace std;

class LearnedCandiSNPhash_Info_Vec
{
private:
	vector<LearnedCandiSNPhash_Info*> learnedCandiSNPhashInfoVec;

public:
	LearnedCandiSNPhash_Info_Vec()
	{}

	void freeMemory()
	{
		for(int tmp = 0; tmp < learnedCandiSNPhashInfoVec.size(); tmp++)
			delete learnedCandiSNPhashInfoVec[tmp];
	}

	void initiateLearnedCandiSNPhashInfoVec(int learnedCandiSNPhashInfoVecSize, int tmpChromNum)
	{
		for(int tmp = 0; tmp < learnedCandiSNPhashInfoVecSize; tmp++)
		{
			LearnedCandiSNPhash_Info* tmpLearnedCandiSNPhashInfo = new LearnedCandiSNPhash_Info();
			tmpLearnedCandiSNPhashInfo->initiate(tmpChromNum);
			learnedCandiSNPhashInfoVec.push_back(tmpLearnedCandiSNPhashInfo);
		}
	}

	void addMismatchInSam(int tmpChrNameInt, int tmpPos, string& tmpMismatchBase, 
		Index_Info* indexInfo, int tmpIndexInLearnedCandiSNPhashInfoVec)
	{
		learnedCandiSNPhashInfoVec[tmpIndexInLearnedCandiSNPhashInfoVec]->addMismatchInSam(
			tmpChrNameInt, tmpPos, tmpMismatchBase, indexInfo);
	}

	void mergeLearnedCandiSNPhashInfoVec2one(LearnedCandiSNPhash_Info& mergedLearnedCandiSNPhashInfo, 
		Index_Info* indexInfo)
	{
		int tmpVecSize = learnedCandiSNPhashInfoVec.size();
		for(int tmp = 0; tmp < tmpVecSize; tmp++)
			mergedLearnedCandiSNPhashInfo.mergeWithAnotherLearnedCandiSNPhash(
				learnedCandiSNPhashInfoVec[tmp], indexInfo);
	}
};
#endif