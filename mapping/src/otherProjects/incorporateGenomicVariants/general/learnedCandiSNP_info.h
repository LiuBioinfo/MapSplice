#ifndef LEARNEDCANDISNP_INFO_H
#define LEARNEDCANDISNP_INFO_H

#include "stdio.h"
#include "stdlib.h"

using namespace std;

int fiveBaseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 4,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};

class LearnedCandiSNP_Info
{
private:
	int chrNameInt;
	int pos;
	string referBase;

	int totalMismatchNum;
	vector<int> mismatchReadCountVec;
public:
	LearnedCandiSNP_Info()
	{
		mismatchReadCountVec.push_back(0); // A
		mismatchReadCountVec.push_back(0); // C
		mismatchReadCountVec.push_back(0); // G
		mismatchReadCountVec.push_back(0); // T
		mismatchReadCountVec.push_back(0); // N
	}

	void mergedWithLearnedCandiSNPinfoFromAnotherHash(int tmpTotalMismatchNum, int tmpMismatchNum_A, 
		int tmpMismatchNum_C, int tmpMismatchNum_G, int tmpMismatchNum_T, int tmpMismatchNum_N)
	{
		totalMismatchNum += tmpTotalMismatchNum;
		mismatchReadCountVec[0] += tmpMismatchNum_A;
		mismatchReadCountVec[1] += tmpMismatchNum_C;
		mismatchReadCountVec[2] += tmpMismatchNum_G;
		mismatchReadCountVec[3] += tmpMismatchNum_T;
		mismatchReadCountVec[4] += tmpMismatchNum_N;
	}

	int returnTotalMismatchNum()
	{
		return totalMismatchNum;
	}

	int returnMismatchReadCountForSpecificBaseIndex(
		int BaseIndex) // 0--A, 1--C, 2--G, 3--T, 4--N
	{
		return mismatchReadCountVec[BaseIndex];
	}

	string learnedCandiSNPinfoStr(Index_Info* indexInfo)
	{
		string tmpChrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string tmpPosStr = int_to_str(pos);
		string totalMismatchNumStr = int_to_str(totalMismatchNum);
		vector<string> mismatchReadCountStrVec;
		mismatchReadCountStrVec.push_back(int_to_str(mismatchReadCountVec[0]));
		mismatchReadCountStrVec.push_back(int_to_str(mismatchReadCountVec[1]));
		mismatchReadCountStrVec.push_back(int_to_str(mismatchReadCountVec[2]));
		mismatchReadCountStrVec.push_back(int_to_str(mismatchReadCountVec[3]));
		mismatchReadCountStrVec.push_back(int_to_str(mismatchReadCountVec[4]));
		string tmpStr = tmpChrNameStr + "\t" + tmpPosStr + "\t" + referBase + "\t"
			+ totalMismatchNumStr + "\t" + mismatchReadCountStrVec[0] + "\t" + mismatchReadCountStrVec[1]
			+ "\t" + mismatchReadCountStrVec[2] + "\t" + mismatchReadCountStrVec[3] + "\t" + mismatchReadCountStrVec[4];
		return tmpStr;
	}

	void initiate_1stMismatch(int tmpChrNameInt, int tmpChrPos, 
		string& mismatchStrBase, Index_Info* indexInfo)
	{
		chrNameInt = tmpChrNameInt;
		pos = tmpChrPos;
		referBase = indexInfo->returnChromStrSubstr(chrNameInt, tmpChrPos, 1);
		char mismatchCharBase = mismatchStrBase[0];
		int tmpIndex = fiveBaseChar2intArray[mismatchCharBase - 'A'];
		mismatchReadCountVec[tmpIndex] = 1;
		totalMismatchNum = 1;
	}

	void initiate_withLearnedCandiSNPinfoFromAnotherHash(int tmpChrNameInt, 
		int tmpChrPos, int tmpTotalMismatchNum, int tmpMismatchNum_A, 
		int tmpMismatchNum_C, int tmpMismatchNum_G, int tmpMismatchNum_T, 
		int tmpMismatchNum_N, Index_Info* indexInfo)
	{
		chrNameInt = tmpChrNameInt;
		pos = tmpChrPos;
		referBase = indexInfo->returnChromStrSubstr(chrNameInt, tmpChrPos, 1);		
		totalMismatchNum = tmpTotalMismatchNum;
		mismatchReadCountVec[0] = tmpMismatchNum_A;
		mismatchReadCountVec[1] = tmpMismatchNum_C;
		mismatchReadCountVec[2] = tmpMismatchNum_G;
		mismatchReadCountVec[3] = tmpMismatchNum_T;
		mismatchReadCountVec[4] = tmpMismatchNum_N;
	}

	void addMismatchReadCount(string& mismatchStrBase)
	{
		char mismatchCharBase = mismatchStrBase[0];
		int tmpIndex = fiveBaseChar2intArray[mismatchCharBase - 'A'];
		mismatchReadCountVec[tmpIndex] ++;
		totalMismatchNum ++;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnPos()
	{
		return pos;
	}

	string returnReferBase()
	{
		return referBase;
	}

	string returnSNPinfoStr(Index_Info* indexInfo)
	{
		string tmpSNPstr;
		string tmpChrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string tmpPosStr = int_to_str(pos);
		tmpSNPstr = tmpChrNameStr + "\t" + tmpPosStr + "\t" + referBase
			+ "\t" + int_to_str(mismatchReadCountVec[0]) + "\t" + int_to_str(mismatchReadCountVec[1])
			+ "\t" + int_to_str(mismatchReadCountVec[2]) + "\t" + int_to_str(mismatchReadCountVec[3])
			+ "\t" + int_to_str(mismatchReadCountVec[4]);
		return tmpSNPstr;
	}
};


#endif