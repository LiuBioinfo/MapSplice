// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef BASECONTIGUOUSREADCOUNT_INFO_VEC_H
#define BASECONTIGUOUSREADCOUNT_INFO_VEC_H

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

#include "baseContiguousReadCount_info.h"

class BaseContiguousReadCount_Info_Vec
{
private:
	vector<BaseContiguousReadCount_Info*> baseContiguousReadCountInfoVec;
public:
	BaseContiguousReadCount_Info_Vec()
	{}

	int return_contiguousReadCount(int tmpChr, int tmpPos, int tmpEle)
	{
		baseContiguousReadCountInfoVec[tmpEle]->return_contiguousReadCount(tmpChr, tmpPos);
	}

	void initiateBaseContiguousReadCountInfoVec(int vecSize, int chromNum)
	{
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			BaseContiguousReadCount_Info* tmpBaseCountInfo = new BaseContiguousReadCount_Info();
			tmpBaseCountInfo->initaite_chromNum(chromNum);
			baseContiguousReadCountInfoVec.push_back(tmpBaseCountInfo);
		}
	}

	void initiateBase_spliceSiteFromAlignInferJuncHash(AlignInferJunctionHash_Info* juncHash)
	{
		int vecSize = baseContiguousReadCountInfoVec.size();
		for(int tmp = 0; tmp < vecSize; tmp++)
			baseContiguousReadCountInfoVec[tmp]->initiateBase_spliceSiteFromAlignInferJuncHash(juncHash);
	}

	void generatePosAreaMap_withStoredTotalBase()
	{
		int vecSize = baseContiguousReadCountInfoVec.size();
		for(int tmp = 0; tmp < vecSize; tmp++)
			baseContiguousReadCountInfoVec[tmp]->generatePosAreaMap_withStoredTotalBase();
	}

	void update_alignmentFile_parallel(string& alignmentFile, Index_Info* indexInfo)
	{
		int tmpThreadNum = baseContiguousReadCountInfoVec.size();

		ifstream sam_ifs(alignmentFile.c_str());
		int normalRecordNum = 500000;
		vector<string> readAlignmentSAMstrVec(normalRecordNum);
		
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;
		int alignmentTotalNum;
		string tmpLineStr;
		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			realRecordNum = normalRecordNum;
			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				getline(sam_ifs, tmpLineStr);
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				readAlignmentSAMstrVec[recordNumTmp] = tmpLineStr;
			}
			alignmentTotalNum += realRecordNum;
			omp_set_num_threads(tmpThreadNum);
			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			{
				if((readAlignmentSAMstrVec[tmpOpenMP] == "")||(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@'))
					continue;
				int threadNO = omp_get_thread_num();
				//cout << "start to process: " << endl
				//	<< readAlignmentSAMstrVec[tmpOpenMP] << endl;
				baseContiguousReadCountInfoVec[threadNO]->update_contiguousReadCount_withSamStr(
					readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
			}
		}
		sam_ifs.close();		
	}

	void mergeBaseContiguousReadCountInfoVec2one(BaseContiguousReadCount_Info* baseCountInfo_merged)
	{
		int vecSize = baseContiguousReadCountInfoVec.size();
		int totalChromNum = baseCountInfo_merged->returnPos2readCountMapVecSize();
		for(int tmpChr = 0; tmpChr < totalChromNum; tmpChr ++)
		{
			for(map<int,int>::iterator tmpIter = ((baseCountInfo_merged->pos2readCountMapVec)[tmpChr]).begin();
				tmpIter != ((baseCountInfo_merged->pos2readCountMapVec)[tmpChr]).end(); tmpIter ++)
			{
				int tmpPos = (tmpIter->first);
				int tmpReadCount = (tmpIter->second);
				for(int tmpEle = 0; tmpEle < vecSize; tmpEle++)
				{
					int tmpReadCount_inTmpEle = this->return_contiguousReadCount(tmpChr, tmpPos, tmpEle);
					tmpReadCount += tmpReadCount_inTmpEle;
				}
				(tmpIter->second) = tmpReadCount;
			}
		}
	}

	void freeMemory()
	{
		for(int tmp = 0; tmp < baseContiguousReadCountInfoVec.size(); tmp++)
			delete baseContiguousReadCountInfoVec[tmp];
	}
};
#endif