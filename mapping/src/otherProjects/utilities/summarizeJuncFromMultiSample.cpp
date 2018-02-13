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

time_t nowtime;
struct tm *local;

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable inputIndexPath outputSummaryFile inputSJ_1 inputSJ_2 ..." << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();

	int juncFileNum = argc - 3;
	vector<string> juncFileVec;
	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < juncFileNum; tmp++)
	{
		string tmpJuncFile = argv[3 + tmp];
		juncFileVec.push_back(tmpJuncFile);
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionHashInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum_includeBackSplice(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}
	AlignInferJunctionHash_Info* juncHash_merged = new AlignInferJunctionHash_Info();
	juncHash_merged->initiateAlignInferJunctionHashInfo(chromNum);
	juncHash_merged->insertJuncFromJuncFileVec_chrNamePosOnly(juncFileVec, indexInfo);

	string outputSummaryFile = argv[2];
	ofstream sum_ofs(outputSummaryFile.c_str()); 
	int totalJuncNum = juncHash_merged->returnAlignInferInfoVecSize();
	for(int tmpJunc = 0; tmpJunc < totalJuncNum; tmpJunc++)
	{
		int tmpSJ_chrNameInt = juncHash_merged->returnAlignInferInfo_chrNameInt(tmpJunc);
		string tmpSJ_chrName = indexInfo->returnChrNameStr(tmpSJ_chrNameInt);
		int tmpSJ_donerEndPos = juncHash_merged->returnAlignInferInfo_donerEndPos(tmpJunc);
		int tmpSJ_acceptorStartPos = juncHash_merged->returnAlignInferInfo_acceptorStartPos(tmpJunc);
		int tmpSJ_totalSupNum = 0;
		vector<int> tmpSJ_supNumVec;
		for(int tmpSample = 0; tmpSample < juncFileNum; tmpSample ++)
		{
			int tmpSample_supNum = juncHashVec[tmpSample]->searchAndReturnAlignInferJuncHashSupNum(
				tmpSJ_chrNameInt, tmpSJ_donerEndPos, tmpSJ_acceptorStartPos);
			tmpSJ_supNumVec.push_back(tmpSample_supNum);
			tmpSJ_totalSupNum += tmpSample_supNum;
		}
		if(tmpSJ_totalSupNum > 0)
		{
			sum_ofs << tmpSJ_chrName << "\t" << tmpSJ_donerEndPos << "\t" << tmpSJ_acceptorStartPos 
				<< "\t" << tmpSJ_totalSupNum;
			for(int tmpSample = 0; tmpSample < juncFileNum; tmpSample ++)
			{
				int tmpSJ_tmpSampleSupNum = tmpSJ_supNumVec[tmpSample];
				sum_ofs << "\t" << tmpSJ_tmpSampleSupNum;
			}
			sum_ofs << endl;
		}
	}
	sum_ofs.close();
	delete indexInfo;
	parameter_ifs.close();
	delete juncHash_merged;
	for(int tmp = 0; tmp < juncFileNum; tmp++)
		delete juncHashVec[tmp];
	return 0;
}