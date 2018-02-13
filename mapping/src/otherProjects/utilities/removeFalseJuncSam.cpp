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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"

time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputFilteredOutJuncFile inputSam outputSam threads_num" << endl;
		exit(1);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to remove alignments containing false junctions ......" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();

	AlignInferJunctionHash_Info* juncHashInfo = new AlignInferJunctionHash_Info();
	juncHashInfo->initiateAlignInferJunctionHashInfo(chromNum);

	string inputFilteredOutJuncFile = argv[2];
	int maximumSupNum = 2;
	juncHashInfo->insertJuncFromJuncFile_chrNamePosOnly_lowSupOnly(inputFilteredOutJuncFile, maximumSupNum, indexInfo);

	string threads_num_str = argv[5];
	int threads_num = atoi(threads_num_str.c_str());
	string inputSam = argv[3];
	string outputSam = argv[4];	
	ifstream sam_ifs(inputSam.c_str());
	ofstream sam_ofs(outputSam.c_str());
	int normalRecordNum = 2000000;
	vector<string> readAlignmentSAMstrVec(normalRecordNum);
	vector<string> refinedSAMstr2outputVec(normalRecordNum);
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
		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
		{
			if(readAlignmentSAMstrVec[tmpOpenMP] == "")
				continue;
			else if(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@')
				refinedSAMstr2outputVec[tmpOpenMP] = readAlignmentSAMstrVec[tmpOpenMP];
			else
			{
				bool juncFound_bool = juncHashInfo->searchJuncFromSAMstr(readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
				if(juncFound_bool)
					refinedSAMstr2outputVec[tmpOpenMP] = "";
				else
					refinedSAMstr2outputVec[tmpOpenMP] = readAlignmentSAMstrVec[tmpOpenMP];
			}
		}
		for(int recordNumTmp = 0; recordNumTmp < realRecordNum; recordNumTmp++)
		{
			string tmpSamStr = refinedSAMstr2outputVec[recordNumTmp];
			if(tmpSamStr != "")
				sam_ofs << tmpSamStr << endl;
		}
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "end of removing alignments containing false junctions ......" << endl;
	delete juncHashInfo;
	parameter_ifs.close();
	delete indexInfo;
	sam_ifs.close();
	sam_ofs.close();
	return 0;
}