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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/alignInferJunctionHash_info_vec.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <threadNum> <outputFolder> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc.txt";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string threadNumStr = argv[2];
	int threadNum_int = atoi(threadNumStr.c_str());
	int alignInferJuncHashInfoVecSize = threadNum_int;
	vector<string> inputSAMfileVec;
	for(int tmp = 4; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}		

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;	
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();		
	indexInfo->initiateChrNameIndexArray(1000);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of initiating chrNameIndexArray" << endl;
	log_ofs << endl << "[" << asctime(local) << "... end of initiating chrNameIndexArray" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initaite merged alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...start to initaite merged alignInferJunctionHashInfo " << endl;	
	int chromNum = indexInfo->returnChromNum();
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_merged = new AlignInferJunctionHash_Info();	
	alignInferJunctionHashInfo_merged->initiateAlignInferJunctionHashInfo(chromNum);
	AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec = new AlignInferJunctionHash_Info_Vec();
	alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(alignInferJuncHashInfoVecSize, chromNum);

	// insert SJ into SJmap
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to insert SJ into SJmap from SAM file" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to insert SJ into SJmap from SAM file" << endl;
	alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize_XM_parallel(
		inputSAMfileVec, indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to merge alignInferJuncHashInfo in vec" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to merge alignInferJuncHashInfo in vec" << endl;
	alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum_anchorSize_XM(
		alignInferJunctionHashInfo_merged, indexInfo);

	

	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	free(chrom);
	delete indexInfo;
	return 0;
}