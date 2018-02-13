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

//#include "../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"
//#include "../general/alignInferJunctionHash_info_vec.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable inputIndexInfo outputMergedFile juncFileNum juncFile_1 (juncFile_2 ...)" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string outputMergedFile = argv[2];
	string juncFileNumStr = argv[3];
	int juncFileNum = atoi(juncFileNumStr.c_str());
	int juncFile_num_in_commandLine = argc - 4;
	if(juncFileNum != juncFile_num_in_commandLine)
	{
		cout << "juncFileNum != juncFile_num_in_commandLine" << endl;
		exit(1);
	}
	vector<string> juncFileVec;
	for(int tmp = 4; tmp < argc; tmp++)
	{
		string tmpJuncFile = argv[tmp];
		cout << "Junc_file\t" << tmp-3 << "\t" << tmpJuncFile << endl; 
		juncFileVec.push_back(tmpJuncFile);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load each junc file" << endl;
	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < juncFileNum; tmp++)
	{
		cout << "start to load juncHash " << tmp + 1 << ": " << juncFileVec[tmp] << endl;
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionHashInfo(chromNum);
		string tmpJuncFile = juncFileVec[tmp];
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load junc files in vec to a merged junc hash" << endl;
	ofstream merged_ofs(outputMergedFile.c_str());

	AlignInferJunctionHash_Info* juncHash_merged = new AlignInferJunctionHash_Info();
	juncHash_merged->initiateAlignInferJunctionHashInfo(chromNum);
	juncHash_merged->insertJuncFromJuncFileVec_chrNamePos_supportNum(juncFileVec, indexInfo);
	int juncVecSize_merged = juncHash_merged->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < juncVecSize_merged; tmp++)
	{
		int tmpJunc_chrNameInt = juncHash_merged->returnAlignInferInfo_chrNameInt(tmp);
		string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
		int tmpJunc_donerEndPos = juncHash_merged->returnAlignInferInfo_donerEndPos(tmp);
		int tmpJunc_acceptorStartPos = juncHash_merged->returnAlignInferInfo_acceptorStartPos(tmp);
		merged_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_donerEndPos << "\t" << tmpJunc_acceptorStartPos << "\t";
		vector<int> tmpJunc_supNumVec_inEachFile;
		int tmpJunc_totalSupNum_inAllFiles = 0;
		for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncFileNum; tmpJuncHashIndex ++)
		{
			int tmpJunc_supNum_inTmpFile = juncHashVec[tmpJuncHashIndex]->searchAndReturnAlignInferJuncHashSupNum(
				tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
			tmpJunc_supNumVec_inEachFile.push_back(tmpJunc_supNum_inTmpFile);
			tmpJunc_totalSupNum_inAllFiles += tmpJunc_supNum_inTmpFile;
		}
		merged_ofs << tmpJunc_totalSupNum_inAllFiles;
		for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncFileNum; tmpJuncHashIndex ++)
			merged_ofs << "\t" << tmpJunc_supNumVec_inEachFile[tmpJuncHashIndex];
		merged_ofs << endl;
	}	
	delete juncHash_merged;
	for(int tmp = 0; tmp < juncFileNum; tmp++)
		delete juncHashVec[tmp];
	merged_ofs.close();
	free(chrom);
	delete indexInfo;
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}