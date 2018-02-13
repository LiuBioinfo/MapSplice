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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 7)
	{
		cout << "Executable inputIndexFolderPath SJsizeMin SJmaxMin outputFolder inputJuncInList2checkCoverage SJ_1 (SJ_2 ...)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "SJsizeMin: " << argv[2] << endl;
	log_ofs << "SJmaxMin: " << argv[3] << endl;
	log_ofs << "outputFolder: " << argv[4] << endl;
	log_ofs << "inputJuncInList2checkCoverage: " << argv[5] << endl;

	string juncCov = outputFolderStr + "juncExistence.txt";
	ofstream juncCov_ofs(juncCov.c_str());
	string invalidJuncFile = outputFolderStr + "invalidJunc.txt";
	ofstream invalidJunc_ofs(invalidJuncFile.c_str());
	vector<string> juncFileVec;
	for(int tmp = 6; tmp <= argc - 1; tmp++)
	{
		string tmpJuncFile = argv[tmp];
		juncFileVec.push_back(tmpJuncFile);
		log_ofs << "totalJuncFile: " << tmpJuncFile << endl;
	}
	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	cout << "assigning parameters ......" << endl;
	log_ofs << "assigning parameters ......" << endl;
	string min_intron_size_str = argv[2];
	int min_intron_size_int = atoi(min_intron_size_str.c_str());
	log_ofs << "min_intron_size_int: " << min_intron_size_int << endl;
	string max_intron_size_str = argv[3];
	int max_intron_size_int = atoi(max_intron_size_str.c_str());	
	log_ofs << "max_intron_size_int: " << max_intron_size_int << endl;

	cout << "start to initiate alignInferJunctionHashInfo_2cmp ...." << endl;
	log_ofs << "start to initiate alignInferJunctionHashInfo_2cmp ...." << endl;
	string juncFile_2cmp = argv[5];
	AlignInferJunctionHash_Info* juncHash_2cmp = new AlignInferJunctionHash_Info();
	juncHash_2cmp->initiateAlignInferJunctionInfo(chromNum);
	juncHash_2cmp->insertJuncFromJuncFile_chrNamePos_supportNum(juncFile_2cmp, indexInfo);
	cout << "start to initiate alignInferJunctionHashInfo_vec ...." << endl;
	log_ofs << "start to initiate alignInferJunctionHashInfo_vec ...." << endl;
	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < juncFileVec.size(); tmp++)
	{
		string tmpJuncFile = juncFileVec[tmp];
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);		
	}
	cout << "start to output ...." << endl;
	log_ofs << "start to output ...." << endl;
	int toCmpJuncNum = juncHash_2cmp->returnAlignInferInfoVecSize();
	int toCmpJuncNum_valid = 0;
	int toCmpJuncNum_invalid = 0;	
	for(int tmp = 0; tmp < toCmpJuncNum; tmp++)
	{
		int tmpJunc_chrNameInt = juncHash_2cmp->returnAlignInferInfo_chrNameInt(tmp);
		int tmpJunc_startPos = juncHash_2cmp->returnAlignInferInfo_donerEndPos(tmp);
		int tmpJunc_endPos = juncHash_2cmp->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpJunc_supNum = juncHash_2cmp->returnAlignInferInfo_supportNum(tmp);
		int tmpJunc_size = tmpJunc_endPos - tmpJunc_startPos - 1;
		if((tmpJunc_size > max_intron_size_int)||(tmpJunc_size < min_intron_size_int))
		{
			invalidJunc_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos << endl;
			toCmpJuncNum_invalid ++;
			continue;
		}
		toCmpJuncNum_valid ++;
		juncCov_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" << tmpJunc_supNum;
		for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncFileVec.size(); tmpJuncHashIndex++)
		{
			int tmpSupNum_in_tmpJuncHash = juncHashVec[tmpJuncHashIndex]->searchAndReturnAlignInferJuncHashSupNum(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
			juncCov_ofs << "\t" << tmpSupNum_in_tmpJuncHash;
		}
		juncCov_ofs << endl;
	}	
	log_ofs << "toCmpJuncNum: " << toCmpJuncNum << endl;
	log_ofs << "toCmpJuncNum_valid: " << toCmpJuncNum_valid << endl;
	log_ofs << "toCmpJuncNum_invalid: " << toCmpJuncNum_invalid << endl;
	for(int tmp = 0; tmp < juncFileVec.size(); tmp++)
		delete juncHashVec[tmp];
	delete juncHash_2cmp;
	juncCov_ofs.close();
	invalidJunc_ofs.close();
	log_ofs.close();
	return 0;
}