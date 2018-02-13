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
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputAnnotatedJuncFile inputIntropolisJuncFile outputFolder" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
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

	string inputAnnotatedJuncFile = argv[2];
	AlignInferJunctionHash_Info* annotation_SJhash = new AlignInferJunctionHash_Info();
	annotation_SJhash->initiateAlignInferJunctionHashInfo(chromNum);
	annotation_SJhash->insertJuncFromJuncFile_chrNamePosOnly(inputAnnotatedJuncFile, indexInfo);
	
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	string annotated_intropolis_file = outputFolderStr + "annotated_junc.txt";
	string unannotated_intropolis_file = outputFolderStr + "unannotated_junc.txt";
	ofstream annotated_intropolis_ofs(annotated_intropolis_file.c_str());
	ofstream unannotated_intropolis_ofs(unannotated_intropolis_file.c_str());

	int juncNum_in_annotation = annotation_SJhash->returnAlignInferInfoVecSize();
	log_ofs << "juncNum_in_annotation: " << juncNum_in_annotation << endl;

	string inputIntropolisJuncFile = argv[3];
	ifstream intropolis_ifs(inputIntropolisJuncFile.c_str());
	int tmpJuncNum = 0;
	int tmpJuncNum_annotated = 0;
	int tmpJuncNum_unannotated = 0;
	while(!intropolis_ifs.eof())
	{
		string tmpStr;
		getline(intropolis_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpJuncNum ++;
		int tmpThousandIndex = tmpJuncNum / 100000;
		if(tmpJuncNum == tmpThousandIndex * 100000)
			cout << "Processed Junc #: " << tmpJuncNum << endl;
		vector<string> tmpJuncFieldVec;
		int tmpStartLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
			tmpJuncFieldVec.push_back(tmpFieldStr);
			tmpStartLoc = tmpTabLoc + 1;
		}
		string tmpJunc_chrName = tmpJuncFieldVec[0];
		string tmpJunc_startPosStr = tmpJuncFieldVec[1];
		string tmpJunc_endPosStr = tmpJuncFieldVec[2];
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		bool tmpJunc_annotatedOrNot_bool = annotation_SJhash->SJexistInAlignInferJuncHash(
			tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
		if(tmpJunc_annotatedOrNot_bool)
		{
			tmpJuncNum_annotated ++;
			annotated_intropolis_ofs << tmpStr << endl;
		}
		else
		{
			tmpJuncNum_unannotated ++;
			unannotated_intropolis_ofs << tmpStr << endl;
		}
	}
	intropolis_ifs.close();

	delete indexInfo;
	delete annotation_SJhash;
	annotated_intropolis_ofs.close();
	unannotated_intropolis_ofs.close();
	free(chrom);
	intropolis_ifs.close();
	log_ofs.close();
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}	