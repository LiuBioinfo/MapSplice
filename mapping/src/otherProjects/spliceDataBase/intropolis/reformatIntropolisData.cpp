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

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder intropolisData outputFolder" << endl;
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

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	string reformatted_file = outputFolderStr + "intropolis_reformatted.txt";
	ofstream reformatted_ofs(reformatted_file.c_str());
	string invalid_file = outputFolderStr + "invalid.txt";
	ofstream invalid_ofs(invalid_file.c_str());
	string intropolisData_file = argv[2];
	ifstream intropolis_ifs(intropolisData_file.c_str());
	int tmpJuncNum = 0;
	int tmpJuncNum_valid = 0;
	int tmpJuncNum_invalid = 0;
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
		vector<string> tmpJuncFieldVec_ori;
		int tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			if(tmpTabLoc != string::npos)
			{	
				string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
				tmpJuncFieldVec_ori.push_back(tmpFieldStr);
				tmpStartLoc = tmpTabLoc + 1;
			}
			else
			{	
				string tmpLastFieldStr = tmpStr.substr(tmpStartLoc);
				tmpJuncFieldVec_ori.push_back(tmpLastFieldStr);
				break;
			}
		}		
		string tmpJunc_chrName_ori = tmpJuncFieldVec_ori[0];
		string tmpJunc_intronStartPosStr_ori = tmpJuncFieldVec_ori[1];
		int tmpJunc_intronStartPos_ori = atoi(tmpJunc_intronStartPosStr_ori.c_str());
		string tmpJunc_intronEndPosStr_ori = tmpJuncFieldVec_ori[2];
		int tmpJunc_intronEndPos_ori = atoi(tmpJunc_intronEndPosStr_ori.c_str());
		string tmpJunc_strand_ori = tmpJuncFieldVec_ori[3];
		//string tmpJunc_fkStr_1_ori = tmpJuncFieldVec_ori
		string tmpJunc_sampleIdStr_ori = tmpJuncFieldVec_ori[6];
		string tmpJunc_sampleCovStr_ori = tmpJuncFieldVec_ori[7];

		string tmpJunc_chrName = tmpJunc_chrName_ori;
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		if(tmpJunc_chrNameInt < 0)
		{
			invalid_ofs << tmpStr << endl;
			tmpJuncNum_invalid ++;
		}
		else
			tmpJuncNum_valid ++;
		int tmpJunc_startPos = tmpJunc_intronStartPos_ori - 1;
		int tmpJunc_endPos = tmpJunc_intronEndPos_ori + 1;
		string tmpJunc_flankString = indexInfo->returnFlankString(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
		string tmpJunc_strand = tmpJunc_strand_ori;
		string tmpJunc_sampleIdStr = tmpJunc_sampleIdStr_ori + ",";
		vector<string> tmpJunc_sampleIdVec;
		tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpCommaLoc = tmpJunc_sampleIdStr.find(",", tmpStartLoc);
			if(tmpCommaLoc != string::npos)
			{	
				string tmpFieldStr = tmpJunc_sampleIdStr.substr(tmpStartLoc, tmpCommaLoc - tmpStartLoc);
				tmpJunc_sampleIdVec.push_back(tmpFieldStr);
				tmpStartLoc = tmpCommaLoc + 1;
			}
			else
				break;
		}
		int tmpJunc_sampleNum = tmpJunc_sampleIdVec.size();
		string tmpJunc_sampleCovStr = tmpJunc_sampleCovStr_ori + ",";
		vector<int> tmpJunc_covVec;
		int tmpJunc_covSum = 0;
		tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpCommaLoc = tmpJunc_sampleCovStr.find(",", tmpStartLoc);
			if(tmpCommaLoc != string::npos)
			{	
				string tmpFieldStr = tmpJunc_sampleCovStr.substr(tmpStartLoc, tmpCommaLoc - tmpStartLoc);
				int tmpJuncCovInTmpSample = atoi(tmpFieldStr.c_str());
				tmpJunc_covVec.push_back(tmpJuncCovInTmpSample);
				tmpJunc_covSum += tmpJuncCovInTmpSample;
				tmpStartLoc = tmpCommaLoc + 1;
			}
			else
				break;
		}
		reformatted_ofs << tmpJunc_chrName << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" 
			<< tmpJunc_strand << "\t" << tmpJunc_flankString << "\t" << tmpJunc_sampleNum << "\t" 
			<< tmpJunc_covSum << "\t" << tmpJunc_sampleIdStr << "\t" << tmpJunc_sampleCovStr << endl;
	}
	log_ofs << "JuncNum: " << tmpJuncNum << endl;
	log_ofs << "JuncNum_valid: " << tmpJuncNum_valid << endl;
	log_ofs << "JuncNum_invalid: " << tmpJuncNum_invalid << endl;
	delete indexInfo;
	free(chrom);
	intropolis_ifs.close();
	reformatted_ofs.close();
	invalid_ofs.close();
	log_ofs.close();
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}