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
#include "general/alu_total_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{	
		cout << "Executable inputIndexFolderPath inputAluInfoFile inputAluIndexMapSamFile outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "initiate indexInfo ..." << endl;
	log_ofs << endl << "[" << asctime(local) << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "end of initiating indexInfo; start to do aluInfo initiating" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of initiating indexInfo; start to do aluInfo initiating" << endl;

	string inputAluInfoFile = argv[2];
	string outputInvalidAluInfoFile = outputFolderStr + "invalid_aluInfo.txt";
	Alu_Total_Info tmpAluTotalInfo;
	tmpAluTotalInfo.initiate_withAluInfoFile(inputAluInfoFile, indexInfo, outputInvalidAluInfoFile);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "end of aluInfo initiating; start to do alu element quantification" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of aluInfo initiating; start to do alu element quantification" << endl;

	string inputAluIndexMapSamFile = argv[3];
	string outputInvalidSamFile = outputFolderStr + "invalid.sam";
	string readAssignmentFile = outputFolderStr + "read_assignment.txt";
	string aluElementQuantificationFile = outputFolderStr + "aluElement_quant.txt";
	tmpAluTotalInfo.estimateAluEleAbundance(inputAluIndexMapSamFile, outputInvalidSamFile, 
		readAssignmentFile, aluElementQuantificationFile, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "end of alu element quantification; All jobs done !" << endl;
	log_ofs << endl << "[" << asctime(local) << "end of alu element quantification; All jobs done !" << endl;

	log_ofs.close();
	free(chrom);
	delete indexInfo;	
	return 0;
}