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
#include <sstream>
#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../../incorporateGenomicVariants/general/SNPhash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputIndexFolderPath" << endl;
		cout << "#2 inputSNPfile_ann" << endl;
		cout << "#3 outputDir" << endl;
		//cout << "#4 SNPmerLength" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputSNPfile_ann = argv[2];
	string outputFolderStr = argv[3];
	// string SNPmerLengthStr = argv[4];
	// int SNPmerLength = atoi(SNPmerLengthStr.c_str());
	vector<int> snpMerLengthVec;
	snpMerLengthVec.push_back(201);
	snpMerLengthVec.push_back(101);
	snpMerLengthVec.push_back(51);
	snpMerLengthVec.push_back(41);
	snpMerLengthVec.push_back(31);
	snpMerLengthVec.push_back(21);

	vector< vector<int> > otherSNPnumVecVec;
	for(int tmp = 0; tmp < snpMerLengthVec.size(); tmp++)
	{
		vector<int> tmpSNPmerLength_otherSNPnumVec;
		int tmpSNPmerLength = snpMerLengthVec[tmp];
		int tmpSNPmerLength_otherSNPnum_max = tmpSNPmerLength - 1;
		for(int tmpSNPnum = 0; tmpSNPnum <= tmpSNPmerLength_otherSNPnum_max; tmpSNPnum++)
			tmpSNPmerLength_otherSNPnumVec.push_back(0);
		otherSNPnumVecVec.push_back(tmpSNPmerLength_otherSNPnumVec);
	}	

	outputFolderStr += "/";
	string cmd_mkdir_outputDir = "mkdir -p " + outputFolderStr;
	system(cmd_mkdir_outputDir.c_str());

	string outputFile_log = outputFolderStr + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	cout << "Command: \n" << argv[0] << endl << argv[1] << endl << argv[2] << endl << argv[3] << endl;
	log_ofs << "Command: \n" << argv[0] << endl << argv[1] << endl << argv[2] << endl << argv[3] << endl;

	cout << "initiate indexInfo ..." << endl;
	log_ofs << "initiate indexInfo ..." << endl;
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	chrom_bit_file_ifs.close();
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	parameter_ifs.close();
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;

	cout << "initiate snpHashInfo ..." << endl;
	log_ofs << "initiate snpHashInfo ..." << endl;
	SNPhash_Info snpHashInfo;
	snpHashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	snpHashInfo.generateSNPhash_formattedSNPfile(inputSNPfile_ann, indexInfo);
	cout << "end of initiating snpHashInfo ..." << endl;
	log_ofs << "end of initiating snpHashInfo ..." << endl;

	int SNP_num = snpHashInfo.returnSNPnum();
	cout << "SNP_num: " << SNP_num << endl;
	log_ofs << "SNP_num: " << SNP_num << endl;
	for(int tmp = 0; tmp < SNP_num; tmp++)
	{
		int tmpSNP_chrNameInt = snpHashInfo.returnSNP_chrNameInt(tmp);
		int tmpSNP_chrPos = snpHashInfo.returnSNP_chrPos(tmp);
		for(int tmpSNPmerIndex = 0; tmpSNPmerIndex < snpMerLengthVec.size(); tmpSNPmerIndex ++)
		{
			int tmpSNPmerLength = snpMerLengthVec[tmpSNPmerIndex];
			int halfLength = (tmpSNPmerLength - 1)/2;
			int tmpRegion_startPos = tmpSNP_chrPos - halfLength;
			int tmpRegion_endPos = tmpSNP_chrPos + halfLength;
			vector<int> tmpSNPposVec;
			snpHashInfo.returnSNPposVecWithinRegion(tmpSNP_chrNameInt, tmpRegion_startPos, 
				tmpRegion_endPos, tmpSNPposVec);
			int tmpOtherSNPnum = tmpSNPposVec.size() - 1;
			(otherSNPnumVecVec[tmpSNPmerIndex])[tmpOtherSNPnum] ++;
		}
	}

	cout << "start to print snp # distribution!" << endl;
	log_ofs << "start to print snp # distribution!" << endl;
	for(int tmpSNPmerIndex = 0; tmpSNPmerIndex < snpMerLengthVec.size(); tmpSNPmerIndex ++)
	{
		int tmpSNPmerLength = snpMerLengthVec[tmpSNPmerIndex];
		string output_snpMer_snpNum_distribution = outputFolderStr + int_to_str(tmpSNPmerLength) + "mer.otherSNPnum.distribution.txt";
		ofstream snpMer_snpNum_ofs(output_snpMer_snpNum_distribution.c_str());
		for(int tmpSNPnum = tmpSNPmerLength - 1; tmpSNPnum >= 0; tmpSNPnum --)
			snpMer_snpNum_ofs << tmpSNPnum << "\t" << (otherSNPnumVecVec[tmpSNPmerIndex])[tmpSNPnum] << endl;
		snpMer_snpNum_ofs.close();
	}

	free(chrom);
	delete indexInfo;
	log_ofs.close();
	return 0;
}