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
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSam" << endl;
		cout << "#2 outputFilePrefix" << endl;
		exit(1);
	}

	string inputSam = argv[1];
	string outputFilePrefix = argv[2];
	string outputUnmapSamFile = outputFilePrefix + ".unmap.sam";
	string outputReadFile_end1 = outputFilePrefix + ".unmap.1.fastq";
	string outputReadFile_end2 = outputFilePrefix + ".unmap.2.fastq";
	ofstream unmapSam_ofs(outputUnmapSamFile.c_str());
	ofstream read_ofs_1(outputReadFile_end1.c_str());
	ofstream read_ofs_2(outputReadFile_end2.c_str());
	ifstream sam_ifs(inputSam.c_str());
	int tmpLineNO = 0;
	while(!sam_ifs.eof())
	{
		string tmpSamStr_1;
		getline(sam_ifs, tmpSamStr_1);
		if(tmpSamStr_1 == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 1000000;
		if(tmpLineNO == tmpThousandIndex * 1000000)
			cout << "Processed Line #: " << tmpLineNO << endl;
		if(tmpSamStr_1.at(0) == '@')
			continue;
		string tmpSamStr_2;
		getline(sam_ifs, tmpSamStr_2);
		
		vector<string> samFieldVec_1;
		vector<string> samFieldVec_2;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpSamStr_1.find("\t", startLoc);
			string tmpSamField = tmpSamStr_1.substr(startLoc, tabLoc-startLoc);
			samFieldVec_1.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpSamStr_2.find("\t", startLoc);
			string tmpSamField = tmpSamStr_2.substr(startLoc, tabLoc-startLoc);
			samFieldVec_2.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		if(!((samFieldVec_1[2] == "*")&&(samFieldVec_2[2] == "*")))
			continue;
		string tmpReadId_1 = samFieldVec_1[0];
		string tmpReadSeq_1 = samFieldVec_1[9];
		string tmpReadQual_1 = samFieldVec_1[10];
		string tmpReadId_2 = samFieldVec_2[0];
		string tmpReadSeq_2 = samFieldVec_2[9];
		string tmpReadQual_2 = samFieldVec_2[10];
		unmapSam_ofs << tmpSamStr_1 << endl << tmpSamStr_2 << endl;
		read_ofs_1 << "@" << tmpReadId_1 << "/1" << endl << tmpReadSeq_1 << endl << "+" << endl << tmpReadQual_1 << endl;
		read_ofs_2 << "@" << tmpReadId_2 << "/2" << endl << tmpReadSeq_2 << endl << "+" << endl << tmpReadQual_2 << endl;
	}
	unmapSam_ofs.close();
	sam_ifs.close();
	read_ofs_2.close();
	read_ofs_1.close();
	return 0;
}