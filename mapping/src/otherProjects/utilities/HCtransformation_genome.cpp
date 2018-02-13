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

//#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputIndex outputHCcompressedFileFolder" << endl;
		exit(1);
	}
	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	string outputFolderStr_seq = outputFolderStr + "seq/";
	string outputFolderStr_pos = outputFolderStr + "pos/";
	string mkdir_seq = "mkdir -p " + outputFolderStr_seq;
	string mkdir_pos = "mkdir -p " + outputFolderStr_pos;
	system(mkdir_seq.c_str());
	system(mkdir_pos.c_str());

	log_ofs << "start to initiate indexInfo" << endl;
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;	
	int tmpHClen = 0;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		cout << "start to do HC transformation for tmpChrName" << endl;
		log_ofs << "start to do HC transformation for tmpChrName" << endl;	
		string tmpChrSeq_HC = outputFolderStr_seq + tmpChrName + ".fa";
		string tmpChrPos_HC = outputFolderStr_pos + tmpChrName + ".pos";
		ofstream tmpChrSeq_HC_ofs(tmpChrSeq_HC.c_str());
		ofstream tmpChrPos_HC_ofs(tmpChrPos_HC.c_str());
		int tmpChrLength = indexInfo->returnChromLength(tmpChr);
		char lastBaseChar = 'M';
		tmpHClen = 0;
		tmpChrSeq_HC_ofs << ">" << tmpChrName << endl;
		for(int tmpBase = 0; tmpBase < tmpChrLength; tmpBase ++)
		{
			char tmpChar = indexInfo->returnOneBaseCharInGenome(tmpChr, tmpBase + 1);
			if(tmpChar != lastBaseChar)
			{
				tmpHClen ++;
				if((tmpHClen/50) * 50 == tmpHClen)
					tmpChrSeq_HC_ofs << tmpChar << endl;
				else
					tmpChrSeq_HC_ofs << tmpChar;				
				tmpChrPos_HC_ofs << tmpBase << ",";
				lastBaseChar = tmpChar;
			}
		}
		tmpChrSeq_HC_ofs.close();
		tmpChrPos_HC_ofs.close();
		log_ofs << "chrName: " << tmpChrName << " -- after HC, Length: " << tmpHClen << endl;
	}

	log_ofs.close();
	delete indexInfo;
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	free(chrom);
	log_ofs.close();
	return 0;
}