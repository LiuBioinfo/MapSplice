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
	if(argc != 6)
	{
		cout << "Executable inputIndex inputFaFolder outputBfBvFolder inputHashFile outputFile" << endl;
		exit(1);
	}	
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
	string inputFaFolder = argv[2];
	string outputBfBvFolder = argv[3];
	string inputHashFile = argv[4];
	string outputScript = argv[5];
	ofstream script_ofs(outputScript.c_str());
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		int tmpChr_subRegion_num = indexInfo->returnSecondLevelIndexPartsNum(tmpChr);
		cout << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;
		for(int tmpSubRegion = 0; tmpSubRegion < tmpChr_subRegion_num; tmpSubRegion ++)
		{		
			string tmpLocalRegionSubSeqReadFaFile = inputFaFolder + "/" 
				+ indexInfo->returnChrNameStr(tmpChr) + "_" + int_to_str(tmpSubRegion + 1) + ".fa";
			string tmpOutputBfBvFile = outputBfBvFolder + "/" + indexInfo->returnChrNameStr(tmpChr) 
				+ "_" + int_to_str(tmpSubRegion + 1) + ".bf.bv";
			string tmp_cmd = "./bt count --cutoff 1 --threads 16 " + inputHashFile + " 3000000000 "
				+ tmpLocalRegionSubSeqReadFaFile + " " + tmpOutputBfBvFile;
			script_ofs << tmp_cmd << endl;
		}
	}
	script_ofs.close();
	free(chrom);
	delete indexInfo;
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	return 0;
}