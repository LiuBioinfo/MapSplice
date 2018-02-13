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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexInfoDir" << endl;
		cout << "#2 SNP_file" << endl;		
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	int chr_seq_fa_line_length = 50;
	string indexFolderPath = argv[1];
	string SNP_file = argv[2];
	string outputDir = argv[3];
	indexFolderPath += "/";
	cout << "start to initiate indexInfo" << endl;
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
	free(chrom);	

	indexInfo->maskSNPwithNinChromStr(SNP_file);

	cout << "start to output snp-masked chromStr" << endl;
	indexInfo->printEachChromStr2oneFile(outputDir, chr_seq_fa_line_length);
	cout << "All jobs done!" << endl;
	return 0;
}