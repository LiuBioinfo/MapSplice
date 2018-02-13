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
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 index_dir" << endl;
		cout << "#2 exon_gtf_file" << endl;
		cout << "#3 output_file" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string exon_gtf_file = argv[2];
	string output_file = argv[3];

	cout << "start to initiate indexInfo" << endl;
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

	int tmpSeqNO = 0;
	ifstream exonGtf_ifs(exon_gtf_file.c_str());
	ofstream exonSeq_ofs(output_file.c_str());
	while(!exonGtf_ifs.eof())
	{
		string tmpStr;
		getline(exonGtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1+1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2+1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3+1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4+1);
		int tabLoc_6 = tmpStr.find("\t", tabLoc_5+1);
		int tabLoc_7 = tmpStr.find("\t", tabLoc_6+1);
		string tmpChrName = tmpStr.substr(0, tabLoc_1);
		string tmpFeatureType = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpStartPosStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		string tmpEndPosStr = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpStrand = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
	
		int tmpChrNameIndex = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameIndex >= 0)
		{
			tmpSeqNO ++;
			string tmpId = ">" + tmpChrName + "_" + tmpStartPosStr + "_" 
				+ tmpEndPosStr + "_";
			string tmpSeq = indexInfo->returnChromStrSubstr(tmpChrNameIndex, tmpStartPos, tmpEndPos - tmpStartPos + 1);
			exonSeq_ofs << tmpId << tmpSeqNO << endl << tmpSeq << endl;
		}
	}
	exonGtf_ifs.close();
	exonSeq_ofs.close();
	delete indexInfo;
	free(chrom);
	return 0;
}
