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
#define CHR_SEQ_FA_LINE_LENGTH 50
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputRawCab3File" << endl;
		cout << "#2 outputReformedCab3ChrSeqDir" << endl;
		exit(1);
	}

	string inputRawCab3File = argv[1];
	string outputReformedCab3ChrSeqDir = argv[2];

	cout << "creating folder ......" << endl;
	string outputFolderStr = outputReformedCab3ChrSeqDir + "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	

	ifstream rawCab3_ifs(inputRawCab3File.c_str());
	while(!rawCab3_ifs.eof())
	{
		string tmpId;
		getline(rawCab3_ifs, tmpId);
		if(tmpId == "")
			break;
		string tmpReformedChrSeqFile = outputFolderStr + tmpId.substr(1) + ".fa";
		ofstream reformedChrSeq_ofs(tmpReformedChrSeqFile.c_str());
		string tmpSeq;
		getline(rawCab3_ifs, tmpSeq);
		int tmpSeqLength = tmpSeq.length();
		int tmpFullBaseLineNum = tmpSeqLength / CHR_SEQ_FA_LINE_LENGTH;
		bool lastLineFullBaseOrNot = (tmpFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH == tmpSeqLength);
		// print id
		reformedChrSeq_ofs << tmpId << endl;
		// print all fullBase lines
		for(int tmp = 0; tmp < tmpFullBaseLineNum; tmp++)
			reformedChrSeq_ofs << tmpSeq.substr(tmp * CHR_SEQ_FA_LINE_LENGTH, 50) << endl;
		// print the last line
		if(!lastLineFullBaseOrNot)
			reformedChrSeq_ofs << tmpSeq.substr(tmpFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH, 
				tmpSeqLength - tmpFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH) << endl;
		reformedChrSeq_ofs.close();
	}
	rawCab3_ifs.close();
	return 0;
}