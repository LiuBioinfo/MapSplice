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
		cout << "#1 inputRawCab3ChrUnFile" << endl;
		cout << "#2 outputReformedCab3ChrUnFile" << endl;
		exit(1);
	}

	string inputRawCab3ChrUnFile = argv[1];
	string outputReformedCab3ChrUnFile = argv[2];

	string chrUnSeq = "";
	ifstream rawCab3ChrUn_ifs(inputRawCab3ChrUnFile.c_str());
	while(!rawCab3ChrUn_ifs.eof())
	{
		string tmpId;
		getline(rawCab3ChrUn_ifs, tmpId);
		if(tmpId == "")
			break;
		string tmpSeq;
		getline(rawCab3ChrUn_ifs, tmpSeq);
		chrUnSeq += tmpSeq;
	}
	rawCab3ChrUn_ifs.close();

	ofstream reformedCab3ChrUn_ofs(outputReformedCab3ChrUnFile.c_str());
	reformedCab3ChrUn_ofs << ">chrUn" << endl;
	int chrUnSeqLength = chrUnSeq.length();
	int chrUnFullBaseLineNum = chrUnSeqLength / CHR_SEQ_FA_LINE_LENGTH;
	bool lastLineFullBaseOrNot = (chrUnFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH == chrUnSeqLength);
	// print all fullBase lines
	for(int tmp = 0; tmp < chrUnFullBaseLineNum; tmp++)
		reformedCab3ChrUn_ofs << chrUnSeq.substr(tmp * CHR_SEQ_FA_LINE_LENGTH, 50) << endl;
	// print the last line
	if(!lastLineFullBaseOrNot)
		reformedCab3ChrUn_ofs << chrUnSeq.substr(chrUnFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH,
			chrUnSeqLength - chrUnFullBaseLineNum * CHR_SEQ_FA_LINE_LENGTH) << endl;
	reformedCab3ChrUn_ofs.close();
	return 0;
}