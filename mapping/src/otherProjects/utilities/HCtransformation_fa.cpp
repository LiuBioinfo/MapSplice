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
		cout << "Executable inputFa outputHCcompressedFilePrefix" << endl;
		exit(1);
	}
	string inputFa = argv[1];
	string outputHCcompressedFilePrefix = argv[2];
	string outputHCcompressedFile_fa = outputHCcompressedFilePrefix + ".HC.fa";
	string outputHCcompressedFile_pos = outputHCcompressedFilePrefix + ".HC.pos";
	string outputHCcompressedFile_len = outputHCcompressedFilePrefix + ".HC.len";
	ifstream fa_ifs(inputFa.c_str());
	ofstream HCfa_ofs(outputHCcompressedFile_fa.c_str());
	ofstream HCpos_ofs(outputHCcompressedFile_pos.c_str());
	ofstream HClen_ofs(outputHCcompressedFile_len.c_str());
	int tmpHClen = 0;
	while(!fa_ifs.eof())
	{
		string tmpId;
		getline(fa_ifs, tmpId);
		if(tmpId == "")
			break;
		string tmpSeq;
		getline(fa_ifs, tmpSeq);
		int tmpSeqLen = tmpSeq.length();
		char lastBaseChar = 'X';
		tmpHClen = 0;
		HCfa_ofs << tmpId << endl;
		for(int tmp = 0; tmp < tmpSeqLen; tmp++)
		{
			char tmpChar = tmpSeq.at(tmp);
			if(tmpChar != lastBaseChar)
			{
				tmpHClen ++;
				HCfa_ofs << tmpChar;
				HCpos_ofs << tmp << ",";
				lastBaseChar = tmpChar;
			}
		}
		HClen_ofs << tmpHClen << endl;
		HCfa_ofs << endl;
		HCpos_ofs << endl;
	}
	fa_ifs.close();
	HCfa_ofs.close();
	HCpos_ofs.close();
	return 0;
}