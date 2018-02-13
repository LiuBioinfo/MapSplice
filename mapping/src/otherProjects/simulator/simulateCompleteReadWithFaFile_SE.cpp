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
#include "../../general/index_info.h"

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputFaFile outputFile readLength" << endl;
		exit(1);
	}
	string inputFaFile = argv[1];
	string outputFile = argv[2];
	string readLengthStr = argv[3];
	int readLength = atoi(readLengthStr.c_str());
	ifstream fa_ifs(inputFaFile.c_str());
	ofstream read_ofs(outputFile.c_str());
	while(!fa_ifs.eof())
	{
		string tmpNameStr, tmpSeqStr;
		getline(fa_ifs, tmpNameStr);
		if(tmpNameStr == "")
			break;
		getline(fa_ifs, tmpSeqStr);
		int tmpSeqLength = tmpSeqStr.length();
		if(tmpSeqLength < readLength)
			continue;
		int tmpStartPos_simulatedReadStartLocInSeq = 1;
		int tmpEndPos_simulatedReadStartLocInSeq = tmpSeqLength - readLength + 1;
		for(int tmp = tmpStartPos_simulatedReadStartLocInSeq; tmp <= tmpEndPos_simulatedReadStartLocInSeq; tmp++)
		{
			string tmpSimulatedRead_name = tmpNameStr + "_" + int_to_str(tmp);
			string tmpSimulatedRead_seq = tmpSeqStr.substr(tmp - 1, readLength);
			read_ofs << tmpSimulatedRead_name << endl << tmpSimulatedRead_seq << endl;
		}
	}
	read_ofs.close();
	fa_ifs.close();
	return 0;
}

