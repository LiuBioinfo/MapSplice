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
#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSailfishResults outputFilePrefix" << endl;
		exit(1);
	}
	string inputSailfishResults = argv[1];
	string outputFilePrefix = argv[2];
	string outputFile_raw = outputFilePrefix + ".raw";
	string outputFile_rounded = outputFilePrefix + ".rounded";

	ifstream sf_ifs(inputSailfishResults.c_str());
	ofstream readCount_ofs_raw(outputFile_raw.c_str());
	ofstream readCount_ofs_rounded(outputFile_rounded.c_str());
	string str_1;
	getline(sf_ifs, str_1); // 1st line
	while(!sf_ifs.eof())
	{
		string tmpSfStr;
		getline(sf_ifs, tmpSfStr);
		if(tmpSfStr == "")
			break;
		int tabLoc_1 = tmpSfStr.find("\t");
		int tabLoc_2 = tmpSfStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpSfStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpSfStr.find("\t", tabLoc_3 + 1);
		string tmpTranscriptFullName = tmpSfStr.substr(0, tabLoc_1);
		int straightLineLoc = tmpTranscriptFullName.find("|");
		string tmpTranscriptName = tmpTranscriptFullName.substr(0, straightLineLoc);
		string tmpReadCountStr = tmpSfStr.substr(tabLoc_4 + 1);
		readCount_ofs_raw << tmpTranscriptName << "\t" << tmpReadCountStr << endl;
		double tmpReadCount_double = atof(tmpReadCountStr.c_str());
		int tmpReadCount_int = (int)(tmpReadCount_double + 0.5);
		readCount_ofs_rounded << tmpTranscriptName << "\t" << tmpReadCount_int << endl;
	}
	sf_ifs.close();
	readCount_ofs_raw.close();
	readCount_ofs_rounded.close();
	return 0;
}