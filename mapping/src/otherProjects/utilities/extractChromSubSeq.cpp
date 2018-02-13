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

#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath chrName startPos endPos" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	//cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	//cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string chrName = argv[2];
	int chrNameInt = indexInfo->convertStringToInt(chrName);
	if(chrNameInt < 0)
	{
		cout << "invalid chrName: " << chrName << endl;
		exit(1);
	}
	int chrSeqLength = indexInfo->returnChromLength(chrNameInt);
	string startPosStr = argv[3];
	int startPos = atoi(startPosStr.c_str());
	string endPosStr = argv[4];
	int endPos = atoi(endPosStr.c_str());
	if(startPos < 1)
	{
		cout << "error! invalid startPos: " << startPos << endl;
		exit(1);
	}
	if(endPos > chrSeqLength)
	{
		cout << "error! endPos > chrSeqLength, endPos: " << endPos << " chrSeqLength: " << chrSeqLength << endl;
		exit(1);
	}
	if(startPos > endPos)
	{
		cout << "error! startPos > endPos" << endl;
		exit(1);
	}
	cout << chrName << ": " << startPos << " ~ " << endPos << " " << indexInfo->returnChromStrSubstr(
		chrNameInt, startPos, endPos - startPos + 1) << endl;
	return 0;
}