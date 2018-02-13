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
#include "general/paralogGene.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputDuplicatedGeneDatabaseFile outputFilePrefix" << endl;
		exit(1);
	}
	cout << "start to load duplicatedGeneDatabaseFile" << endl;
	string inputDuplicatedGeneDatabaseFile = argv[1];
	ParalogGeneGroupVec tmpParalogGeneGroupVecInfo;
	tmpParalogGeneGroupVecInfo.initiate_duplicatedGeneDatabaseFile(inputDuplicatedGeneDatabaseFile);
	cout << "start to output paralogGene files" << endl;
	string outputFilePrefix = argv[2];
	string outputFile_duplicatedGeneDatabaseFile = outputFilePrefix + ".duplicatedGeneDatabase";
	string outputFile_reformattedParalogGeneGroupFile = outputFilePrefix + ".reformattedParalogGeneGroup";
	tmpParalogGeneGroupVecInfo.outputDuplicatedGeneDatabaseFile(outputFile_duplicatedGeneDatabaseFile);
	tmpParalogGeneGroupVecInfo.outputReformattedParalogGeneGroupFile(outputFile_reformattedParalogGeneGroupFile);
	return 0;
}