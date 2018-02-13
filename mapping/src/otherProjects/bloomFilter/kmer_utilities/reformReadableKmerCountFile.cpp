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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputReadableKmerCountFile_ori outputReformedKmerCountFilePrefix" << endl;
		exit(1);
	}
	string inputReadableKmerCountFile_ori = argv[1];
	string outputReformedKmerCountFilePrefix = argv[2];
	string outputReformedKmerCountFilePrefix_seqCount = outputReformedKmerCountFilePrefix + ".seqCount";
	string outputReformedKmerCountFilePrefix_stats = outputReformedKmerCountFilePrefix + ".stats";
	string extractSeqCount_cmd = "head -n -2 " + inputReadableKmerCountFile_ori + " > "
		+ outputReformedKmerCountFilePrefix_seqCount;
	string extractStats_cmd = "tail -2 " + inputReadableKmerCountFile_ori + " > "
		+ outputReformedKmerCountFilePrefix_stats;
	system(extractSeqCount_cmd.c_str());
	system(extractStats_cmd.c_str());
	return 0;
}