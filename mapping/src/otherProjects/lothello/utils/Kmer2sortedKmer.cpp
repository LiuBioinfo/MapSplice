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
#include "../../../general/read_block_test.h"
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"
#include "../general/Kmer2sortedKmer_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable temporaryDir threads_num inputRawKmerFile outputSortedKmerFile" << endl;
		exit(1);
	}
	string temporaryDir = argv[1];
	temporaryDir += "/";
	string threads_num_str = argv[2];
	int threads_num = atoi(threads_num_str.c_str());
	string inputRawKmerFile = argv[3];
	string outputSortedKmerFile = argv[4];

	Kmer2sortedKmer_Info kmer2sortedKmerInfo;
	kmer2sortedKmerInfo.initiate_mkTmpDir(inputRawKmerFile, outputSortedKmerFile, temporaryDir);
	kmer2sortedKmerInfo.sort_Kmer_keepInterDir(threads_num);
	return 0;
}