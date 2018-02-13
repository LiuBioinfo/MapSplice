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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSeqFile" << endl;
		cout << "#2 outputLengthFile" << endl;
		exit(1);
	}
	string inputSeqFile = argv[1];
	string outputLengthFile = argv[2];
	ifstream seq_ifs(inputSeqFile.c_str());
	ofstream len_ofs(outputLengthFile.c_str());
	while(!seq_ifs.eof())
	{
		string tmpStr;
		getline(seq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int seqLength = tmpStr.length();
		len_ofs << seqLength << endl;
	}
	seq_ifs.close();
	len_ofs.close();
	return 0;
}