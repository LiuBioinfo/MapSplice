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
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputFa" << endl;
		cout << "#2 outputSeq" << endl;
		cout << "#3 minLength" << endl;
		exit(1);
	}
	string inputFa = argv[1];
	string outputSeq = argv[2];
	string minLengthStr = argv[3];
	int minLength = atoi(minLengthStr.c_str());
	ifstream fa_ifs(inputFa.c_str());
	ofstream seq_ofs(outputSeq.c_str());
	while(!fa_ifs.eof())
	{
		string tmpId;
		getline(fa_ifs, tmpId);
		if(tmpId == "")
			break;
		string tmpSeq;
		getline(fa_ifs, tmpSeq);
		if(tmpSeq.length() >= minLength)
			seq_ofs << tmpSeq << endl;
	}
	fa_ifs.close();
	seq_ofs.close();
	return 0;
}