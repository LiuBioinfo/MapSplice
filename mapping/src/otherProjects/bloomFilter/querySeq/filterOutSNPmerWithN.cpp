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
		cout << "Executable inputKmerFile_gt4 outputFile" << endl;
		exit(1);
	}
	string inputKmerFile_gt4 = argv[1];
	string outputFile = argv[2];
	string outputFile_filteredOut = outputFile + ".filterOut_N";
	ifstream gt4_ifs(inputKmerFile_gt4.c_str());
	ofstream kept_ofs(outputFile.c_str());
	ofstream filteredOut_ofs(outputFile_filteredOut.c_str());
	while(!gt4_ifs.eof())
	{
		string tmpStr;
		getline(gt4_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string kmerStr = tmpStr.substr(0, tabLoc);
		if(kmerStr.find_first_not_of("ATCG") == string::npos)
			kept_ofs << tmpStr << endl;
		else
			filteredOut_ofs << tmpStr << endl;
	}
	gt4_ifs.close();
	kept_ofs.close();
	filteredOut_ofs.close();
	return 0;
}