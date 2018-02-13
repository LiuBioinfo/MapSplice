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
#include <sys/types.h>    
#include <dirent.h>    
#include <stdio.h>    
#include <errno.h>
#include <set>

#include "build_snpMerIndex.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputMergedFaFile outputIndexFolder" << endl;
		exit(1);
	}

	string inputMergedFaFile = argv[1];
	string outputIndexFolder = argv[2];

	Build_snpMerIndex buildSnpMerIndexInfo;
	//cout << "MAX: " << MAX << endl;
	bool snpMerIndexBuild_success_bool = buildSnpMerIndexInfo.build_snpMerIndex(inputMergedFaFile, outputIndexFolder);
	//cout << "MAX: " << MAX << endl;
	string outputResultsFile = outputIndexFolder + "/results.txt";
	ofstream results_ofs(outputResultsFile.c_str());
	if(snpMerIndexBuild_success_bool)
		results_ofs << "success" << endl;
	else
		results_ofs << "fail" << endl;
	results_ofs.close();
	return 0;
}