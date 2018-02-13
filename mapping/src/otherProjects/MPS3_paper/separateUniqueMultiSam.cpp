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
//#include <omp.h>
// #include "read_block_test.h"
// #include "bwtmap_info.h"
// #include "DoubleAnchorScore.h"
// #include "sbndm.h"
// #include "splice_info.h"
// #include "jumpCode_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 3)
	{
		cout << "Executable <InputSam> <OutputFolder>" << endl;
		exit(0);
	}
	
	string inputSamFileStr = argv[1];
	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	//string outputLogFile = outputFolderStr + "log.txt";
	//ofstream log_ofs(outputLogFile.c_str());

	string outputSamFileStr_unique = outputFolderStr + "unique.sam";
	string outputSamFileStr_multi = outputFolderStr + "multi.sam";
	string outputSamFileStr_unmap = outputFolderStr + "unmap.sam";

	ofstream outputSamFile_unique_ofs(outputSamFileStr_unique.c_str());
	ofstream outputSamFile_multi_ofs(outputSamFileStr_multi.c_str());
	ofstream outputSamFile_unmap_ofs(outputSamFileStr_unmap.c_str());
	
	ifstream sam_ifs(inputSamFileStr.c_str());
	
	while(!sam_ifs.eof())
	{   
		string tmpStr;
		getline(sam_ifs, tmpStr);
		if(tmpStr.substr(0,1) == "@")
			continue;
		int tmpCandiAlignmentNum = 0;
		int NH_loc = tmpStr.find("NH:i:");
		int IH_loc = tmpStr.find("IH:i:");
		string tmpCandiAlignmentNumStr;
		if((NH_loc != string::npos)&&(IH_loc != string::npos))
		{
			cout << "error: NH_loc: " << NH_loc << " IH_loc: " << IH_loc << endl;
			exit(1);
		}
		else if(NH_loc != string::npos)
		{
			int nextTabLoc = tmpStr.find("\t", NH_loc + 1);
			if(nextTabLoc != string::npos)
				tmpCandiAlignmentNumStr = tmpStr.substr(NH_loc + 5, nextTabLoc - 1 - NH_loc - 5 + 1);
			else
				tmpCandiAlignmentNumStr = tmpStr.substr(NH_loc + 5);
		}
		else if(IH_loc != string::npos)
		{
			int nextTabLoc = tmpStr.find("\t", IH_loc + 1);
			if(nextTabLoc != string::npos)
				tmpCandiAlignmentNumStr = tmpStr.substr(IH_loc + 5, nextTabLoc - 1 - IH_loc - 5 + 1);
			else
				tmpCandiAlignmentNumStr = tmpStr.substr(IH_loc + 5);
		}
		else
			continue;
		tmpCandiAlignmentNum = atoi(tmpCandiAlignmentNumStr.c_str());
		if(tmpCandiAlignmentNum == 0)
			outputSamFile_unmap_ofs << tmpStr << endl;
		else if(tmpCandiAlignmentNum == 1)
			outputSamFile_unique_ofs << tmpStr << endl;
		else
			outputSamFile_multi_ofs << tmpStr << endl;
	}	
	sam_ifs.close();
	outputSamFile_unmap_ofs.close();
	outputSamFile_multi_ofs.close();
	outputSamFile_unique_ofs.close();
	return 0;	
}