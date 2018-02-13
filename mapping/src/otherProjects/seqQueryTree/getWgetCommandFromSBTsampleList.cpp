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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSBTsampleList" << endl;
		cout << "#2 outputWgetCommandFile" << endl;
		exit(1);
	}
	string inputSBTsampleList = argv[1];
	string outputWgetCommandFile = argv[2];
	vector<string> SRRidVec;
	ifstream SBTsampleList_ifs(inputSBTsampleList.c_str());
	while(!SBTsampleList_ifs.eof())
	{
		string tmpStr;
		getline(SBTsampleList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int space_loc_1 = tmpStr.find(" ");
		int space_loc_2 = tmpStr.find(" ", space_loc_1 + 1);
		string tmpSRRid = tmpStr.substr(space_loc_1 + 1, space_loc_2 - space_loc_1 - 1);
		SRRidVec.push_back(tmpSRRid);
	}
	SBTsampleList_ifs.close();

	ofstream wget_ofs(outputWgetCommandFile.c_str());
	for(int tmp = 0; tmp < SRRidVec.size(); tmp++)
	{
		string tmpSRRid = SRRidVec[tmp];
		wget_ofs << "https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/" 
			<< tmpSRRid.substr(0,6) << "/" << tmpSRRid << "/" << tmpSRRid << ".sra" << endl;
	}
	wget_ofs.close();
	return 0;
}