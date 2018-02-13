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
		cout << "Executable inputRawSortedKmerFile outputNewSortedKmerFile" << endl;
		exit(1);
	}
	string inputRawSortedKmerFile = argv[1];
	string outputNewSortedKmerFile = argv[2];
	ifstream raw_ifs(inputRawSortedKmerFile.c_str());
	ofstream new_ofs(outputNewSortedKmerFile.c_str());
	while(!raw_ifs.eof())
	{
		string tmpStr;
		getline(raw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpKmerStr = tmpStr.substr(0, tabLoc);
		new_ofs << tmpKmerStr << "\t35" << endl; 
	}
	raw_ifs.close();
	new_ofs.close();
	return 0;
}