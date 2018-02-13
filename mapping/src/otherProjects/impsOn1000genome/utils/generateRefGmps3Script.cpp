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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputDataIdListFile" << endl;
		cout << "#2 inputDataDir" << endl;
		cout << "#3 outputResultsDir" << endl;
		cout << "#4 outputScript" << endl;
		exit(1);
	}
	string inputDataIdListFile = argv[1];
	string inputDataDir = argv[2]; inputDataDir += "/";
	string outputResultsDir = argv[3]; outputResultsDir += "/";
	string outputScript = argv[4];
	ofstream script_ofs(outputScript.c_str());

	ifstream idList_ifs(inputDataIdListFile.c_str());
	int id_num = 0;
	while(!idList_ifs.eof())
	{
		string tmpIdStr;
		getline(idList_ifs, tmpIdStr);
		if(tmpIdStr == "")
			break;
		id_num ++;	

		string tmpInputReadFile_end1 = 
		string tmpInputReadFile_end2 = 
		string tmpOutputResultsDir = 

	}
	script_ofs.close();
	return 0;
}