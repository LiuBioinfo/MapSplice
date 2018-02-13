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
#include "../../../general/index_info.h"
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 oneLinePerRead_file" << endl;
		cout << "#2 fa_1_file" << endl;
		cout << "#3 fa_2_file" << endl;
		exit(1);
	}
	string input_oneLinePerRead_file = argv[1];
	string output_fa_1_file = argv[2];
	string output_fa_2_file = argv[3];
	
	ifstream oneLinePerRead_ifs(input_oneLinePerRead_file.c_str());
	ofstream fa_1_ofs(output_fa_1_file.c_str());
	ofstream fa_2_ofs(output_fa_2_file.c_str());
	
	while(!oneLinePerRead_ifs.eof())
	{
		string tmpStr;
		getline(oneLinePerRead_ifs, tmpStr); 
		if(tmpStr == "")
			break;
		string tmpId_1, tmpId_2, tmpSeq_1, tmpSeq_2;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		tmpId_1 = tmpStr.substr(0, tabLoc_1);
		tmpSeq_1 = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		tmpId_2 = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		tmpSeq_2 = tmpStr.substr(tabLoc_3 + 1);
		fa_1_ofs << tmpId_1 << endl << tmpSeq_1 << endl;
		fa_2_ofs << tmpId_2 << endl << tmpSeq_2 << endl;
	}

	fa_1_ofs.close();
	fa_2_ofs.close();
	oneLinePerRead_ifs.close();
	return 0;
}