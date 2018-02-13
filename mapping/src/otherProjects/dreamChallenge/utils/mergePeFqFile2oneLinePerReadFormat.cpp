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
		cout << "#1 fq_1_file" << endl;
		cout << "#2 fq_2_file" << endl;
		cout << "#3 oneLinePerRead_file" << endl;
		exit(1);
	}
	string input_fq_1_file = argv[1];
	string input_fq_2_file = argv[2];
	string output_oneLinePerRead_file = argv[3];
	ifstream fq_1_ifs(input_fq_1_file.c_str());
	ifstream fq_2_ifs(input_fq_2_file.c_str());
	ofstream oneLinePerRead_ofs(output_oneLinePerRead_file.c_str());
	while((!fq_1_ifs.eof())&&(!fq_2_ifs.eof()))
	{
		string tmpId_1, tmpId_2;
		getline(fq_1_ifs, tmpId_1); 
		getline(fq_2_ifs, tmpId_2);
		if((tmpId_1 == "")||(tmpId_2 == ""))
			break;
		string tmpSeq_1, tmpSeq_2, tmpCom_1, tmpCom_2, tmpQua_1, tmpQua_2;
		getline(fq_1_ifs, tmpSeq_1); 
		getline(fq_2_ifs, tmpSeq_2);
		getline(fq_1_ifs, tmpCom_1); 
		getline(fq_2_ifs, tmpCom_2);
		getline(fq_1_ifs, tmpQua_1); 
		getline(fq_2_ifs, tmpQua_2);						
		oneLinePerRead_ofs << tmpId_1 << "\t" << tmpSeq_1 << "\t" << tmpId_2 << "\t" << tmpSeq_2 << endl;
	}
	oneLinePerRead_ofs.close();
	fq_1_ifs.close();
	fq_2_ifs.close();
	return 0;
}