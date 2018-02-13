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
		cout << "#1 inputReadIdList" << endl;
		cout << "#2 inputFq_1" << endl;
		cout << "#3 inputFq_2" << endl;
		cout << "#4 outputPrefix" << endl;
		exit(1);
	}
	string inputReadIdList = argv[1];
	string inputFq_1 = argv[2];
	string inputFq_2 = argv[3];
	string outputPrefix = argv[4];
	string outputFq_1 = outputPrefix + ".1.fq";
	string outputFq_2 = outputPrefix + ".2.fq";
	
	ifstream id_ifs(inputReadIdList.c_str());
	set<string> idSet;
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		if(tabLoc == string::npos)
			idSet.insert(tmpStr);
		else
			idSet.insert(tmpStr.substr(0, tabLoc));
	}
	id_ifs.close();

	ifstream fq_1_ifs(inputFq_1.c_str());
	ifstream fq_2_ifs(inputFq_2.c_str());
	ofstream fq_1_ofs(outputFq_1.c_str());
	ofstream fq_2_ofs(outputFq_2.c_str());
	while((!fq_1_ifs.eof())&&(!fq_2_ifs.eof()))
	{
		string tmpStr_1_1, tmpStr_2_1;
		getline(fq_1_ifs, tmpStr_1_1);
		getline(fq_2_ifs, tmpStr_2_1);
		if((tmpStr_1_1 == "")||(tmpStr_2_1 == ""))
			break;
		string tmpStr_1_2, tmpStr_2_2, tmpStr_1_3, tmpStr_2_3, tmpStr_1_4, tmpStr_2_4;
		getline(fq_1_ifs, tmpStr_1_2);
		getline(fq_1_ifs, tmpStr_1_3);
		getline(fq_1_ifs, tmpStr_1_4);
		getline(fq_2_ifs, tmpStr_2_2);
		getline(fq_2_ifs, tmpStr_2_3);
		getline(fq_2_ifs, tmpStr_2_4);
		if(idSet.find(tmpStr_1_1) != idSet.end()) // found
		{
			fq_1_ofs << tmpStr_1_1 << endl << tmpStr_1_2 << endl << tmpStr_1_3 << endl << tmpStr_1_4 << endl;
			fq_2_ofs << tmpStr_2_1 << endl << tmpStr_2_2 << endl << tmpStr_2_3 << endl << tmpStr_2_4 << endl;
		}
	}
	fq_1_ifs.close();
	fq_2_ifs.close();
	fq_1_ofs.close();
	fq_2_ofs.close();	
	return 0;
}