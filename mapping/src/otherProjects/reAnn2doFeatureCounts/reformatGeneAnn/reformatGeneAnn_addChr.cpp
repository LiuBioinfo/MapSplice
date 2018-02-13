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

time_t nowtime;
struct tm *local;

bool reformatChrName(string& raw_chrName, string& new_chrName)
{
	if((raw_chrName == "1")||(raw_chrName == "2")||(raw_chrName == "3")||(raw_chrName == "4")||(raw_chrName == "5")
		||(raw_chrName == "6")||(raw_chrName == "7")||(raw_chrName == "8")||(raw_chrName == "9")||(raw_chrName == "10")										
		||(raw_chrName == "11")||(raw_chrName == "12")||(raw_chrName == "13")||(raw_chrName == "14")||(raw_chrName == "15")
		||(raw_chrName == "16")||(raw_chrName == "17")||(raw_chrName == "18")||(raw_chrName == "19")||(raw_chrName == "20")
		||(raw_chrName == "21")||(raw_chrName == "22")||(raw_chrName == "23")||(raw_chrName == "24")||(raw_chrName == "25")
		||(raw_chrName == "26")||(raw_chrName == "27")||(raw_chrName == "28")||(raw_chrName == "29")||(raw_chrName == "30")
		||(raw_chrName == "31")||(raw_chrName == "X")||(raw_chrName == "Y"))
	{
		new_chrName = "chr" + raw_chrName;
		return true;
	}
	else
		return false;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputGTF outputReformattedGTF" << endl;
		exit(1);
	}
	string inputGTF = argv[1];
	string outputReformattedGTF = argv[2];
	ifstream gtf_ifs(inputGTF.c_str());
	ofstream reformattedGtf_ofs(outputReformattedGTF.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1st = tmpStr.find("\t");
		if(tabLoc_1st == string::npos)
			continue;
		string tmpChrStr_raw = tmpStr.substr(0, tabLoc_1st);
		string tmpGtfStr_others_raw = tmpStr.substr(tabLoc_1st + 1);
		string tmpChrStr_new;
		bool reformatChrName_success_bool = reformatChrName(tmpChrStr_raw, tmpChrStr_new);
		if(reformatChrName_success_bool)
			reformattedGtf_ofs << tmpChrStr_new << "\t" << tmpGtfStr_others_raw << endl;
	}
	gtf_ifs.close();
	reformattedGtf_ofs.close();
	return 0;
}