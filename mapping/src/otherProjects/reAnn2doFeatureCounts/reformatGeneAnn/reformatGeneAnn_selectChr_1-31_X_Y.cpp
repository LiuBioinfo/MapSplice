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

bool checkChrName(string& raw_chrName)
{
	if((raw_chrName == "chr1")||(raw_chrName == "chr2")||(raw_chrName == "chr3")||(raw_chrName == "chr4")||(raw_chrName == "chr5")
		||(raw_chrName == "chr6")||(raw_chrName == "chr7")||(raw_chrName == "chr8")||(raw_chrName == "chr9")||(raw_chrName == "chr10")
		||(raw_chrName == "chr11")||(raw_chrName == "chr12")||(raw_chrName == "chr13")||(raw_chrName == "chr14")||(raw_chrName == "chr15")
		||(raw_chrName == "chr16")||(raw_chrName == "chr17")||(raw_chrName == "chr18")||(raw_chrName == "chr19")||(raw_chrName == "chr20")
		||(raw_chrName == "chr21")||(raw_chrName == "chr22")||(raw_chrName == "chr23")||(raw_chrName == "chr24")||(raw_chrName == "chr25")
		||(raw_chrName == "chr26")||(raw_chrName == "chr27")||(raw_chrName == "chr28")||(raw_chrName == "chr29")||(raw_chrName == "chr30")
		||(raw_chrName == "chr31")||(raw_chrName == "chrX")||(raw_chrName == "chrY"))
		return true;
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
		bool checkChrName_success_bool = checkChrName(tmpChrStr_raw);
		if(checkChrName_success_bool)
			reformattedGtf_ofs << tmpStr << endl;
	}
	gtf_ifs.close();
	reformattedGtf_ofs.close();
	return 0;
}