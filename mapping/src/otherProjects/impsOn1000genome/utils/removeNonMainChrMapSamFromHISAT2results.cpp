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
		cout << "#0 Executable" << endl;
		cout << "#1 inputSam_raw" << endl;
		cout << "#2 outputSam" << endl;
		exit(1);
	}
	string inputSam_raw = argv[1];
	string outputSam = argv[2];
	ifstream raw_ifs(inputSam_raw.c_str());
	ofstream sam_ofs(outputSam.c_str());
	while(!raw_ifs.eof())
	{
		string tmpStr_1;
		getline(raw_ifs, tmpStr_1);
		if(tmpStr_1 == "")
			break;
		if(tmpStr_1.at(0) == '@') // header
		{
			if((tmpStr_1.substr(0,3) == "@HD")||(tmpStr_1.substr(0,3) == "@PG"))
				sam_ofs << tmpStr_1 << endl;
			else if(tmpStr_1.substr(0,9) == "@SQ\tSN:GL")
			{}
			else if(tmpStr_1.substr(0,9) == "@SQ\tSN:MT")
				sam_ofs << "@SQ\tSN:chrM" << tmpStr_1.substr(9) << endl;
			else
				sam_ofs << tmpStr_1.substr(0,7) << "chr" << tmpStr_1.substr(7) << endl;
		}
		else // alignment
		{
			string tmpStr_2;
			getline(raw_ifs, tmpStr_2);
			if(tmpStr_2 == "")
				break;
			int tabLoc_1_1 = tmpStr_1.find("\t");
			int tabLoc_1_2 = tmpStr_1.find("\t", tabLoc_1_1 + 1);
			int tabLoc_1_3 = tmpStr_1.find("\t", tabLoc_1_2 + 1);
			int tabLoc_2_1 = tmpStr_2.find("\t");
			int tabLoc_2_2 = tmpStr_2.find("\t", tabLoc_2_1 + 1);
			int tabLoc_2_3 = tmpStr_2.find("\t", tabLoc_2_2 + 1);			
			string tmpChrName_1 = tmpStr_1.substr(tabLoc_1_2 + 1, tabLoc_1_3 - tabLoc_1_2 - 1);
			string tmpChrName_2 = tmpStr_2.substr(tabLoc_2_2 + 1, tabLoc_2_3 - tabLoc_2_2 - 1);
			if(tmpChrName_1 != tmpChrName_2)
			{
				cout << "(tmpChrName_1 != tmpChrName_2)" << endl;
				cout << "tmpStr_1: " << endl << tmpStr_1 << endl << "tmpStr_2: " << endl << tmpStr_2 << endl;
				exit(1);
			}
			if(tmpChrName_1.length() >= 2)
			{
				if(tmpChrName_1.substr(0,2) == "GL")
				{}
				else if(tmpChrName_1.substr(0,2) == "MT")
				{
					sam_ofs << tmpStr_1.substr(0, tabLoc_1_2 + 1) << "chrM" << tmpStr_1.substr(tabLoc_1_3) << endl;
					sam_ofs << tmpStr_2.substr(0, tabLoc_2_2 + 1) << "chrM" << tmpStr_2.substr(tabLoc_2_3) << endl;					
				}
				else
				{
					sam_ofs << tmpStr_1.substr(0, tabLoc_1_2 + 1) << "chr" << tmpStr_1.substr(tabLoc_1_2 + 1) << endl;
					sam_ofs << tmpStr_2.substr(0, tabLoc_2_2 + 1) << "chr" << tmpStr_2.substr(tabLoc_2_2 + 1) << endl; 					
				}
			}
			else
			{
				sam_ofs << tmpStr_1.substr(0, tabLoc_1_2 + 1) << "chr" << tmpStr_1.substr(tabLoc_1_2 + 1) << endl;
				sam_ofs << tmpStr_2.substr(0, tabLoc_2_2 + 1) << "chr" << tmpStr_2.substr(tabLoc_2_2 + 1) << endl; 
			}
		}
	}
	raw_ifs.close();
	sam_ofs.close();
	return 0;
}