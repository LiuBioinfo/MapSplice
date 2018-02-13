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
		cout << "#1 BeersSNPfile" << endl;
		cout << "#2 hisat2formatSNPfile" << endl;
		exit(1);
	}
	string BeersSNPfile = argv[1];
	string hisat2formatSNPfile = argv[2];
	ifstream beersSNP_ifs(BeersSNPfile.c_str());
	ofstream hisat2SNP_ofs(hisat2formatSNPfile.c_str());
	int snp_num = 0;
	while(!beersSNP_ifs.eof())
	{
		string tmpStr;
		getline(beersSNP_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpChrName = tmpStr.substr(0, tabLoc_1);
		string tmpChrPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tmpChrPos = atoi(tmpChrPosStr.c_str());
		string tmpRefBase = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpAltBase = tmpStr.substr(tabLoc_3 + 1);
		snp_num ++;
		hisat2SNP_ofs << "SNP_" << snp_num << "\tsingle\t" << tmpChrName << "\t" << tmpChrPos - 1 << "\t" << tmpAltBase << endl;
	}
	beersSNP_ifs.close();
	hisat2SNP_ofs.close();
	return 0;
}