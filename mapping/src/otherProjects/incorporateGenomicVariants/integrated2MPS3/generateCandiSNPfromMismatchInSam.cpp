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

void parseMismatchInfo(string& tmpStr, string& tmpChrName, int& tmpChrPos, string& tmpRefBase, int& tmpTotalMismatchNum,
		int& tmpMismatchNum_A, int& tmpMismatchNum_C, int& tmpMismatchNum_G, int& tmpMismatchNum_T, int& tmpMismatchNum_N)
{
	int tabLoc_1 = tmpStr.find("\t");
	int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
	int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
	int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
	int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
	int tabLoc_8 = tmpStr.find("\t", tabLoc_7 + 1);
	tmpChrName = tmpStr.substr(0, tabLoc_1);
	string tmpChrPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	tmpChrPos = atoi(tmpChrPosStr.c_str());
	tmpRefBase = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	string tmpTotalMismatchNumStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	tmpTotalMismatchNum = atoi(tmpTotalMismatchNumStr.c_str());
	string tmpMismatchNum_A_str = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
	//cout << "tmpMismatchNum_A_str: " << tmpMismatchNum_A_str << endl;
	string tmpMismatchNum_C_str = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
	string tmpMismatchNum_G_str = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
	string tmpMismatchNum_T_str = tmpStr.substr(tabLoc_7 + 1, tabLoc_8 - tabLoc_7 - 1);
	string tmpMismatchNum_N_str = tmpStr.substr(tabLoc_8 + 1);
	tmpMismatchNum_A = atoi(tmpMismatchNum_A_str.c_str());
	//cout << "tmpMismatchNum_A: " << tmpMismatchNum_A << endl;
	tmpMismatchNum_C = atoi(tmpMismatchNum_C_str.c_str());
	tmpMismatchNum_G = atoi(tmpMismatchNum_G_str.c_str());
	tmpMismatchNum_T = atoi(tmpMismatchNum_T_str.c_str());
	tmpMismatchNum_N = atoi(tmpMismatchNum_N_str.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputMismatchFromSamFile outputCandiSNPfile supNumMin ratioMin" << endl;
		exit(1);
	}
	string supNumMinStr = argv[3];
	string ratioMinStr = argv[4];
	int supNumMin = atoi(supNumMinStr.c_str());
	double ratioMin = atof(ratioMinStr.c_str());

	string inputMismatchFromSamFile = argv[1];
	string outputCandiSNPfile = argv[2];
	ifstream mismatch_ifs(inputMismatchFromSamFile.c_str());
	ofstream candiSNP_ofs(outputCandiSNPfile.c_str());
	while(!mismatch_ifs.eof())
	{
		string tmpStr;
		getline(mismatch_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpChrName;
		int tmpChrPos;
		string tmpRefBase;
		int tmpTotalMismatchNum;
		int tmpMismatchNum_A;
		int tmpMismatchNum_C;
		int tmpMismatchNum_G;
		int tmpMismatchNum_T;		
		int tmpMismatchNum_N;
		parseMismatchInfo(tmpStr, tmpChrName, tmpChrPos, tmpRefBase, tmpTotalMismatchNum,
			tmpMismatchNum_A, tmpMismatchNum_C, tmpMismatchNum_G, tmpMismatchNum_T, tmpMismatchNum_N);
		//cout << "tmpMismatchNum_A: " << tmpMismatchNum_A << endl;
		double tmpMismatchRatio_A = (double)tmpMismatchNum_A / (double)tmpTotalMismatchNum;
		double tmpMismatchRatio_C = (double)tmpMismatchNum_C / (double)tmpTotalMismatchNum;
		double tmpMismatchRatio_G = (double)tmpMismatchNum_G / (double)tmpTotalMismatchNum;
		double tmpMismatchRatio_T = (double)tmpMismatchNum_T / (double)tmpTotalMismatchNum;
		//cout << "tmpMismatchNum_A: " << tmpMismatchNum_A << endl;
		//cout << "tmpMismatchRatio_A: " << tmpMismatchRatio_A << endl;
		if((tmpMismatchNum_A >= supNumMin)&&(tmpMismatchRatio_A >= ratioMin))
			candiSNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tA" << endl;
		else if((tmpMismatchNum_C >= supNumMin)&&(tmpMismatchRatio_C >= ratioMin))
			candiSNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tC" << endl;
		else if((tmpMismatchNum_G >= supNumMin)&&(tmpMismatchRatio_G >= ratioMin))
			candiSNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tG" << endl;
		else if((tmpMismatchNum_T >= supNumMin)&&(tmpMismatchRatio_T >= ratioMin))
			candiSNP_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\tT" << endl;				
		else
		{}
	}
	candiSNP_ofs.close();
	mismatch_ifs.close();
	return 0;
}