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

void parseDreamChallengeBedPeStr(string& tmpDreamChallengeBedPeStr, 
	string& chrName_A, int& pos_A, string& chrName_B, int& pos_B,
	string& geneId_A, string& geneId_B, string& strand_A, string& strand_B)
{
	vector<string> tmpFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 9; tmp++)
	{
		int tabLoc = tmpDreamChallengeBedPeStr.find("\t", startLoc);
		string tmpField = tmpDreamChallengeBedPeStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpDreamChallengeBedPeStr.substr(startLoc));
	chrName_A = "chr" + tmpFieldVec[0];
	chrName_B = "chr" + tmpFieldVec[3];
	string genePairStr = tmpFieldVec[6];
	int lineLoc = genePairStr.find("-");
	geneId_A = genePairStr.substr(0, lineLoc);
	geneId_B = genePairStr.substr(lineLoc + 1);
	strand_A = tmpFieldVec[8];
	strand_B = tmpFieldVec[9];
	
	string posStr_A_1 = tmpFieldVec[1];
	string posStr_A_2 = tmpFieldVec[2];
	string posStr_B_1 = tmpFieldVec[4];
	string posStr_B_2 = tmpFieldVec[5];
	if(strand_A == "+")
		pos_A = atoi(posStr_A_2.c_str());
	else if(strand_A == "-")
		pos_A = atoi(posStr_A_1.c_str());
	else
	{
		cout << "invalid strand_A: " << strand_A << endl;
		exit(1);
	}
	if(strand_B == "+")
		pos_B = atoi(posStr_B_1.c_str());
	else if(strand_B == "-")
		pos_B = atoi(posStr_B_2.c_str());
	else
	{
		cout << "invalid strand_B: " << strand_B << endl;
		exit(1);
	}
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputDreamChallengeBedPeFile outputReformattedFile" << endl;
		exit(1);
	}
	string inputDreamChallengeBedPeFile = argv[1];
	string outputReformattedFile = argv[2];

	ifstream dreamChallengeBedPe_ifs(inputDreamChallengeBedPeFile.c_str());
	ofstream reformatted_ofs(outputReformattedFile.c_str());
	while(!dreamChallengeBedPe_ifs.eof())
	{
		string tmpStr;
		getline(dreamChallengeBedPe_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string chrName_A, chrName_B, geneId_A, geneId_B, strand_A, strand_B;
		int pos_A, pos_B;
		parseDreamChallengeBedPeStr(tmpStr, chrName_A, pos_A, chrName_B, pos_B,
			geneId_A, geneId_B, strand_A, strand_B);
		reformatted_ofs << chrName_A << "\t" << pos_A << "\t" << chrName_B 
			<< "\t" << pos_B << "\t" << geneId_A << "\t" << geneId_B << "\t" 
			<< strand_A << "\t" << strand_B << endl;
	}
	reformatted_ofs.close();
	dreamChallengeBedPe_ifs.close();
	return 0;
}