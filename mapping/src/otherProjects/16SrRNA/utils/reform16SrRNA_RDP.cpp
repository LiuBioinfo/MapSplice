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

// string refromRDP16SrRNAname(string& tmpName)
// {
// 	int firstTabLoc = tmpName.find(" ");
// 	int familyLoc = tmpName.find(";family;");
// 	int genusLoc = tmpName.find(";genus");
// 	if(firstTabLoc == string::npos)
// 	{
// 		cout << "firstTabLoc: " << firstTabLoc << endl;
// 		cout << "exiting ......" << endl;
// 		cout << "tmpName: " << tmpName << endl;
// 		exit(1);
// 	}
// 	string tmpSeqID = tmpName.substr(0,firstTabLoc);
// 	string tmpGenusName;
// 	if(familyLoc == string::npos)
// 	{
// 		//cout << "familyLoc: " << familyLoc << endl;
// 		//cout << "exiting ......" << endl;
// 		//cout << "tmpName: " << tmpName << endl;
// 		//exit(1);
// 		tmpGenusName = "NULL";
// 	}
// 	else
// 	{	
// 		if(genusLoc != string::npos)
// 			tmpGenusName = tmpName.substr(familyLoc + 8, genusLoc - 1 - familyLoc - 8 + 1);
// 		else
// 		{
// 			int tmpNameLength = tmpName.length();
// 			tmpGenusName = tmpName.substr(familyLoc + 8, tmpNameLength - 1 - familyLoc - 8 + 1);
// 			if(tmpGenusName == "")
// 			{
// 				cout << "exiting ......" << endl;
// 				cout << "tmpName: " << tmpName << endl;
// 				exit(1);
// 			}
// 		}
// 	}
// 	string tmpReformedName = tmpSeqID + ":" + tmpGenusName;
// 	return tmpReformedName;
// }

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputOriRDP16SrRNAfasta outputFilePrefix" << endl;
		exit(1);
	}
	unsigned int totalLength = 0;

	string inputOriRDP16SrRNAfasta = argv[1];
	string outputFilePrefix = argv[2];
	ifstream oriRDP16SrRNA_ifs(inputOriRDP16SrRNAfasta.c_str());
	string outputFastaFile = outputFilePrefix + "_seq.fa";
	string outputInfoFile = outputFilePrefix + "_info.txt";
	ofstream reformed_ofs(outputFastaFile.c_str());
	ofstream info_ofs(outputInfoFile.c_str());
	string name;
	string seq;
	getline(oriRDP16SrRNA_ifs, name);
	getline(oriRDP16SrRNA_ifs, seq);
	int tmpSeqNum = 0;
	while(!oriRDP16SrRNA_ifs.eof())
	{
		string tmpStr;
		getline(oriRDP16SrRNA_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == ">")
		{
			tmpSeqNum ++;
			//string reformedNameStr = refromRDP16SrRNAname(name);
			//reformed_ofs << reformedNameStr << endl;
			reformed_ofs << ">S_" << tmpSeqNum << endl;
			//name = tmpStr;
			reformed_ofs << seq << endl;
			info_ofs << "S_" << tmpSeqNum << endl;
			info_ofs << name << endl;
			unsigned int tmpLength = seq.length();
			totalLength += tmpLength;
			seq = "";
		}
		else
			seq += tmpStr;
	}
	cout << "totalLength: " << totalLength << endl;
	reformed_ofs.close();
	info_ofs.close();
	oriRDP16SrRNA_ifs.close();
	return 0;
}