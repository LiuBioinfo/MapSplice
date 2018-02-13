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
	if(argc != 4)
	{
		cout << "Executable inputChrSeqFilePath outputChrSeqFilePath renamedChrName" << endl;
		exit(1);
	}
	string inputChrSeqFilePath = argv[1];
	string outputChrSeqFilePath = argv[2];
	string renamedChrName = argv[3];
	ofstream chrSeq_ofs(outputChrSeqFilePath.c_str());
	ifstream chrSeq_ifs(inputChrSeqFilePath.c_str());
	string tmpChrNameStr;
	getline(chrSeq_ifs, tmpChrNameStr);
	chrSeq_ofs << ">" << renamedChrName << endl;
	while(!chrSeq_ifs.eof())
	{
		string tmpStr;
		getline(chrSeq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		chrSeq_ofs << tmpStr << endl;
	}
	chrSeq_ifs.close();
	chrSeq_ofs.close();
	return 0;
}