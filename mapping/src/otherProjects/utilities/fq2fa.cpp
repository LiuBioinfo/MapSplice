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

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFqFile outputFaFile" << endl;
		exit(1);
	}
	string inputFqFile = argv[1];
	string outputFaFile = argv[2];
	ifstream fq_ifs(inputFqFile.c_str());	
	ofstream fa_ofs(outputFaFile.c_str());
	while(!fq_ifs.eof())
	{	
		string tmpName;
		getline(fq_ifs, tmpName);
		if(tmpName == "")
			break;
		string tmpSeq;
		getline(fq_ifs, tmpSeq);
		fa_ofs << ">" << tmpName.substr(1) << endl << tmpSeq << endl;
		string tmpComment, tmpQual;
		getline(fq_ifs, tmpComment);
		getline(fq_ifs, tmpQual);
	}
	fa_ofs.close();
	fq_ifs.close();
	return 0;
}