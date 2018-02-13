#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFastqFile outputReadInfo" << endl;
		exit(1);
	}
	int readLengthMax = 100;
	vector<int> readNumVecWithDiffLength;
	for(int tmp = 0; tmp < readLengthMax; tmp++)
		readNumVecWithDiffLength.push_back(0);

	string inputFastqFile = argv[1];
	string outputReadInfoFile = argv[2];
	ifstream fq_ifs(inputFastqFile.c_str());
	ofstream readInfo_ofs(outputReadInfoFile.c_str());
	int readNum = 0;
	while(!fq_ifs.eof())
	{
		string str_1, str_2, str_3, str_4;
		getline(fq_ifs, str_1);
		if(str_1 == "")
			break;
		readNum ++;
		getline(fq_ifs, str_2);
		getline(fq_ifs, str_3);
		getline(fq_ifs, str_4);
		int tmpReadLength = str_2.size();
		if(tmpReadLength > readLengthMax)
			readNumVecWithDiffLength[readLengthMax - 1] ++;
		else
			readNumVecWithDiffLength[tmpReadLength - 1] ++;
	}
	readInfo_ofs << "TotalReadNum: " << readNum << endl;
	for(int tmp = 0; tmp < readLengthMax; tmp++)
		readInfo_ofs << "readLength[" << tmp + 1 << "]:\t" << readNumVecWithDiffLength[tmp] << endl;
	fq_ifs.close();
	readInfo_ofs.close();
	return 0;
}