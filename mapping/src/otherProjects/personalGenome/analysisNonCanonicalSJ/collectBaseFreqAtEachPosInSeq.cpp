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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputSeqFile outputBaseFreqFile seqLength" << endl;
		exit(1);
	}
	string seqLengthStr = argv[3];
	int seqLength = atoi(seqLengthStr.c_str());
	vector<int> baseFreqVec_A;
	vector<int> baseFreqVec_C;
	vector<int> baseFreqVec_G;
	vector<int> baseFreqVec_T;
	vector<int> baseFreqVec_others;
	for(int tmp = 0; tmp < seqLength; tmp++)
	{
		baseFreqVec_A.push_back(0);
		baseFreqVec_C.push_back(0);
		baseFreqVec_G.push_back(0);
		baseFreqVec_T.push_back(0);
		baseFreqVec_others.push_back(0);		
	}
	string inputSeqFile = argv[1];
	string outputBaseFreqFile = argv[2];
	ifstream seq_ifs(inputSeqFile.c_str());
	ofstream baseFreq_ofs(outputBaseFreqFile.c_str());
	int seqNum = 0;
	while(!seq_ifs.eof())
	{
		string tmpStr;
		getline(seq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		seqNum ++;
		for(int tmp = 0; tmp < seqLength; tmp++)
		{
			char tmpBaseChar = tmpStr.at(tmp);
			if(tmpBaseChar == 'A')
				baseFreqVec_A[tmp] ++;
			else if(tmpBaseChar == 'C')
				baseFreqVec_C[tmp] ++;
			else if(tmpBaseChar == 'G')
				baseFreqVec_G[tmp] ++;
			else if(tmpBaseChar == 'T')
				baseFreqVec_T[tmp] ++;
			else
				baseFreqVec_others[tmp] ++;
		}
	}
	cout << "seqNum: " << seqNum << endl;
	for(int tmp = 0; tmp < seqLength; tmp++)
	{
		double tmpRelativeFreq_A_double = (double)baseFreqVec_A[tmp]*100/(double)seqNum;
		int tmpRelativeFreq_A_int = tmpRelativeFreq_A_double;
		double tmpRelativeFreq_C_double = (double)baseFreqVec_C[tmp]*100/(double)seqNum;
		int tmpRelativeFreq_C_int = tmpRelativeFreq_C_double;		
		double tmpRelativeFreq_G_double = (double)baseFreqVec_G[tmp]*100/(double)seqNum;
		int tmpRelativeFreq_G_int = tmpRelativeFreq_G_double;
		double tmpRelativeFreq_T_double = (double)baseFreqVec_T[tmp]*100/(double)seqNum;
		int tmpRelativeFreq_T_int = tmpRelativeFreq_T_double;
		double tmpRelativeFreq_others_double = (double)baseFreqVec_others[tmp]*100/(double)seqNum;
		int tmpRelativeFreq_others_int = tmpRelativeFreq_others_double;
		baseFreq_ofs << tmp + 1 << "\t" << tmpRelativeFreq_A_int << "\t" << tmpRelativeFreq_C_int << "\t" 
			<< tmpRelativeFreq_G_int << "\t" << tmpRelativeFreq_T_int << "\t" << tmpRelativeFreq_others_int << endl;
	}
	seq_ifs.close();
	baseFreq_ofs.close();
	return 0;
}