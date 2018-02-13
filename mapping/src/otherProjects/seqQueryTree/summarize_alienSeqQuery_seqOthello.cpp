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
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSeqOthelloResults" << endl;
		cout << "#2 outputSummary" << endl;
		cout << "#3 hitNumPercMin" << endl;
		exit(1);
	}
	string inputSeqOthelloResults = argv[1];
	string outputSummary = argv[2];
	string hitNumPercMinStr = argv[3];
	double hitNumPercMin = atof(hitNumPercMinStr.c_str());
	ifstream seqOthello_ifs(inputSeqOthelloResults.c_str());
	ofstream summary_ofs(outputSummary.c_str());
	int hit_query = 0;
	int hit_number = 0;
	while(!seqOthello_ifs.eof())
	{
		string tmpSeq;
		getline(seqOthello_ifs, tmpSeq);
		if(tmpSeq == "")
			break;
		int tmpSeqLength = tmpSeq.length();
		vector<double> tmpHitPercVec;
		for(int tmp = 0; tmp < 152; tmp++)
		{
			string tmpKmerHitNumStr;
			getline(seqOthello_ifs, tmpKmerHitNumStr);
			int tmpKmerHitNum = atoi(tmpKmerHitNumStr.c_str());
			double tmpKmerHitPerc = (double)tmpKmerHitNum/(double)tmpSeqLength;
			tmpHitPercVec.push_back(tmpKmerHitPerc);
		}
		bool tmp_query_hit_or_not_bool = false;
		for(int tmp = 0; tmp < 148; tmp++)
		{
			if(tmpHitPercVec[tmp] >= hitNumPercMin)
			{
				tmp_query_hit_or_not_bool = true;
				hit_number ++;
			}
		}
		if(tmp_query_hit_or_not_bool)
			hit_query ++;
	}

	summary_ofs << "hit_query: " << hit_query << endl;
	summary_ofs << "hit_number: " << hit_number << endl;
	summary_ofs.close();
	seqOthello_ifs.close();
	return 0;
}