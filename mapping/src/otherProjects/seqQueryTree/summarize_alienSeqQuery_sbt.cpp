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
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSbtResults" << endl;
		cout << "#2 outputSummary" << endl;
		exit(1);
	}
	string inputSbtResults = argv[1];
	string outputSummary = argv[2];
	ifstream sbtResults_ifs(inputSbtResults.c_str());
	ofstream summary_ofs(outputSummary.c_str());
	int hit_query = 0;
	int hit_num = 0;
	while(!sbtResults_ifs.eof())
	{
		string tmpId;
		getline(sbtResults_ifs, tmpId);
		if(tmpId == "")
			break;
		int tabLoc = tmpId.find(" ");
		int tmpHitNum = atoi((tmpId.substr(tabLoc + 1)).c_str());
		if(tmpHitNum > 0)
		{
			hit_query ++;
			hit_num += tmpHitNum;
			for(int tmp = 0; tmp < tmpHitNum; tmp++)
			{
				string tmpHitId;
				getline(sbtResults_ifs, tmpHitId);
			}
		}
	}
	summary_ofs << "hit_query: " << hit_query << endl;
	summary_ofs << "hit_num: " << hit_num;
	sbtResults_ifs.close();
	summary_ofs.close();
	return 0;
}