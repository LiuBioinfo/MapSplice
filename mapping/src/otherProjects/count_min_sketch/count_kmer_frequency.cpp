// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
# include <iostream>
# include <cstdlib>
# include <cmath>
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
# include "count_min_sketch/count_min_sketch.hpp"
using namespace std;

time_t nowtime;
struct tm *local;

int getKmerCount(string& tmpStr)
{
	string countStr = tmpStr.substr(1);
	return atoi(countStr.c_str());
}
// run basic tests on CountMinSketch
int main(int argc, char **argv) 
{
	if(argc != 3)
	{
		cout << "Executable inputJellyFishKmerCountFile querySequence" << endl;
		exit(1);
	}
  
	CountMinSketch c(0.01, 0.01);
	string inputJellyFishFile = argv[1];
	ifstream jellyfish_ifs(inputJellyFishFile.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to build Kmer sketch......" << endl << endl;
	int tmpLineNO = 0;
	while(!jellyfish_ifs.eof())
	{
		string tmpStr_id;
		getline(jellyfish_ifs, tmpStr_id);
		if(tmpStr_id == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 100000;
		if(tmpLineNO == tmpThousandIndex * 100000)
			cout << "Processed Line #: " << tmpLineNO << endl;	
		string tmpStr_seq;
		getline(jellyfish_ifs, tmpStr_seq);
		int tmpKmerCount = getKmerCount(tmpStr_id);
		cout << "tmpKmerCount: " << tmpKmerCount << endl;
		c.update(tmpStr_seq.c_str(), tmpKmerCount);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of building Kmer sketch......" << endl << endl; 	
	string query = argv[2];
	cout << "query: " << query << endl;
	cout << "esitimated frequency: " << c.estimate(query.c_str()) << endl;

	// c.update("AAAAA", 3);
	// c.update("CCCCC", 2);
	// c.update("GGGGG", 1);
	// cout << "c.estimate('AAAAAA'): " << c.estimate("AAAAAA") << endl;
	// cout << "c.estimate('AAAAA'): " << c.estimate("AAAAA") << endl;
	// cout << "c.estimate('AAAVTC'): " << c.estimate("AAGTV") << endl;

	jellyfish_ifs.close();
    return 0;
}