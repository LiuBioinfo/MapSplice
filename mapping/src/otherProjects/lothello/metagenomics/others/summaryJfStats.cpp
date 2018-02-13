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
		cout << "Executable inputStatsFileListFile outputFile" << endl;
		exit(1);
	}
	string inputStatsFileListFile = argv[1];
	string output_stats = argv[2];

	vector<string> statsFileVec;
	ifstream statsFile_ifs(inputStatsFileListFile.c_str());
	while(!statsFile_ifs.eof())
	{
		string tmpStr;
		getline(statsFile_ifs, tmpStr);
		if(tmpStr == "")
			break;
		statsFileVec.push_back(tmpStr);
	}
	statsFile_ifs.close();

	ofstream summary_ofs(output_stats.c_str());
	summary_ofs << "Stats_file\tUnique\tDistinct\tTotal\tMax_count" << endl;
	for(int tmp = 0; tmp < statsFileVec.size(); tmp++)
	{
		string tmpStatsFile = statsFileVec[tmp];
		ifstream tmpStats_ifs(tmpStatsFile.c_str());
		string unique_info_str, distinct_info_str, total_info_str, maxCount_info_str;
		getline(tmpStats_ifs, unique_info_str);
		string unique_num_str = unique_info_str.substr(11);
		getline(tmpStats_ifs, distinct_info_str);
		string distinct_num_str = distinct_info_str.substr(11);
		getline(tmpStats_ifs, total_info_str);
		string total_num_str = total_info_str.substr(11);
		getline(tmpStats_ifs, maxCount_info_str);
		string maxCount_num_str = maxCount_info_str.substr(11);
		summary_ofs << tmpStatsFile << "\t" << unique_num_str << "\t" << distinct_num_str 
			<< "\t" << total_num_str << "\t" << maxCount_num_str << endl;
		tmpStats_ifs.close();
	}
	summary_ofs.close();
	return 0;
}