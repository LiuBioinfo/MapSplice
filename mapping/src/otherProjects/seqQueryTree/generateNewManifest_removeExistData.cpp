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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

void get_strVec_bamFileNameVec_from_manifestFile(string& manifestFile, 
	vector<string>& strVec, vector<string>& bamFileNameVec)
{
	ifstream manifest_ifs(manifestFile.c_str());
	while(!manifest_ifs.eof())
	{
		string tmpStr;
		getline(manifest_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		string tmpFileName = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		strVec.push_back(tmpStr);
		bamFileNameVec.push_back(tmpFileName);
	}
	manifest_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_ori_manifest_file" << endl;
		cout << "#2 input_exist_bam_filename_list (filename! not Path!)" << endl;
		cout << "#3 output_new_manifest_file" << endl;
		exit(1);
	}

	string input_ori_manifest_file = argv[1];
	string input_exist_bam_filename_list = argv[2];
	string output_new_manifest_file = argv[3];

	vector<string> oriManifestStrVec;
	vector<string> oriManifestStrVec_bamFileNameVec;
	get_strVec_bamFileNameVec_from_manifestFile(input_ori_manifest_file, 
		oriManifestStrVec, oriManifestStrVec_bamFileNameVec);

	vector<string> existBamFileVec;
	ifstream exist_ifs(input_exist_bam_filename_list.c_str());
	while(!exist_ifs.eof())
	{
		string tmpStr;
		getline(exist_ifs, tmpStr);
		if(tmpStr == "")
			break;
		existBamFileVec.push_back(tmpStr);
	}
	exist_ifs.close();

	ofstream new_ofs(output_new_manifest_file.c_str());
	for(int tmp = 0; tmp < oriManifestStrVec.size(); tmp++)
	{
		string tmp_oriBamFile = oriManifestStrVec_bamFileNameVec[tmp];
		bool tmp_exist_bool = false;
		for(int tmp2 = 0; tmp2 < existBamFileVec.size(); tmp2++)
		{
			if(tmp_oriBamFile == existBamFileVec[tmp2])
			{
				tmp_exist_bool = true;
				break;
			}
		}
		if(!tmp_exist_bool)
			new_ofs << oriManifestStrVec[tmp] << endl;
	}
	new_ofs.close();
	return 0;
}