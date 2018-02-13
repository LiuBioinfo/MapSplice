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
#include "../../../general/read_block_test.h"
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable outputDir KmerSetNum inputSortedKmer_union inputSortedKmer_set1 (inputSortedKmer_set2 ...)" << endl;
		exit(1);
	}
	string outputDir = argv[1];
	string KmerSetNumStr = argv[2];
	int KmerSetNum = atoi(KmerSetNumStr.c_str());
	string inputSortedKmer_union = argv[3];
	vector<string> inputSingleSetSortedKmerFileVec;
	for(int tmp = 4; tmp <= argc - 1; tmp++)
	{	
		string tmpInputSingleSetSortedKmerFile = argv[tmp];
		inputSingleSetSortedKmerFileVec.push_back(tmpInputSingleSetSortedKmerFile);
	}
	if(KmerSetNum != inputSingleSetSortedKmerFileVec.size())
	{
		cout << "error! (KmerSetNum != inputSingleSetSortedKmerFileVec.size())" << endl;
		exit(1);
	}

	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());

	string output_KmerClassProfile = outputDir + "/Kmer.class";
	string output_KmerAlienProfile = outputDir + "/Kmer.alien";
	string output_KmerRepetitiveProfile = outputDir + "/Kmer.repetitive";
	string output_KmerSetSpecificProfile = outputDir + "/Kmer.setSpecific";
	ifstream mergedKmer_ifs(inputSortedKmer_union.c_str());
	ofstream KmerClass_ofs(output_KmerClassProfile.c_str());
	ofstream KmerAlien_ofs(output_KmerAlienProfile.c_str());
	ofstream KmerRepetitive_ofs(output_KmerRepetitiveProfile.c_str());
	ofstream KmerSetSpecific_ofs(output_KmerSetSpecificProfile.c_str());

	cout << "start to generate subRegion specific kmers" << endl;
	cout << "start to initiate localRegionKmerIfsVec" << endl;
	vector<ifstream*> KmerReadableIfsVec;
	for(int tmp = 0; tmp < KmerSetNum; tmp++)
	{
		string tmp_kmer_sorted_file = inputSingleSetSortedKmerFileVec[tmp];
		ifstream *tmpSortedKmer_ifs = new ifstream(tmp_kmer_sorted_file.c_str());
		KmerReadableIfsVec.push_back(tmpSortedKmer_ifs);
	}
	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < KmerSetNum; tmp ++)
		kmer_file_end_bool_vec.push_back(false);
	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < KmerSetNum; tmp ++)
		currentKmerStrVec.push_back("");
	
	//cout << "head lines" << endl;
	for(int tmp = 0; tmp < KmerSetNum; tmp ++)
	{
		string tmpKmerFileStr;
		getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
		currentKmerStrVec[tmp] = tmpKmerFileStr;
	}

	while(!mergedKmer_ifs.eof())
	{
		string tmpMergedStr;
		getline(mergedKmer_ifs, tmpMergedStr);
		if(tmpMergedStr == "")
			break;
		int tmpKmerExistence_fileNum = 0;
		int tmpKmerExistence_lastFileId = -1;
		vector<bool> tmp_kmer_exist_in_sample_bool_vec;
		for(int tmp = 0; tmp < KmerSetNum; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		//int tmpKmerExistingClassFlag = 0;
		for(int tmp = 0; tmp < KmerSetNum; tmp ++)
		{
			if(tmpMergedStr == currentKmerStrVec[tmp])
			{
				tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
				tmpKmerExistence_fileNum ++;
				tmpKmerExistence_lastFileId = tmp + 1;
				//tmpKmerExistingClassFlag += pow(2, tmp);
				if(!kmer_file_end_bool_vec[tmp])
				{
					if((*KmerReadableIfsVec[tmp]).eof())
						kmer_file_end_bool_vec[tmp] = true;
					else
					{	
						string tmpKmerFileStr;
						getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
						if(tmpKmerFileStr == "")
							kmer_file_end_bool_vec[tmp] = true;
						else
							currentKmerStrVec[tmp] = tmpKmerFileStr;
					}
				}
			}
			else
				tmp_kmer_exist_in_sample_bool_vec[tmp] = false;
		}
		if(tmpKmerExistence_fileNum == 0) // alien Kmer
		{
			KmerClass_ofs << tmpMergedStr << "\t0" << endl;
			KmerAlien_ofs << tmpMergedStr << "\t0" << endl;
		}
		else if(tmpKmerExistence_fileNum == 1) // local region specific kmers
		{
			KmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistence_lastFileId << endl;
			KmerSetSpecific_ofs << tmpMergedStr << "\t" << tmpKmerExistence_lastFileId << endl;
		}
		else if(tmpKmerExistence_fileNum > 1) // multi local region kmers
		{
			KmerClass_ofs << tmpMergedStr << "\t" << KmerSetNum + 1 << endl;
			KmerRepetitive_ofs << tmpMergedStr << "\t" << KmerSetNum + 1 << endl;
		}
		else
		{
			cout << "error ! tmpKmerExistence_fileNum < 0 !" << endl;
			exit(1); 
		}
	}

	KmerAlien_ofs.close();
	KmerRepetitive_ofs.close();
	KmerSetSpecific_ofs.close();
	KmerClass_ofs.close();
	mergedKmer_ifs.close();
	for(int tmp = 0; tmp < KmerSetNum; tmp++)
	{
		(*KmerReadableIfsVec[tmp]).close();
		delete KmerReadableIfsVec[tmp];
	}	
	return 0;
}