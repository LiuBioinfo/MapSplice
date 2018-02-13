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
#include "../../../../general/read_block_test.h"
#include "../../../../general/bwtmap_info.h"
#include "../../../../general/DoubleAnchorScore.h"
#include "../../../../general/sbndm.h"
#include "../../../../general/splice_info.h"
#include "../../../../general/index_info.h"
time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable outputFile inputSortedKmerFile_union inputSortedKmerFile_1 inputSortedKmerFile_2 (inputSortedKmerFile_3 ...) " << endl;
		exit(1);
	}

	int fileNum = argc - 3;
	vector<string> fileVec;
	for(int tmp = 0; tmp < fileNum; tmp++)
		fileVec.push_back(argv[3 + tmp]);
	vector<ifstream*> KmerReadableIfsVec;
	for(int tmp = 0; tmp < fileNum; tmp++)
	{
		ifstream *tmpKmer_ifs = new ifstream(fileVec[tmp].c_str());
		KmerReadableIfsVec.push_back(tmpKmer_ifs);
	}

	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < fileNum; tmp ++)
		kmer_file_end_bool_vec.push_back(false);
	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < fileNum; tmp ++)
		currentKmerStrVec.push_back("");
	for(int tmp = 0; tmp < fileNum; tmp ++)
	{
		string tmpKmerFileStr;
		getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
		int tmpTabLoc = tmpKmerFileStr.find("\t");
		currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
	}

	string outputFile = argv[1];
	string inputSortedKmerFile_union = argv[2];	
	ifstream mergedKmer_ifs(inputSortedKmerFile_union.c_str());
	ofstream kmerClass_ofs(outputFile.c_str());
	unsigned int alienKmerNum = 0;
	unsigned int uniqueKmerNum = 0;
	unsigned int repetitiveKmerNum = 0;
	while(!mergedKmer_ifs.eof())
	{
		string tmpStr;
		getline(mergedKmer_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpMergedStr = tmpStr.substr(0, tabLoc);
		int tmpKmerExistence_fileNum = 0;
		int tmpKmerExistence_lastFileId = -1;
		vector<bool> tmp_kmer_exist_in_sample_bool_vec;
		for(int tmp = 0; tmp < fileNum; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		//int tmpKmerExistingClassFlag = 0;
		for(int tmp = 0; tmp < fileNum; tmp ++)
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
						{
							int tmpTabLoc = tmpKmerFileStr.find("\t");
							currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
						}
					}
				}
			}
			else
				tmp_kmer_exist_in_sample_bool_vec[tmp] = false;
		}
		if(tmpKmerExistence_fileNum == 0) // aliens
		{
			kmerClass_ofs << tmpMergedStr << "\t" << 0 << endl;
			alienKmerNum ++;
		}
		else if(tmpKmerExistence_fileNum == 1) // local region specific kmers
		{
			kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistence_lastFileId << endl;
			uniqueKmerNum ++;
		}
		else if(tmpKmerExistence_fileNum > 1) // multi local region kmers
		{
			kmerClass_ofs << tmpMergedStr << "\t" << fileNum + 1 << endl;
			repetitiveKmerNum ++;
		}
		else
		{
			cout << "error ! tmpKmerExistence_fileNum < 0 !" << endl;
			exit(1); 
		}
	}
	cout << "alienKmerNum: " << alienKmerNum << endl;
	cout << "uniqueKmerNum: " << uniqueKmerNum << endl;
	cout << "repetitiveKmerNum: " << repetitiveKmerNum << endl;
	cout << endl;
	cout << "solid kmerNum: " << uniqueKmerNum + repetitiveKmerNum << endl;
	cout << "total KmerNum: " << alienKmerNum + uniqueKmerNum + repetitiveKmerNum << endl;
	for(int tmp = 0; tmp < fileNum; tmp++)
	{
		(*KmerReadableIfsVec[tmp]).close();
		delete KmerReadableIfsVec[tmp];
	}
	mergedKmer_ifs.close();
	kmerClass_ofs.close();
	return 0;
}