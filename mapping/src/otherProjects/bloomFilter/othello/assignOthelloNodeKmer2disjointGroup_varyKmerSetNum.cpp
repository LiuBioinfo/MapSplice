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
	if(argc < 6)
	{
		cout << "Executable inputMergedKmerReadableFile outputMergedKmerClassInfoFile BranchNum KmerReadableFileNum ";
		cout << "KmerReadableFile_1 KmerReadableFile_2 ..." << endl;
		exit(1);
	}
	string branchNumStr = argv[3];
	int branchNum = atoi(branchNumStr.c_str());
	string KmerReadableFileNumStr = argv[4];
	int KmerReadableFileNum = atoi(KmerReadableFileNumStr.c_str());
	int providedKmerReadableFileNum = argc - 5;
	if((branchNum > 64)||(branchNum < KmerReadableFileNum)||(KmerReadableFileNum != providedKmerReadableFileNum))
	{
		cout << "error ! ((branchNum > 64)||(branchNum < KmerReadableFileNum)||(KmerReadableFileNum != providedKmerReadableFileNum)) " << endl;
		exit(1);
	}

	vector<string> KmerReadableFilePathVec;
	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		KmerReadableFilePathVec.push_back(argv[5 + tmp]);
	vector<ifstream*> KmerReadableIfsVec;
	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
	{
		ifstream tmpKmerReadableIfs(KmerReadableFilePathVec[tmp].c_str());
		KmerReadableIfsVec.push_back(&tmpKmerReadableIfs);
	}

	string inputMergedKmerReadableFile = argv[1];
	string outputMergedKmerClassInfoFile = argv[2];

	ofstream kmerClass_ofs(outputMergedKmerClassInfoFile.c_str());
	ifstream mergedKmer_ifs(inputMergedKmerReadableFile.c_str());	

	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
		kmer_file_end_bool_vec.push_back(false);

	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
		currentKmerStrVec.push_back("");

	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
	{
		string tmpKmerFileStr;
		getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
		int tmpTabLoc = tmpKmerFileStr.find("\t");
		currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
	}

	while(!mergedKmer_ifs.eof())
	{
		string tmpStr;
		getline(mergedKmer_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//cout << "tmpStr: " << tmpStr << endl;
		int tabLoc = tmpStr.find("\t");
		string tmpMergedStr = tmpStr.substr(0, tabLoc);
		//cout << "tmpMergedStr: " << tmpMergedStr << endl;
		vector<bool> tmp_kmer_exist_in_sample_bool_vec;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		int tmpKmerExistingClassFlag = 0;
		for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp ++)
		{
			if(tmpMergedStr == currentKmerStrVec[tmp])
			{
				tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
				tmpKmerExistingClassFlag += pow(2, tmp);
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

		//cout << "tmpKmerExistingClassFlag: " << tmpKmerExistingClassFlag << endl;
		if(tmpKmerExistingClassFlag > 0)
			kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistingClassFlag << endl;
	}

	mergedKmer_ifs.close();
	kmerClass_ofs.close();
	for(int tmp = 0; tmp < providedKmerReadableFileNum; tmp++)
		*(KmerReadableIfsVec[tmp]).close();
	return 0;
}