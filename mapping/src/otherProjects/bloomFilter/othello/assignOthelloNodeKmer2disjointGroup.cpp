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
	if(argc != 11)
	{
		cout << "Executable inputMergedKmerReadableFile outputMergedKmerClassInfoFile KmerReadableFile_1 KmerReadableFile_2";
		cout << " KmerReadableFile_3 KmerReadableFile_4 KmerReadableFile_5 KmerReadableFile_6 KmerReadableFile_7 KmerReadableFile_8" << endl;
		exit(1);
	}
	int kmerReadableFileNum = 8;

	string inputMergedKmerReadableFile = argv[1];
	string outputMergedKmerClassInfoFile = argv[2];
	vector<string> kmerReadableFileVec;
	for(int tmp = 0; tmp < kmerReadableFileNum; tmp ++)
		kmerReadableFileVec.push_back(argv[3 + tmp]);

	ofstream kmerClass_ofs(outputMergedKmerClassInfoFile.c_str());
	ifstream mergedKmer_ifs(inputMergedKmerReadableFile.c_str());
	ifstream kmer_0_ifs(kmerReadableFileVec[0].c_str());
	ifstream kmer_1_ifs(kmerReadableFileVec[1].c_str());
	ifstream kmer_2_ifs(kmerReadableFileVec[2].c_str());
	ifstream kmer_3_ifs(kmerReadableFileVec[3].c_str());
	ifstream kmer_4_ifs(kmerReadableFileVec[4].c_str());
	ifstream kmer_5_ifs(kmerReadableFileVec[5].c_str());
	ifstream kmer_6_ifs(kmerReadableFileVec[6].c_str());
	ifstream kmer_7_ifs(kmerReadableFileVec[7].c_str());	

	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < kmerReadableFileNum; tmp ++)
		kmer_file_end_bool_vec.push_back(false);

	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < kmerReadableFileNum; tmp ++)
		currentKmerStrVec.push_back("");

	int tmpTabLoc;
	string tmpKmerFileStr;
	getline(kmer_0_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[0] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_1_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[1] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_2_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[2] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_3_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[3] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_4_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[4] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_5_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[5] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_6_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[6] = tmpKmerFileStr.substr(0, tmpTabLoc);
	getline(kmer_7_ifs, tmpKmerFileStr);
	tmpTabLoc = tmpKmerFileStr.find("\t");
	currentKmerStrVec[7] = tmpKmerFileStr.substr(0, tmpTabLoc);

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
		for(int tmp = 0; tmp < kmerReadableFileNum; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		int tmpKmerExistingClassFlag = 0;
		// file 1
		//cout << "currentKmerStrVec[0]: " << currentKmerStrVec[0] << endl;
		if(tmpMergedStr == currentKmerStrVec[0])
		{
			tmp_kmer_exist_in_sample_bool_vec[0] = true;
			tmpKmerExistingClassFlag += 128;
			if(!kmer_file_end_bool_vec[0])
			{
				string tmpKmerFileStr;
				getline(kmer_0_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[0] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[0] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[0] = false;

		// file 2
		//cout << "currentKmerStrVec[1]: " << currentKmerStrVec[1] << endl;
		if(tmpMergedStr == currentKmerStrVec[1])
		{
			tmp_kmer_exist_in_sample_bool_vec[1] = true;
			tmpKmerExistingClassFlag += 64;
			if(!kmer_file_end_bool_vec[1])
			{
				string tmpKmerFileStr;
				getline(kmer_1_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[1] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[1] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[1] = false;

		// file 3
		//cout << "currentKmerStrVec[2]: " << currentKmerStrVec[2] << endl;
		if(tmpMergedStr == currentKmerStrVec[2])
		{
			tmp_kmer_exist_in_sample_bool_vec[2] = true;
			tmpKmerExistingClassFlag += 32;
			if(!kmer_file_end_bool_vec[2])
			{
				string tmpKmerFileStr;
				getline(kmer_2_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[2] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[2] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[2] = false;		

		// file 4
		//cout << "currentKmerStrVec[3]: " << currentKmerStrVec[3] << endl;
		if(tmpMergedStr == currentKmerStrVec[3])
		{
			tmp_kmer_exist_in_sample_bool_vec[3] = true;
			tmpKmerExistingClassFlag += 16;
			if(!kmer_file_end_bool_vec[3])
			{
				string tmpKmerFileStr;
				getline(kmer_3_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[3] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[3] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[3] = false;

		// file 5
		//cout << "currentKmerStrVec[4]: " << currentKmerStrVec[4] << endl;
		if(tmpMergedStr == currentKmerStrVec[4])
		{
			tmp_kmer_exist_in_sample_bool_vec[4] = true;
			tmpKmerExistingClassFlag += 8;
			if(!kmer_file_end_bool_vec[4])
			{
				string tmpKmerFileStr;
				getline(kmer_4_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[4] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[4] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[4] = false;

		// file 6
		//cout << "currentKmerStrVec[5]: " << currentKmerStrVec[5] << endl;
		if(tmpMergedStr == currentKmerStrVec[5])
		{
			tmp_kmer_exist_in_sample_bool_vec[5] = true;
			tmpKmerExistingClassFlag += 4;
			if(!kmer_file_end_bool_vec[5])
			{
				string tmpKmerFileStr;
				getline(kmer_5_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[5] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[5] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[5] = false;

		// file 7
		//cout << "currentKmerStrVec[6]: " << currentKmerStrVec[6] << endl;
		if(tmpMergedStr == currentKmerStrVec[6])
		{
			tmp_kmer_exist_in_sample_bool_vec[6] = true;
			tmpKmerExistingClassFlag += 2;
			if(!kmer_file_end_bool_vec[6])
			{
				string tmpKmerFileStr;
				getline(kmer_6_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[6] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[6] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[6] = false;

		// file 8
		//cout << "currentKmerStrVec[7]: " << currentKmerStrVec[7] << endl;
		if(tmpMergedStr == currentKmerStrVec[7])
		{
			tmp_kmer_exist_in_sample_bool_vec[7] = true;
			tmpKmerExistingClassFlag += 1;
			if(!kmer_file_end_bool_vec[7])
			{
				string tmpKmerFileStr;
				getline(kmer_7_ifs, tmpKmerFileStr);
				if(tmpKmerFileStr == "")
					kmer_file_end_bool_vec[7] = true;
				else
				{
					int tmpTabLoc = tmpKmerFileStr.find("\t");
					currentKmerStrVec[7] = tmpKmerFileStr.substr(0, tmpTabLoc);
				}
			}
		}
		else
			tmp_kmer_exist_in_sample_bool_vec[7] = false;

		//cout << "tmpKmerExistingClassFlag: " << tmpKmerExistingClassFlag << endl;
		kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistingClassFlag << endl;
	}
	kmer_0_ifs.close();
	kmer_1_ifs.close();
	kmer_2_ifs.close();
	kmer_3_ifs.close();
	kmer_4_ifs.close();
	kmer_5_ifs.close();
	kmer_6_ifs.close();
	kmer_7_ifs.close();
	mergedKmer_ifs.close();
	kmerClass_ofs.close();
	return 0;
}