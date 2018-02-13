// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef SORTEDKMERVEC2KMEROCCURRENCEPROFILE_INFO_H
#define SORTEDKMERVEC2KMEROCCURRENCEPROFILE_INFO_H
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

class SortedKmerVec2KmerOccurrenceProfile_Info
{
private:
	// input
	string inputSortedKmerFile_union;
	vector<string> inputSortedKmerFileVec;
	vector<string> input2AddIdVec;
	string temporaryDir;
	// output
	string outputKmerOccurrenceProfile;
public:
	SortedKmerVec2KmerOccurrenceProfile_Info()
	{}

	void initiate_noMkTmpDir(string& tmp_inputSortedKmerFile_union, vector<string>& tmp_inputSortedKmerFileVec, 
		vector<string>& tmp_input2AddIdVec, string& tmp_tmpDir, string& tmp_outputKmerOccurrenceProfile)
	{
		// input
		inputSortedKmerFile_union = tmp_inputSortedKmerFile_union;
		for(int tmp = 0; tmp < tmp_inputSortedKmerFileVec.size(); tmp++)
		{	
			inputSortedKmerFileVec.push_back(tmp_inputSortedKmerFileVec[tmp]);
			input2AddIdVec.push_back(tmp_input2AddIdVec[tmp]);
		}
		temporaryDir = tmp_tmpDir;
		temporaryDir += "/";
		// output
		outputKmerOccurrenceProfile = tmp_outputKmerOccurrenceProfile;
	}	

	void initiate_noMkTmpDir(string& tmp_inputSortedKmerFile_union, vector<string>& tmp_inputSortedKmerFileVec,
		string& tmp_tmpDir, string& tmp_outputKmerOccurrenceProfile)
	{
		vector<string> tmp_input2AddIdVec;
		for(int tmp = 0; tmp < tmp_inputSortedKmerFileVec.size(); tmp++)
			tmp_input2AddIdVec.push_back(int_to_str(tmp));
		this->initiate_noMkTmpDir(tmp_inputSortedKmerFile_union, tmp_inputSortedKmerFileVec,
			tmp_input2AddIdVec, tmp_tmpDir, tmp_outputKmerOccurrenceProfile);
	}

	void initiate_mkTmpDir(string& tmp_inputSortedKmerFile_union, vector<string>& tmp_inputSortedKmerFileVec, 
		vector<string>& tmp_input2AddIdVec, string& tmp_tmpDir, string& tmp_outputKmerOccurrenceProfile)
	{
		this->initiate_noMkTmpDir(tmp_inputSortedKmerFile_union,
			tmp_inputSortedKmerFileVec, tmp_input2AddIdVec,
			tmp_tmpDir, tmp_outputKmerOccurrenceProfile);
		string cmd_mkdir_temporaryDir = "mkdir " + temporaryDir;
		system(cmd_mkdir_temporaryDir.c_str());
	}

	void initiate_mkTmpDir(string& tmp_inputSortedKmerFile_union, vector<string>& tmp_inputSortedKmerFileVec,
		string& tmp_tmpDir, string& tmp_outputKmerOccurrenceProfile)
	{
		vector<string> tmp_input2AddIdVec;
		for(int tmp = 0; tmp < tmp_inputSortedKmerFileVec.size(); tmp++)
			tmp_input2AddIdVec.push_back(int_to_str(tmp));
		this->initiate_mkTmpDir(tmp_inputSortedKmerFile_union, tmp_inputSortedKmerFileVec,
			tmp_input2AddIdVec, tmp_tmpDir, tmp_outputKmerOccurrenceProfile);
	}	

	void generateKmerOccurrenceProfile()
	{
		for(int tmp = 0; tmp < inputSortedKmerFileVec.size(); tmp++)
			this->getIntermediateKmerOccurrenceProfile(tmp);
	}

	void getIntermediateKmerOccurrenceProfile(int tmpSortedKmerFileIndexInVec) // cmp union file 2 each kmer file, then generate single Kmer occurrence profile
	{
		//string toAddFileIndex = tmpSortedKmerFileIndexInVec;
		string tmpSortedKmerFile = inputSortedKmerFileVec[tmpSortedKmerFileIndexInVec];
		string tmpInputSortedKmerOccurrenceFile_ori;
		string tmpOutputSortedKmerOccurrenceFile_updated;

		if(tmpSortedKmerFileIndexInVec == 0)
			tmpInputSortedKmerOccurrenceFile_ori = inputSortedKmerFile_union;
		else
			tmpInputSortedKmerOccurrenceFile_ori = temporaryDir + "/KmerOccurrence." + int_to_str(tmpSortedKmerFileIndexInVec - 1);
		
		if(tmpSortedKmerFileIndexInVec != inputSortedKmerFileVec.size() - 1)
			tmpOutputSortedKmerOccurrenceFile_updated = temporaryDir + "/KmerOccurrence." + int_to_str(tmpSortedKmerFileIndexInVec);
		else // tmpSortedKmerFileIndexInVec == inputSortedKmerFileVec.size() - 1
			tmpOutputSortedKmerOccurrenceFile_updated = outputKmerOccurrenceProfile;
	
		ifstream sortedKmer_ifs(tmpSortedKmerFile.c_str());
		ifstream occurrence_ori_ifs(tmpInputSortedKmerOccurrenceFile_ori.c_str());
		ofstream occurrence_updated_ofs(tmpOutputSortedKmerOccurrenceFile_updated.c_str());
		
		bool sortedKmer_file_end_bool = false;
		
		string tmp1stSortedKmerFileStr;
		getline(sortedKmer_ifs, tmp1stSortedKmerFileStr);
		int currentKmerTabLoc = tmp1stSortedKmerFileStr.find("\t");
		string currentKmerStr = tmp1stSortedKmerFileStr.substr(0, currentKmerTabLoc);		
		
		while(!occurrence_ori_ifs.eof())
		{
			string tmpStr;
			getline(occurrence_ori_ifs, tmpStr);
			if(tmpStr )
		}
		sortedKmer_ifs.close();
		occurrence_ori_ifs.close();
		occurrence_updated_ofs.close();
	}

	void updateUnionKmerOccurrenceProfile(string& inputOriUnionKmerOccurrenceProfileFile,
		string& outputUpdatedUnionKmerOccurrenceProfileFile, string& toAddSampleId)
	{
		ifstream ori_ifs(inputOriUnionKmerOccurrenceProfileFile.c_str());
		ofstream updated_ofs(outputUpdatedUnionKmerOccurrenceProfileFile.c_str());
		while(!ori_ifs.eof())
		{
			string tmpOriStr;
			getline(ori_ifs, tmpOriStr);
			if(tmpOriStr == "")
				break;
		}
		ori_ifs.close();
		updated_ofs.close();
	}

	void rm_tmpDir()
	{
		string cmd_rm_tmpDir = "rm -r " + temporaryDir;
		system(cmd_rm_tmpDir.c_str());
	}

	void printoutSortedKmerVec2KmerOccurrenceProfileInfo(string& infoFile)
	{
		ofstream info_ofs(infoFile.c_str());
		info_ofs << "Kmer_Occurrence_Profile:" << endl << outputKmerOccurrenceProfile << endl; 
		info_ofs << "Sorted_Kmer_Union:" << endl << inputSortedKmerFile_union << endl;
		info_ofs << "Sorted_Kmer_Vec:" << endl;
		for(int tmp = 0; tmp < inputSortedKmerFileVec.size(); tmp++)
			info_ofs << tmp << "\t" << inputSortedKmerFileVec[tmp] << endl;
		info_ofs.close();
	}
};
#endif