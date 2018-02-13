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
#include "../general/fa2jf_info.h"
#include "../general/jf2Kmer_info.h"
#include "../general/Kmer2groupedKmer_info.h"
#include "../general/KmerVec2KmerOccurrence_info.h"

using namespace std;

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

string getMinCovStr(string& filePath)
{
	long long fileSize = filesize(filePath.c_str());
	long long fileSize_M = fileSize/1048576;
	string tmp_mincov_str;
	if(fileSize_M <= 300)
		tmp_mincov_str = "1";
	else if(fileSize_M <= 500)
		tmp_mincov_str = "3";
	else if(fileSize_M <= 1024)
		tmp_mincov_str = "10";
	else if(fileSize_M <= 3072)
		tmp_mincov_str = "20";
	else
		tmp_mincov_str = "50";
	return tmp_mincov_str;
}

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 jellyFishiBin" << endl;
		cout << "#2 threads_num" << endl;
		cout << "#3 Kmer_length" << endl;
		cout << "#4 bf_size" << endl;
		cout << "#5 outputDir" << endl;
		cout << "#6 withToAddIdOrNot_bool" << endl;
		cout << "#7 FaFile_num" << endl;
		cout << "#8 FaFileListFile" << endl;
		exit(1);
	}

	string jellyFishiBin = argv[1];
	string threads_num_str = argv[2];
	string Kmer_length_str = argv[3];
	string bf_size_str = argv[4];
	string outputDir = argv[5];
	outputDir += "/";

	bool withToAddIdOrNot_bool;
	string withToAddIdOrNot_bool_str = argv[6];
	if(withToAddIdOrNot_bool_str == "Y")
		withToAddIdOrNot_bool = true;
	else if(withToAddIdOrNot_bool_str == "N")
		withToAddIdOrNot_bool = false;
	else
	{
		cout << "error! withToAddIdOrNot_bool_str should be Y or N" << endl;
		exit(1); 
	}

	string FaFile_num_str = argv[7];
	int FaFile_num = atoi(FaFile_num_str.c_str());
	vector<string> FaFileVec;
	vector<string> toAddIdVec;
	string FaFileListFile = argv[8];
	ifstream FaFileList_ifs(FaFileListFile.c_str());
	int tmpToAddId = 0;
	while(!FaFileList_ifs.eof())
	{
		string tmpStr;
		getline(FaFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(!withToAddIdOrNot_bool)
		{	
			FaFileVec.push_back(tmpStr);
			toAddIdVec.push_back(int_to_str(tmpToAddId));
		}
		else
		{
			int tabLoc = tmpStr.find("\t");
			string tmpFaFile = tmpStr.substr(0, tabLoc);
			FaFileVec.push_back(tmpFaFile);
			toAddIdVec.push_back(tmpStr.substr(tabLoc + 1));
		}
		tmpToAddId ++;
	}	
	FaFileList_ifs.close();
	if(FaFile_num != FaFileVec.size()){
		cout << "error! (FaFile_num != FaFileVec.size())" << endl;
		exit(1);
	}

	cout << "start to initiate dir" << endl;
	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());
	string outputDir_jf = outputDir + "jf/";
	string outputDir_Kmer = outputDir + "Kmer/";
	string cmd_mkdir_jf = "mkdir " + outputDir_jf;
	string cmd_mkdir_Kmer = "mkdir " + outputDir_Kmer;
	system(cmd_mkdir_jf.c_str());
	system(cmd_mkdir_Kmer.c_str());

	cout << "start to do faVec 2 jfVec and jfVec 2 KmerVec" << endl;
	// faVec 2 jfVec
	for(int tmp = 0; tmp < FaFile_num; tmp++)
	{
		cout << endl << "tmpFileIndex: " << tmp << endl;
		string tmpInputFaFile = FaFileVec[tmp];
		cout << "tmp_fa_file: " << tmpInputFaFile << endl;
		string tmp_count_min_str = getMinCovStr(tmpInputFaFile);
		cout << "tmp_count_min: " << tmp_count_min_str << endl;
		string tmpOutputJfFile = outputDir_jf + toAddIdVec[tmp] + ".jf";
		cout << "tmp_jf_file: " << tmpOutputJfFile << endl; 
		Fa2jf_Info tmpFa2jfInfo;
		tmpFa2jfInfo.initaite_withMinCount(tmpInputFaFile, tmpOutputJfFile, jellyFishiBin,
			threads_num_str, tmp_count_min_str, Kmer_length_str, bf_size_str);
		tmpFa2jfInfo.fa2jf();
		string tmpInputJfFile = outputDir_jf + toAddIdVec[tmp] + ".jf";
		string tmpOutputKmerFile = outputDir_Kmer + toAddIdVec[tmp] + ".Kmer";
		cout << "tmp_kmer_file: " << tmpOutputKmerFile << endl;
		Jf2Kmer_Info jf2KmerInfo;
		jf2KmerInfo.initiate(tmpInputJfFile, tmpOutputKmerFile, jellyFishiBin);
		jf2KmerInfo.jf2Kmer();
	}
	cout << "All jobs done!" << endl;
	return 0;
}