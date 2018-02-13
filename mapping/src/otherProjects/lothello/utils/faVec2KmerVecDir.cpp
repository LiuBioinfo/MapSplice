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

int main(int argc, char** argv)
{
	if(argc != 10)
	{
		cout << "Executable#0 jellyFishiBin#1 threads_num#2 count_min#3 Kmer_length#4 bf_size#5" << endl;
		cout << "outputDir#6 withToAddIdOrNot_bool#7" << endl;
		cout << "FaFile_num#8 FaFileListFile#9" << endl;
		exit(1);
	}

	string jellyFishiBin = argv[1];
	string threads_num_str = argv[2];
	string count_min_str = argv[3];
	string Kmer_length_str = argv[4];
	string bf_size_str = argv[5];
	string outputDir = argv[6];
	outputDir += "/";

	bool withToAddIdOrNot_bool;
	string withToAddIdOrNot_bool_str = argv[7];
	if(withToAddIdOrNot_bool_str == "Y")
		withToAddIdOrNot_bool = true;
	else if(withToAddIdOrNot_bool_str == "N")
		withToAddIdOrNot_bool = false;
	else
	{
		cout << "error! withToAddIdOrNot_bool_str should be Y or N" << endl;
		exit(1); 
	}

	string FaFile_num_str = argv[8];
	int FaFile_num = atoi(FaFile_num_str.c_str());
	vector<string> FaFileVec;
	vector<string> toAddIdVec;
	string FaFileListFile = argv[9];
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
		cout << "tmpFileIndex: " << tmp << endl;
		string tmpInputFaFile = FaFileVec[tmp];
		string tmpOutputJfFile = outputDir_jf + toAddIdVec[tmp] + ".jf"; 
		Fa2jf_Info tmpFa2jfInfo;
		tmpFa2jfInfo.initaite_withMinCount(tmpInputFaFile, tmpOutputJfFile, jellyFishiBin,
			threads_num_str, count_min_str, Kmer_length_str, bf_size_str);
		tmpFa2jfInfo.fa2jf();
		string tmpInputJfFile = outputDir_jf + toAddIdVec[tmp] + ".jf";
		string tmpOutputKmerFile = outputDir_Kmer + toAddIdVec[tmp] + ".Kmer";
		Jf2Kmer_Info jf2KmerInfo;
		jf2KmerInfo.initiate(tmpInputJfFile, tmpOutputKmerFile, jellyFishiBin);
		jf2KmerInfo.jf2Kmer();
	}
	cout << "All jobs done!" << endl;
	return 0;
}