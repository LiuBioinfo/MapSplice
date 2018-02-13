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
#include "../general/Kmer2groupedKmer_info.h"
#include "../general/KmerVec2KmerOccurrence_info.h"
//#include "../general/KmerVec2setSpecificKmerProfile_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << endl << "Executable#0 total_Kmer_length#1 grouped_Kmer_prefix_length#2" << endl;
		cout << "outputFolder#3 withToAddIdOrNot_bool#4 KmerFileNum#5" << endl;
		cout << "outputSolidKmerCountOrNot_bool#6 KmerFileListFile#7" << endl;
		exit(1);
	}
	string total_Kmer_length_str = argv[1];
	string grouped_Kmer_prefix_length_str = argv[2];
	int total_Kmer_length = atoi(total_Kmer_length_str.c_str());
	int grouped_Kmer_prefix_length = atoi(grouped_Kmer_prefix_length_str.c_str());
	int toKeepCount_Kmer_length = total_Kmer_length - grouped_Kmer_prefix_length;
	int binNum = pow(4, grouped_Kmer_prefix_length);
	unsigned long long int toKeepCount_bin_Kmer_num = pow(4, toKeepCount_Kmer_length);
	if((total_Kmer_length != 20)||(grouped_Kmer_prefix_length != 3)||(toKeepCount_Kmer_length != 17)
		||(binNum != 64)||(toKeepCount_bin_Kmer_num != 17179869184))
	{
		cout << "error! (total_Kmer_length != 20)||(grouped_Kmer_prefix_length != 3)||(toKeepCount_Kmer_length != 17)" << endl;
		cout << "error! (binNum != 64)||(toKeepCount_bin_Kmer_num != 17179869184)" << endl;
		cout << "total_Kmer_length: " << total_Kmer_length << endl;
		cout << "grouped_Kmer_prefix_length: " << grouped_Kmer_prefix_length << endl;
		cout << "binNum: " << binNum << endl;
		cout << "toKeepCount_Kmer_length: " << toKeepCount_Kmer_length << endl;
		cout << "toKeepCount_bin_Kmer_num: " << toKeepCount_bin_Kmer_num << endl;
		exit(1);
	}

	bool withToAddIdOrNot_bool;
	string withToAddIdOrNot_bool_str = argv[4];
	if(withToAddIdOrNot_bool_str == "Y")
		withToAddIdOrNot_bool = true;
	else if(withToAddIdOrNot_bool_str == "N")
		withToAddIdOrNot_bool = false;
	else
	{
		cout << "error! withToAddIdOrNot_bool_str should be Y or N" << endl;
		exit(1); 
	}

	bool outputSolidKmerCountOrNot_bool;
	string outputSolidKmerCountOrNot_bool_str = argv[6];
	if(outputSolidKmerCountOrNot_bool_str == "Y")
		outputSolidKmerCountOrNot_bool = true;
	else if(outputSolidKmerCountOrNot_bool_str == "N")
		outputSolidKmerCountOrNot_bool = false;
	else
	{
		cout << "error! outputSolidKmerCountOrNot_bool_str should be Y or N" << endl;
		exit(1); 
	}

	string KmerFile_num_str = argv[5];
	int KmerFile_num = atoi(KmerFile_num_str.c_str());
	vector<string> KmerFileVec;
	vector<string> toAddIdVec;
	string KmerFileListFile = argv[7];
	ifstream KmerFileList_ifs(KmerFileListFile.c_str());
	int tmpToAddId = 0;
	while(!KmerFileList_ifs.eof())
	{
		string tmpStr;
		getline(KmerFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(!withToAddIdOrNot_bool)
		{	
			KmerFileVec.push_back(tmpStr);
			toAddIdVec.push_back(int_to_str(tmpToAddId));
		}
		else
		{
			int tabLoc = tmpStr.find("\t");
			string tmpKmerFile = tmpStr.substr(0, tabLoc);
			KmerFileVec.push_back(tmpKmerFile);
			toAddIdVec.push_back(tmpStr.substr(tabLoc + 1));
		}
		tmpToAddId ++;
	}
	KmerFileList_ifs.close();
	if(KmerFile_num != KmerFileVec.size()){
		cout << "error! (KmerFile_num != KmerFileVec.size())" << endl;
		exit(1);
	}

	cout << "mkdir output" << endl;
	string outputFolder = argv[3]; outputFolder += "/";
	string cmd_mkdir = "mkdir " + outputFolder;
	system(cmd_mkdir.c_str());

	cout << "mkdir outputFolder_prefixGroupedKmer" << endl;
	string outputFolder_prefixGroupedKmer = outputFolder + "prefixGroupedKmer/";
	string cmd_mkdir_prefixGroupedKmer = "mkdir " + outputFolder_prefixGroupedKmer;
	system(cmd_mkdir_prefixGroupedKmer.c_str());

	cout << "start to divide raw kmer files" << endl;
	vector<Kmer2groupedKmer_Info> Kmer2groupedKmerInfoVec;
	for(int tmp = 0; tmp < KmerFile_num; tmp++)
	{
		string tmpKmerFile = KmerFileVec[tmp];
		string tmpToAddId = toAddIdVec[tmp];
		string tmpDir = outputFolder_prefixGroupedKmer + tmpToAddId + "/";
		Kmer2groupedKmer_Info tmpKmer2groupedKmerInfo;
		tmpKmer2groupedKmerInfo.initiate_mkdir(tmpKmerFile, tmpDir, grouped_Kmer_prefix_length);
		tmpKmer2groupedKmerInfo.groupKmer();
		Kmer2groupedKmerInfoVec.push_back(tmpKmer2groupedKmerInfo);
	}

	cout << "start to get KmerVec2KmerOccurrenceNum_info" << endl;
	string outputFolder_KmerOccurrenceNum = outputFolder + "occurrenceNum/";
	string cmd_mkdir_KmerOccurrenceNum = "mkdir " + outputFolder_KmerOccurrenceNum;
	system(cmd_mkdir_KmerOccurrenceNum.c_str());
	string outputFolder_setSpecificKmer = outputFolder + "setSpecificKmer/";
	string cmd_mkdir_setSpecificKmer = "mkdir " + outputFolder_setSpecificKmer;
	system(cmd_mkdir_setSpecificKmer.c_str());
	for(int tmpBin = 0; tmpBin < binNum; tmpBin ++)
	{
		string tmpSetSpecificKmerDir_bin = outputFolder_setSpecificKmer + "bin_" + int_to_str(tmpBin);
		string cmd_mkdir_tmpSetSpecificKmerDir_bin = "mkdir " + tmpSetSpecificKmerDir_bin;
		system(cmd_mkdir_tmpSetSpecificKmerDir_bin.c_str()); 
	}

	cout << "start to initiate KmerVec2KmerOccurrenceInfo" << endl;
	KmerVec2KmerOccurrence_Info KmerVec2KmerOccurrenceInfo;
	KmerVec2KmerOccurrenceInfo.initiate(total_Kmer_length, grouped_Kmer_prefix_length);

	/////////////////////  TEST  //////////////////////////
	/////////////////////  TEST  //////////////////////////
	/////////////////////  TEST  //////////////////////////
	// string tmpStr = "AAAAAAAAAAAAAAAAAAAA";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "AAAAAAAAAAAAAAAAAAAC";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "ACTAAAAAAAAAAAAAAAAG";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "AAAAAAAAAAAAAAAAAAAT";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "AGAAAAAAAAAAAAAAAACA";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "AAAAAAAAAAAAAAAAAACC";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "ACCAAAAAAAAAAAAAAACG";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;
	// tmpStr = "ATTAAAAAAAAAAAAAAACT";
	// cout << "tmpStr:" << endl << tmpStr << endl << "tmpKmer2int:" << KmerVec2KmerOccurrenceInfo.toKeepCount_Kmer2int(tmpStr) << endl;		

	// unsigned long long tmpInt = 0;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 10) << endl;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 20) << endl; 
	// tmpInt = 1;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 10) << endl;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 20) << endl; 
	// tmpInt = 3;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 10) << endl;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 20) << endl;
	// tmpInt = 4;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 10) << endl;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 20) << endl;
	// tmpInt = 12;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 10) << endl;
	// cout << "tmpInt: " << endl << tmpInt << endl << "tmpInt2Kmer:" << KmerVec2KmerOccurrenceInfo.int2Kmer(tmpInt, 20) << endl;
	// exit(1);
	/////////////////////  TEST  //////////////////////////
	/////////////////////  TEST  //////////////////////////
	/////////////////////  TEST  //////////////////////////

	char* toKeepKmerCount_charArray = new char[toKeepCount_bin_Kmer_num];
	for(unsigned long long int tmp = 0; tmp < toKeepCount_bin_Kmer_num; tmp++)
		toKeepKmerCount_charArray[tmp] = 0;

	cout << "binNum: " << binNum << endl;
	for(int tmpBin = 0; tmpBin < binNum; tmpBin ++)
	{	
		cout << "tmpBin: " << tmpBin << endl;
	
		///////////////////////////////////////
		// start to obtain Kmer occurrence num
		///////////////////////////////////////
		cout << "start to get Kmer occurrenceNum" << endl;		
		for(int tmpKmerFileIndex = 0; tmpKmerFileIndex < KmerFile_num; tmpKmerFileIndex++)
		{
			cout << "tmp updateKmerOccurrenceNum, tmpKmerFileIndex: " << tmpKmerFileIndex << endl;
			string tmpGroupedKmerFile = Kmer2groupedKmerInfoVec[tmpKmerFileIndex].returnTmpBinFile(tmpBin);
			cout << "tmpGroupedKmerFile: " << tmpGroupedKmerFile << endl;
			KmerVec2KmerOccurrenceInfo.updateKmerOccurrenceNum(toKeepKmerCount_charArray, tmpGroupedKmerFile);
		}
		//////////////////////////////////////
		// start to print Kmer occurrence num
		//////////////////////////////////////
		if(outputSolidKmerCountOrNot_bool)
		{	
			cout << "start to print Kmer occurrence num" << endl;
			string tmpBinOccurrenceNumFile = outputFolder_KmerOccurrenceNum + "bin." + int_to_str(tmpBin) + ".solid.Kmer";
			KmerVec2KmerOccurrenceInfo.printoutKmerOccurrenceNum(toKeepKmerCount_charArray, 
				tmpBinOccurrenceNumFile, tmpBin, false);
		}
		//////////////////////////////////////////
		// start to print set specific Kmer profile
		///////////////////////////////////////////
		cout << "start to print set specific kmer profile" << endl;
		for(int tmpKmerFileIndex = 0; tmpKmerFileIndex < KmerFile_num; tmpKmerFileIndex++)
		{
			cout << "tmp printoutSetSpecificKmer, tmpKmerFileIndex: " << tmpKmerFileIndex << endl;
			string tmpGroupedKmerFile = Kmer2groupedKmerInfoVec[tmpKmerFileIndex].returnTmpBinFile(tmpBin);
			cout << "tmpGroupedKmerFile: " << tmpGroupedKmerFile << endl;
			string tmpOutputSetSpecificKmerFile = outputFolder_setSpecificKmer + "bin_" + int_to_str(tmpBin)
				+ "/" + toAddIdVec[tmpKmerFileIndex] + ".setSpecific.Kmer";
			cout << "tmpOutputSetSpecificKmerFile: " << tmpOutputSetSpecificKmerFile << endl;
			KmerVec2KmerOccurrenceInfo.printoutSetSpecificKmer(toKeepKmerCount_charArray, 
				toAddIdVec[tmpKmerFileIndex], tmpGroupedKmerFile, tmpOutputSetSpecificKmerFile);
		}
		memset(toKeepKmerCount_charArray, 0, toKeepCount_bin_Kmer_num * sizeof(char));
		//for(unsigned long long int tmp = 0; tmp < toKeepCount_bin_Kmer_num; tmp++)
		//	toKeepKmerCount_charArray[tmp] = 0;		
	}
	delete []toKeepKmerCount_charArray;
	return 0;
}