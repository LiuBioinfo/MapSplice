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
	if(argc != 7)
	{
		cout << "Executable inputJellyFishBinPath inputIndexPath inputGenomeSubSeqDirPath outputFile outputIntermediateDir thread_num" << endl;
		exit(1);
	}
	cout << "produceGenomeLocalRegionSpecificKmerProfile starts ......" << endl;
	cout << "start to initiate indexInfo" << endl;
	// string AllowAlienOrNot_str = argv[7];
	// bool AllowAlienOrNot_bool;
	// if(AllowAlienOrNot_str == "Y")
	// 	AllowAlienOrNot_bool = true;
	// else if(AllowAlienOrNot_str == "N")
	// 	AllowAlienOrNot_bool = false;
	// else
	// {
	// 	cout << "invalid AllowAlienOrNot option, should be Y or N" << endl;
	// 	exit(1);
	// }

	string thread_num_str = argv[6];
	int thread_num = atoi(thread_num_str.c_str());
	string inputJellyFishBinPath = argv[1];
	string indexFolderPath = argv[2];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();

	cout << "start to generate jf kmer kmer_sort files for wholeGenomeSubSeqRead" << endl;
	string inputGenomeSubSeqDirPath = argv[3];
	inputGenomeSubSeqDirPath += "/";
	string wholeGenomeSubSeqReadFile = inputGenomeSubSeqDirPath + "wholeGenomeSubSeqRead.fa";

	string outputIntermediateDir = argv[5];
	outputIntermediateDir += "/";
	string outputIntermediateDir_sort = outputIntermediateDir + "tmpSortDir/";
	string mkdir_outputIntermediateDir_sort = "mkdir " + outputIntermediateDir_sort;
	system(mkdir_outputIntermediateDir_sort.c_str());
	string wholeGenomeSubSeqReadFile_jf = outputIntermediateDir + "wholeGenomeSubSeqRead.jf";
	string wholeGenomeSubSeqReadFile_kmer = outputIntermediateDir + "wholeGenomeSubSeqRead.kmer";
	string wholeGenomeSubSeqReadFile_kmer_sorted = outputIntermediateDir + "wholeGenomeSubSeqReadFile_kmer.sorted";
	
	// manipulate wholeGenomeSubSeqRead file
	string cmd_jf_count_wholeGenomeSubSeqRead = inputJellyFishBinPath + "/jellyfish count -o " 
		+ wholeGenomeSubSeqReadFile_jf + " -m 20 -t " + thread_num_str + " -s 3000000000 -C " + wholeGenomeSubSeqReadFile;
	string cmd_jf_dump_wholeGenomeSubSeqRead = inputJellyFishBinPath + "/jellyfish dump -t -c -o " 
		+ wholeGenomeSubSeqReadFile_kmer + " " + wholeGenomeSubSeqReadFile_jf;
	string cmd_sort_wholeGenomeSubSeqRead_kmer = "sort -k1 -T " + outputIntermediateDir_sort + " "
		+ wholeGenomeSubSeqReadFile_kmer + " > " + wholeGenomeSubSeqReadFile_kmer_sorted;
	// jf commands
	cout << "start to generate jf file for wholeGenomeSubSeqRead" << endl;
	system(cmd_jf_count_wholeGenomeSubSeqRead.c_str());
	cout << "start to dump jf file for wholeGenomeSubSeqRead" << endl;
	system(cmd_jf_dump_wholeGenomeSubSeqRead.c_str());
	cout << "start to sort wholeGenomeSubSeqRead kmers" << endl; 
	system(cmd_sort_wholeGenomeSubSeqRead_kmer.c_str());

	// manipulate each local region file
	cout << "start to generate jf kmer kmer_sort files for all local regions" << endl;
	string localRegionSubSeqFolder = inputGenomeSubSeqDirPath + "localRegionSubSeq/";
	string localRegionSubSeqFolder_kmer = outputIntermediateDir + "localRegionSubSeq/";
	string mkdir_localRegionSubSeqFolder_kmer = "mkdir " + localRegionSubSeqFolder_kmer;
	system(mkdir_localRegionSubSeqFolder_kmer.c_str());
	int secondLevelIndexNormalSize = indexInfo->returnSecondLevelIndexNormalSize();
	cout << "localRegionSize: " << secondLevelIndexNormalSize << endl;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		cout << "tmpChr: " << tmpChr << endl;
 		int tmpChr_subRegion_num = indexInfo->returnSecondLevelIndexPartsNum(tmpChr);
 		cout << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;
		for(int tmpSubRegion = 0; tmpSubRegion < tmpChr_subRegion_num; tmpSubRegion ++)
		{
			cout << "tmpSubRegion: " << tmpSubRegion << endl;
			string tmpSubRegionReadFile = localRegionSubSeqFolder + indexInfo->returnChrNameStr(tmpChr)
				+ "_" + int_to_str(tmpSubRegion + 1) + ".fa";
			string tmpSubRegionReadFile_jf = localRegionSubSeqFolder_kmer + indexInfo->returnChrNameStr(tmpChr)
				+ "_" + int_to_str(tmpSubRegion + 1) + ".jf";
			string tmpSubRegionReadFile_kmer = localRegionSubSeqFolder_kmer + indexInfo->returnChrNameStr(tmpChr)
				+ "_" + int_to_str(tmpSubRegion + 1) + ".kmer";
			string tmpSubRegionReadFile_kmer_sorted = localRegionSubSeqFolder_kmer + indexInfo->returnChrNameStr(tmpChr)
				+ "_" + int_to_str(tmpSubRegion + 1) + ".kmer.sorted";
			string cmd_jf_count_subRegion = inputJellyFishBinPath + "/jellyfish count -o "
				+ tmpSubRegionReadFile_jf + " -m 20 -t " + thread_num_str + " -s 3000000000 -C " + tmpSubRegionReadFile;
			string cmd_jf_dump_subRegion = inputJellyFishBinPath + "/jellyfish dump -t -c -o "
				+ tmpSubRegionReadFile_kmer + " " + tmpSubRegionReadFile_jf;
			string cmd_sort_subRegion = "sort -k1 -T " + outputIntermediateDir_sort + " "
				+ tmpSubRegionReadFile_kmer + " > " + tmpSubRegionReadFile_kmer_sorted;
			cout << "start to do jf count for tmp local region" << endl;
			system(cmd_jf_count_subRegion.c_str());
			cout << "start to do dump for tmp local region" << endl;
			system(cmd_jf_dump_subRegion.c_str());
			cout << "start to do sorting for tmp local region" << endl;
			system(cmd_sort_subRegion.c_str());
		}
	}

	// generate local region specific kmer
	cout << "start to generate local region specific kmers" << endl;
	cout << "start to initiate localRegionKmerIfsVec" << endl;
	vector<ifstream*> KmerReadableIfsVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		int tmpChr_subRegion_num = indexInfo->returnSecondLevelIndexPartsNum(tmpChr);
		for(int tmpSubRegion = 0; tmpSubRegion < tmpChr_subRegion_num; tmpSubRegion ++)
		{
			string tmpSubRegionReadFile_kmer_sorted = localRegionSubSeqFolder_kmer 
				+ indexInfo->returnChrNameStr(tmpChr) + "_" + int_to_str(tmpSubRegion + 1) + ".kmer.sorted";
			ifstream *tmpSubRegionKmer_ifs = new ifstream(tmpSubRegionReadFile_kmer_sorted.c_str());
			KmerReadableIfsVec.push_back(tmpSubRegionKmer_ifs);
		}
	}
	
	vector<bool> kmer_file_end_bool_vec;
	for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp ++)
		kmer_file_end_bool_vec.push_back(false);
	vector<string> currentKmerStrVec;
	for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp ++)
		currentKmerStrVec.push_back("");
	//cout << "head lines" << endl;
	for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp ++)
	{
		string tmpKmerFileStr;
		getline((*KmerReadableIfsVec[tmp]), tmpKmerFileStr);
		//cout << "tmpKmerFileStr: " << tmpKmerFileStr << endl;
		int tmpTabLoc = tmpKmerFileStr.find("\t");
		currentKmerStrVec[tmp] = tmpKmerFileStr.substr(0, tmpTabLoc);
		//cout << "tmpCurrentKmerStr: " << tmp << endl << currentKmerStrVec[tmp] << endl;
	}

	string outputFile = argv[4];
	ifstream mergedKmer_ifs(wholeGenomeSubSeqReadFile_kmer_sorted.c_str());
	ofstream kmerClass_ofs(outputFile.c_str());
	while(!mergedKmer_ifs.eof())
	{
		string tmpStr;
		getline(mergedKmer_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpMergedStr = tmpStr.substr(0, tabLoc);
		int tmpKmerExistence_fileNum = 0;
		int tmpKmerExistence_lastFileIndex = -1;
		vector<bool> tmp_kmer_exist_in_sample_bool_vec;
		for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp ++)
			tmp_kmer_exist_in_sample_bool_vec.push_back(false);

		//int tmpKmerExistingClassFlag = 0;
		for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp ++)
		{
			if(tmpMergedStr == currentKmerStrVec[tmp])
			{
				tmp_kmer_exist_in_sample_bool_vec[tmp] = true;
				tmpKmerExistence_fileNum ++;
				tmpKmerExistence_lastFileIndex = tmp;
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
		if(tmpKmerExistence_fileNum == 0) // crossing local region boundaries
			kmerClass_ofs << tmpMergedStr << "\t" << secondLevelIndexNormalSize << endl;
		else if(tmpKmerExistence_fileNum == 1) // local region specific kmers
			kmerClass_ofs << tmpMergedStr << "\t" << tmpKmerExistence_lastFileIndex << endl;
		else if(tmpKmerExistence_fileNum > 1) // multi local region kmers
			kmerClass_ofs << tmpMergedStr << "\t" << secondLevelIndexNormalSize << endl;
		else
		{
			cout << "error ! tmpKmerExistence_fileNum < 0 !" << endl;
			exit(1); 
		}
	}

	for(int tmp = 0; tmp < secondLevelIndexNormalSize; tmp++)
	{
		(*KmerReadableIfsVec[tmp]).close();
		delete KmerReadableIfsVec[tmp];
	}
	mergedKmer_ifs.close();
	kmerClass_ofs.close();

	delete indexInfo;
	return 0;
}