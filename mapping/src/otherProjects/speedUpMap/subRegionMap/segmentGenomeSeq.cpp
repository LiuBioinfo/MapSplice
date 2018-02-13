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
//#include <omp.h>
#include "../../../general/read_block_test.h"
#include "../../../general/bwtmap_info.h"
#include "../../../general/DoubleAnchorScore.h"
#include "../../../general/sbndm.h"
#include "../../../general/splice_info.h"
#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputIndexFolderPath outputFolder" << endl;
		exit(1);
	}
	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	log_ofs << "start to initiate indexInfo" << endl;
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;

	// string subSeqLengthStr = argv[2];
	// int subSeqLength = atoi(subSeqLengthStr.c_str());
	// log_ofs << "subSeqLength: " << subSeqLength << endl;
	// int jumpLength = subSeqLength/2;

	// string wholeGenomeSubSeqReadFile = outputFolderStr + "wholeGenomeSubSeqRead.fa";
	// ofstream wholeGenomeSubSeqRead_ofs(wholeGenomeSubSeqReadFile.c_str());
	// string localRegionSubSeqFolder = outputFolderStr + "/localRegionSubSeq/";
	// string mkdir_localRegionSubSeq = "mkdir -p " + localRegionSubSeqFolder;
	// system(mkdir_localRegionSubSeq.c_str());
	int secondLevelIndexNormalSize = indexInfo->returnSecondLevelIndexNormalSize();
	cout << "localRegionSize: " << secondLevelIndexNormalSize << endl;
	log_ofs << "localRegionSize: " << secondLevelIndexNormalSize << endl;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		int tmpChrLength = indexInfo->returnChromLength(tmpChr);
		cout << "tmpChr: " << tmpChrName << endl;
		log_ofs << "tmpChr: " << tmpChrName << endl;
		string tmpChrFolder = outputFolderStr + "/" + tmpChrName;
		string mkdir_tmpChrFolder = "mkdir -p " + tmpChrFolder;
		system(mkdir_tmpChrFolder.c_str());

		int tmpChr_subRegion_num = indexInfo->returnSecondLevelIndexPartsNum(tmpChr);
		cout << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;
		log_ofs << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;

		for(int tmpSubRegion = 0; tmpSubRegion < tmpChr_subRegion_num; tmpSubRegion ++)
		{
			string tmpChrSubRegionFolder = outputFolderStr + "/" + tmpChrName + "/" + int_to_str(tmpSubRegion + 1);
			string mkdir_tmpChrSubRegionFolder = "mkdir -p " + tmpChrSubRegionFolder;
			system(mkdir_tmpChrSubRegionFolder.c_str());	

			string tmpSubRegionSeqFile = tmpChrSubRegionFolder + "/" + tmpChrName + "_" + int_to_str(tmpSubRegion + 1) + ".fa";
			ofstream tmpSubRegionSeq_ofs(tmpSubRegionSeqFile.c_str());
			int tmpSubRegionSeqBaseNum = 0;
			for(int tmpBase = 0; tmpBase < secondLevelIndexNormalSize; tmpBase ++)
			{
				int tmpBasePosInChr = tmpBase + secondLevelIndexNormalSize * tmpSubRegion + 1;
				if(tmpBasePosInChr <= tmpChrLength)
				{	
					char tmpBaseChar = indexInfo->returnOneBaseCharInGenome(tmpChr, tmpBasePosInChr);
					tmpSubRegionSeqBaseNum ++;
					if((tmpSubRegionSeqBaseNum/50) * 50 == tmpSubRegionSeqBaseNum)
						tmpSubRegionSeq_ofs << tmpBaseChar << endl;
					else
						tmpSubRegionSeq_ofs << tmpBaseChar;
				}
			}
			tmpSubRegionSeq_ofs << endl;
			tmpSubRegionSeq_ofs.close();
		}
	}
	cout << "End of job !" << endl;
	log_ofs << "End of job !" << endl;
	//wholeGenomeSubSeqRead_ofs.close();
	log_ofs.close();
	return 0;
}