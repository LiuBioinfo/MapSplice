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
#include "../../../../general/read_block_test.h"
#include "../../../../general/bwtmap_info.h"
#include "../../../../general/DoubleAnchorScore.h"
#include "../../../../general/sbndm.h"
#include "../../../../general/splice_info.h"
#include "../../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath subSeqLength outputFolder" << endl;
		exit(1);
	}
	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[3];
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

	string subSeqLengthStr = argv[2];
	int subSeqLength = atoi(subSeqLengthStr.c_str());
	log_ofs << "subSeqLength: " << subSeqLength << endl;
	int jumpLength = subSeqLength/2;

	string wholeGenomeSubSeqReadFile = outputFolderStr + "wholeGenomeSubSeqRead.fa";
	ofstream wholeGenomeSubSeqRead_ofs(wholeGenomeSubSeqReadFile.c_str());
	string localRegionSubSeqFolder = outputFolderStr + "/localRegionSubSeq/";
	string mkdir_localRegionSubSeq = "mkdir -p " + localRegionSubSeqFolder;
	system(mkdir_localRegionSubSeq.c_str());
	int secondLevelIndexNormalSize = indexInfo->returnSecondLevelIndexNormalSize();
	cout << "localRegionSize: " << secondLevelIndexNormalSize << endl;
	log_ofs << "localRegionSize: " << secondLevelIndexNormalSize << endl;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		cout << "tmpChr: " << indexInfo->returnChrNameStr(tmpChr) << endl;
		log_ofs << "tmpChr: " << indexInfo->returnChrNameStr(tmpChr) << endl;
		int tmpChr_length = indexInfo->returnChromLength(tmpChr);
		cout << "tmpChr_length: " << tmpChr_length << endl;
		log_ofs << "tmpChr_length: " << tmpChr_length << endl;
		for(int tmpBase = 0; tmpBase < tmpChr_length; tmpBase += jumpLength)
		{
			int tmpSeq_startPos = tmpBase + 1;
			int tmpSeq_endPos = tmpSeq_startPos + subSeqLength - 1;
			if(tmpSeq_endPos > tmpChr_length)
				break;
			string tmpChrSubSeq = indexInfo->returnChromStrSubstr(tmpChr, tmpSeq_startPos, subSeqLength);
			if(tmpChrSubSeq.find("X") != string::npos)
				break;
			wholeGenomeSubSeqRead_ofs << ">" << indexInfo->returnChrNameStr(tmpChr) << "_"
				<< tmpSeq_startPos << "_" << tmpSeq_endPos << endl << tmpChrSubSeq << endl;
		}
		int tmpChr_subRegion_num = indexInfo->returnSecondLevelIndexPartsNum(tmpChr);
		cout << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;
		log_ofs << "tmpChr_subRegion_num: " << tmpChr_subRegion_num << endl;
		for(int tmpSubRegion = 0; tmpSubRegion < tmpChr_subRegion_num; tmpSubRegion ++)
		{
			string tmpSubRegionReadFile = localRegionSubSeqFolder + indexInfo->returnChrNameStr(tmpChr)
				+ "_" + int_to_str(tmpSubRegion + 1) + ".fa";
			ofstream tmpSubRegionRead_ofs(tmpSubRegionReadFile.c_str());
			for(int tmpBase = 0; tmpBase < secondLevelIndexNormalSize; tmpBase += jumpLength)
			{
				int tmpSeq_startPos = tmpSubRegion * secondLevelIndexNormalSize + tmpBase + 1;
				int tmpSeq_endPos = tmpSeq_startPos + subSeqLength - 1;
				if(tmpSeq_endPos > tmpChr_length)
				{
					tmpSeq_endPos = tmpChr_length - 2;
					tmpSeq_startPos = tmpSeq_endPos - subSeqLength + 1;
					if(tmpSeq_startPos < 1)
						break;
				}
				string tmpChrSubSeq = indexInfo->returnChromStrSubstr(tmpChr, tmpSeq_startPos, subSeqLength);
				if(tmpChrSubSeq.find("X") != string::npos)
					break;
				tmpSubRegionRead_ofs << ">" << indexInfo->returnChrNameStr(tmpChr) << "_"
					<< tmpSeq_startPos << "_" << tmpSeq_endPos << endl << tmpChrSubSeq << endl;
			}
			tmpSubRegionRead_ofs.close();
		}
	}
	cout << "End of job !" << endl;
	log_ofs << "End of job !" << endl;
	wholeGenomeSubSeqRead_ofs.close();
	log_ofs.close();
	return 0;
}