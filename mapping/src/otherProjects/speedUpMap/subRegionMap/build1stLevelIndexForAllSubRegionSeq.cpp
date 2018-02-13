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
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputSubRegionSeqFolder outputSubRegionIndexFolder outputScript" << endl;
		exit(1);
	}
	cout << "creating folder ......" << endl;
	string inputSubRegionSeqFolder = argv[2];
	string outputFolderStr = argv[3];
	string outputScript = argv[4];
	ofstream script_ofs(outputScript.c_str());
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
			string tmpChrSubRegionIndexFolder = outputFolderStr + "/" + tmpChrName + "/" + int_to_str(tmpSubRegion + 1);
			string mkdir_tmpChrSubRegionIndexFolder = "mkdir -p " + tmpChrSubRegionIndexFolder;
			system(mkdir_tmpChrSubRegionIndexFolder.c_str());
			string tmpChrSubRegionSeqFolder = inputSubRegionSeqFolder + "/" + tmpChrName + "/" + int_to_str(tmpSubRegion + 1) + "/";
			script_ofs << "./buildWholeGenome " << tmpChrSubRegionSeqFolder << " \\\n" << tmpChrSubRegionIndexFolder << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.1" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.2" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.3" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.4" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.5" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.6" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.7" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.8" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.9" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.10" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.11" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.12" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.13" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_preIndexString.14" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_lcp" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_up" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_down" << endl;
			script_ofs << "rm " << tmpChrSubRegionIndexFolder << "/_next" << endl;
			script_ofs << endl;
		}
	}
	cout << "End of job !" << endl;
	log_ofs << "End of job !" << endl;
	//wholeGenomeSubSeqRead_ofs.close();
	script_ofs.close();
	log_ofs.close();
	return 0;
}