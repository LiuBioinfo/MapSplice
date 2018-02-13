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
#include "../general/geneAnnEntryHash.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexPath inputGeneAnnEntryFile outputIndexFile" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	//string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	//ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;
	indexInfo->initiate_withoutLoadingSeq();
	cout << "end of initiating indexInfo" << endl;

	string inputGeneAnnEntryFile = argv[2];
	cout << "start to load simplified gene annotation file" << endl;
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);
	cout << "start to output index" << endl;
	string outputIndexFile = argv[3];
	tmpGeneAnnHashInfo.outputIndex(outputIndexFile, indexInfo);

	delete indexInfo;
	parameter_ifs.close();
	//chrom_bit_file_ifs.close();
	//free(chrom);
	return 0;
}