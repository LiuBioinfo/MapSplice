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
#include <sstream>
#include "../../../../general/index_info.h"
#include "../general/gene_transcript_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexPath inputEnsemblGTF outputFolder outputGeneRegionFaOrNot outputTranscriptSeqFaOrNot" << endl;
		exit(1);
	}
	string outputGeneRegionFaOrNot = argv[4];
	string outputTranscriptSeqFaOrNot = argv[5];	
	bool printOutGeneRegionFa_bool;
	if(outputGeneRegionFaOrNot == "Y")
		printOutGeneRegionFa_bool = true;
	else if(outputGeneRegionFaOrNot == "N")
		printOutGeneRegionFa_bool = false;
	else
	{
		cout << "invalid option for printOutGeneRegionFa_bool: " << printOutGeneRegionFa_bool << endl;
		exit(1);
	}
	bool printOutTranscriptSeqFa_bool;
	if(outputTranscriptSeqFaOrNot == "Y")
		printOutTranscriptSeqFa_bool = true;
	else if(outputTranscriptSeqFaOrNot == "N")
		printOutTranscriptSeqFa_bool = false;
	else
	{
		cout << "invalid option for printOutTranscriptSeqFa_bool: " << printOutTranscriptSeqFa_bool << endl;
		exit(1);
	}
	if((!printOutGeneRegionFa_bool)&&(!printOutTranscriptSeqFa_bool))
	{
		cout << "error! (!printOutGeneRegionFa_bool)&&(!printOutTranscriptSeqFa_bool)" << endl;
		exit(1);
	}

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

	string inputEnsemblGTF = argv[2];
	string outputFolder = argv[3];
	cout << "inputEnsemblGTF: " << inputEnsemblGTF << endl;
	Gene_Transcript_Vec_Info geneTranscriptVecInfo;
	geneTranscriptVecInfo.initaite_ensemblGeneAnn(inputEnsemblGTF, indexInfo);
	if(printOutGeneRegionFa_bool && printOutTranscriptSeqFa_bool)
		geneTranscriptVecInfo.mkDir_printoutFa_gene_transcript(outputFolder, indexInfo);
	else if(printOutGeneRegionFa_bool)
		geneTranscriptVecInfo.mkDir_printoutFa_gene(outputFolder, indexInfo);
	else if(printOutTranscriptSeqFa_bool)
		geneTranscriptVecInfo.mkDir_printoutFa_transcript(outputFolder, indexInfo);
	else
	{
		cout << "error ! (!printOutGeneRegionFa_bool)&&(!printOutTranscriptSeqFa_bool)" << endl;
		exit(1);
	}

	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	return 0;
}