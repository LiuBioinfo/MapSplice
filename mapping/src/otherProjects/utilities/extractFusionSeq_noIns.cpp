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

#include "../../general/otherFunc.h"
#include "../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Executable inputIndexPath chrName_1 chrPos_1 strand_1 chrName_2 chrPos_2 strand_2 seqLenOnBothSides " << endl;
		exit(1);
	}
	cout << "initiate indexInfo ..." << endl;	
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

	string chrName_1 = argv[2];
	string chrName_2 = argv[5];
	int chrNameInt_1 = indexInfo->convertStringToInt(chrName_1);
	int chrNameInt_2 = indexInfo->convertStringToInt(chrName_2);
	if((chrNameInt_1 < 0)||(chrNameInt_2 < 0))
	{
		cout << "invalid chrName_1: " << chrName_1 << endl;
		cout << "invalid chrName_2: " << chrName_2 << endl;
		exit(1);
	}

	string chrPosStr_1 = argv[3];
	string chrPosStr_2 = argv[6];
	int chrPos_1 = atoi(chrPosStr_1.c_str());
	int chrPos_2 = atoi(chrPosStr_2.c_str());
	string strand_1 = argv[4];
	string strand_2 = argv[7];
	if(((strand_1 != "+")&&(strand_1 != "-"))
		||((strand_2 != "+")&&(strand_2 != "-")))
	{
		cout << "invalid strand_1: " << strand_1 << endl;
		cout << "invalid strand_2: " << strand_2 << endl;
		exit(1);
	}

	string seqLenOnBothSidesStr = argv[8];
	int seqLenOnBothSides = atoi(seqLenOnBothSidesStr.c_str());

	string seq_1, seq_2, seq_1_rcm, seq_2_rcm;
	int startPos_1, startPos_2;
	string flankString_1, flankString_2;
	if(strand_1 == "+")
	{	
		startPos_1 = chrPos_1 - seqLenOnBothSides + 1;
		seq_1 = indexInfo->returnChromStrSubstr(chrNameInt_1, startPos_1, seqLenOnBothSides);
		flankString_1 = indexInfo->returnChromStrSubstr(chrNameInt_1, startPos_1 + 1, 2);
	}
	else
	{
		startPos_1 = chrPos_1;
		seq_1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_1, startPos_1, seqLenOnBothSides));		
		flankString_1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_1, startPos_1 - 2, 2));
	}
	if(strand_2 == "+")
	{
		startPos_2 = chrPos_2;
		seq_2 = indexInfo->returnChromStrSubstr(chrNameInt_2, startPos_2, seqLenOnBothSides);
		flankString_2 = indexInfo->returnChromStrSubstr(chrNameInt_2, startPos_2 - 2, 2);
	}
	else
	{
		startPos_2 = chrPos_2 - seqLenOnBothSides + 1;
		seq_2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_2, startPos_2, seqLenOnBothSides));
		flankString_2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(chrNameInt_2, startPos_2 + 1, 2));
	}
	seq_1_rcm = convertStringToReverseComplement(seq_1);
	seq_2_rcm = convertStringToReverseComplement(seq_2);
	cout << "GENE_1 " << chrName_1 << " : " <<  startPos_1 << " ~ " << startPos_1 + seqLenOnBothSides - 1 << " " + strand_1 << endl;
	cout << seq_1 << endl << seq_1_rcm << endl;
	cout << "GENE_2 " << chrName_2 << " : " <<  startPos_2 << " ~ " << startPos_2 + seqLenOnBothSides - 1 << " " + strand_2 << endl;
	cout << seq_2 << endl << seq_2_rcm << endl;
	cout << "fankString: " << flankString_1 << flankString_2 << endl;
	return 0;
}