#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"
using namespace std;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexPath inputMPS3sj outputFilePrefix anchorSeqLength" << endl;
		exit(1);
	}
	cout << "start to initiate indexInfo" << endl;
	string anchorSeqLengthStr = argv[4];
	int anchorSeqLength = atoi(anchorSeqLengthStr.c_str());
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

	string inputMPS3sj = argv[2];
	string outputFilePrefix = argv[3];
	string output_canSJ = outputFilePrefix + "_canSJ.junc";
	ofstream canSJ_ofs(output_canSJ.c_str());
	string output_nonCanSJ_singleMutated = outputFilePrefix + "_nonCanSJ_singleMutated.junc";
	ofstream nonCanSJ_singleMutated_ofs(output_nonCanSJ_singleMutated.c_str());
	string output_nonCanSJ_others = outputFilePrefix + "_nonSJ_others.junc";
	ofstream nonCanSJ_others_ofs(output_nonCanSJ_others.c_str());
	string output_canSJ_anchorSeq = outputFilePrefix + "_canSJ.anchorSeq";
	ofstream canSJ_ofs_anchorSeq(output_canSJ_anchorSeq.c_str());
	string output_nonCanSJ_singleMutated_anchorSeq = outputFilePrefix + "_nonCanSJ_singleMutated.anchorSeq";
	ofstream nonCanSJ_singleMutated_ofs_anchorSeq(output_nonCanSJ_singleMutated_anchorSeq.c_str());
	string output_nonCanSJ_others_anchorSeq = outputFilePrefix + "_nonSJ_others.anchorSeq";
	ofstream nonCanSJ_others_ofs_anchorSeq(output_nonCanSJ_others_anchorSeq.c_str());

	ifstream mps3_SJ_ifs(inputMPS3sj.c_str());
	while(!mps3_SJ_ifs.eof())
	{
		string tmpStr;
		getline(mps3_SJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpChrName = tmpFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		string tmpStartPosStr = tmpFieldVec[1];
		int tmpStartPos = atoi(tmpStartPosStr.c_str());
		string tmpEndPosStr = tmpFieldVec[2];
		int tmpEndPos = atoi(tmpEndPosStr.c_str());
		string tmpStrand = tmpFieldVec[3];
		string tmpFlankString_raw = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpStartPos + 1, 2)
			+ indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpEndPos - 2, 2);
		string tmpAnchorSeq_raw = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpStartPos - anchorSeqLength + 1, anchorSeqLength * 2)
			+ indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpEndPos - anchorSeqLength, anchorSeqLength * 2);
		string tmpFlankString_formatted, tmpAnchorSeq_formatted;
		if(tmpStrand == "+")
		{
			tmpFlankString_formatted = tmpFlankString_raw;
			tmpAnchorSeq_formatted = tmpAnchorSeq_raw;
		}
		else if(tmpStrand == "-")
		{
			tmpFlankString_formatted = convertStringToReverseComplement(tmpFlankString_raw);
			tmpAnchorSeq_formatted = convertStringToReverseComplement(tmpAnchorSeq_raw);
		}
		else
		{
			cout << "invalid strand: " << tmpStrand << endl;
			exit(1);
		}
		if(tmpFlankString_formatted == "GTAG") // canonical
		{
			canSJ_ofs << tmpStr << endl;
			canSJ_ofs_anchorSeq << tmpAnchorSeq_formatted << endl;
		}
		else if((tmpFlankString_formatted == "GCAG")||(tmpFlankString_formatted == "ATAC")) // semi canonical
		{}
		else // noncanonical
		{
			int tmpMutatedBase = 0;
			if(tmpFlankString_formatted.at(0) != 'G')
				tmpMutatedBase ++;
			if(tmpFlankString_formatted.at(1) != 'T')
				tmpMutatedBase ++;
			if(tmpFlankString_formatted.at(2) != 'A')
				tmpMutatedBase ++;
			if(tmpFlankString_formatted.at(3) != 'G')
				tmpMutatedBase ++;
			if(tmpMutatedBase == 1)
			{
				nonCanSJ_singleMutated_ofs << tmpStr << endl;
				nonCanSJ_singleMutated_ofs_anchorSeq << tmpAnchorSeq_formatted << endl;
			}	
			else
			{
				nonCanSJ_others_ofs << tmpStr << endl;
				nonCanSJ_others_ofs_anchorSeq << tmpAnchorSeq_formatted << endl;
			}					
		}
	}
	mps3_SJ_ifs.close();
	canSJ_ofs.close();
	nonCanSJ_singleMutated_ofs.close();
	nonCanSJ_others_ofs.close();
	canSJ_ofs_anchorSeq.close();
	nonCanSJ_singleMutated_ofs_anchorSeq.close();
	nonCanSJ_others_ofs_anchorSeq.close();	
	return 0;
}