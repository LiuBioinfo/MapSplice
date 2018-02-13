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
#include "../../../general/index_info.h"
using namespace std;

// string covertCharToReverseComplement(const string& Ori_Char)
// {
// 	if(Ori_Char == "A")
// 	{
// 		return "T";
// 	}
// 	else if(Ori_Char == "T")
// 	{
// 		return "A";
// 	}
// 	else if(Ori_Char == "G")
// 	{
// 		return "C";
// 	}
// 	else if(Ori_Char == "C")
// 	{
// 		return "G";
// 	}
// 	else if(Ori_Char == "N")
// 	{
// 		return "N";
// 	}
// 	else
// 	{
// 		cout << "incorrect Ori_Char in covertCharToReverseComplement" << endl;
// 		exit(1);
// 		return "X";
// 	}
// }

// string convertStringToReverseComplement(const string& originalString)
// {
// 	int stringLength = originalString.size();
// 	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
// 	for (int tmp = 1; tmp < stringLength; tmp++)
// 	{
// 		resultString = resultString + covertCharToReverseComplement(
// 			originalString.substr(stringLength-1-tmp, 1));
// 	}
// 	return resultString;
// }

int main(int argc, char** argv)
{
	int anchorSeqLength = 30;
	//int halfSeqLength = seqLength/2;
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexDir" << endl;
		cout << "#2 inputBedpeFile" << endl;
		cout << "#3 outputBreakpointAnchorSeqFile" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputBedpeFile = argv[2];
	string outputBreakpointAnchorSeqFile = argv[3];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);	
	free(chrom);

	vector<string> chrNameVec_gene1;
	vector<string> chrNameVec_gene2;
	vector<int> chrPosVec_gene1;
	vector<int> chrPosVec_gene2;
	vector<bool> forOrRevVec_gene1;
	vector<bool> forOrRevVec_gene2;
	vector<string> rawFusionStrVec;
	ifstream bedpe_ifs(inputBedpeFile.c_str());
	while(!bedpe_ifs.eof())
	{
		string tmpStr;
		getline(bedpe_ifs, tmpStr);
		if(tmpStr == "")
			break;
		rawFusionStrVec.push_back(tmpStr);
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
		int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
		int tabLoc_8 = tmpStr.find("\t", tabLoc_7 + 1);
		int tabLoc_9 = tmpStr.find("\t", tabLoc_8 + 1);
		string tmpStr_gene1_chrName = "chr" + tmpStr.substr(0, tabLoc_1);
		string tmpStr_gene1_pos1 = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpStr_gene1_pos2 = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpStr_gene2_chrName = "chr" + tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string tmpStr_gene2_pos1 = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string tmpStr_gene2_pos2 = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		string tmpStr_geneIdPair = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
		string tmpStr_score = tmpStr.substr(tabLoc_7 + 1, tabLoc_8 - tabLoc_7 - 1);
		string tmpStr_gene1_strandId = tmpStr.substr(tabLoc_8 + 1, tabLoc_9 - tabLoc_8 - 1);
		string tmpStr_gene2_strandId = tmpStr.substr(tabLoc_9 + 1);
		bool tmpStrand_gene1, tmpStrand_gene2;
		if(tmpStr_gene1_strandId == "1")
			tmpStrand_gene1 = true;
		else if(tmpStr_gene1_strandId == "-1")
			tmpStrand_gene1 = false;
		else
		{
			cout << "error! invalid strand" << endl;
			exit(1);
		}

		if(tmpStr_gene2_strandId == "1")
			tmpStrand_gene2 = true;
		else if(tmpStr_gene2_strandId == "-1")
			tmpStrand_gene2 = false;
		else
		{
			cout << "error! invalid strand" << endl;
			exit(1);
		}
		int tmpPos_gene1, tmpPos_gene2;
		if(tmpStrand_gene1)
			tmpPos_gene1 = atoi(tmpStr_gene1_pos2.c_str());
		else
			tmpPos_gene1 = atoi(tmpStr_gene1_pos1.c_str());
		if(tmpStrand_gene2)
			tmpPos_gene2 = atoi(tmpStr_gene2_pos1.c_str());
		else
			tmpPos_gene2 = atoi(tmpStr_gene2_pos2.c_str());
		
		chrNameVec_gene1.push_back(tmpStr_gene1_chrName);
		chrNameVec_gene2.push_back(tmpStr_gene2_chrName);
		chrPosVec_gene1.push_back(tmpPos_gene1);
		chrPosVec_gene2.push_back(tmpPos_gene2);
		forOrRevVec_gene1.push_back(tmpStrand_gene1);
		forOrRevVec_gene2.push_back(tmpStrand_gene2);
	}
	bedpe_ifs.close();

	ofstream seq_ofs(outputBreakpointAnchorSeqFile.c_str());
	for(int tmp = 0; tmp < chrNameVec_gene1.size(); tmp++)
	{
		seq_ofs << rawFusionStrVec[tmp] << endl;
		string tmpChrName_gene1 = chrNameVec_gene1[tmp];
		int tmpChrNameInt_gene1 = indexInfo->convertStringToInt(tmpChrName_gene1);
		int tmpChrPos_gene1 = chrPosVec_gene1[tmp];
		bool tmpStrandBool_gene1 = forOrRevVec_gene1[tmp];
		if(tmpStrandBool_gene1)
		{	
			string tmpStrandStr_gene1 = "+";
			int tmpStartPos_gene1 = tmpChrPos_gene1 - anchorSeqLength + 1;
			int tmpEndPos_gene1 = tmpChrPos_gene1;
			seq_ofs << tmpChrName_gene1 << "\t" << tmpStartPos_gene1 << "\t" << tmpEndPos_gene1 << "\t" << tmpStrandStr_gene1 << endl;	
			string tmpSeq_for_gene1 = indexInfo->returnChromStrSubstr(tmpChrNameInt_gene1, tmpStartPos_gene1, anchorSeqLength);
			string tmpSeq_rev_gene1 = convertStringToReverseComplement(tmpSeq_for_gene1);	
			seq_ofs << tmpSeq_for_gene1 << endl << tmpSeq_rev_gene1 << endl;	
		}
		else
		{
			string tmpStrandStr_gene1 = "-";
			int tmpStartPos_gene1 = tmpChrPos_gene1;
			int tmpEndPos_gene1 = tmpChrPos_gene1 + anchorSeqLength - 1;
			seq_ofs << tmpChrName_gene1 << "\t" << tmpStartPos_gene1 << "\t" << tmpEndPos_gene1 << "\t" << tmpStrandStr_gene1 << endl;
			string tmpSeq_for_gene1 = indexInfo->returnChromStrSubstr(tmpChrNameInt_gene1, tmpStartPos_gene1, anchorSeqLength);
			string tmpSeq_rev_gene1 = convertStringToReverseComplement(tmpSeq_for_gene1);	
			seq_ofs << tmpSeq_for_gene1 << endl << tmpSeq_rev_gene1 << endl;
		}

		
		string tmpChrName_gene2 = chrNameVec_gene2[tmp];
		int tmpChrNameInt_gene2 = indexInfo->convertStringToInt(tmpChrName_gene2);
		int tmpChrPos_gene2 = chrPosVec_gene2[tmp];
		bool tmpStrandBool_gene2 = forOrRevVec_gene2[tmp];
		if(tmpStrandBool_gene2)
		{	
			string tmpStrandStr_gene2 = "+";
			int tmpStartPos_gene2 = tmpChrPos_gene2;// 
			int tmpEndPos_gene2 = tmpChrPos_gene2 + anchorSeqLength - 1;
			seq_ofs << tmpChrName_gene2 << "\t" << tmpStartPos_gene2 << "\t" << tmpEndPos_gene2 << "\t" << tmpStrandStr_gene2 << endl;
			string tmpSeq_for_gene2 = indexInfo->returnChromStrSubstr(tmpChrNameInt_gene2, tmpStartPos_gene2, anchorSeqLength);
			string tmpSeq_rev_gene2 = convertStringToReverseComplement(tmpSeq_for_gene2);
			seq_ofs << tmpSeq_for_gene2 << endl << tmpSeq_rev_gene2 << endl;
		}
		else
		{
			string tmpStrandStr_gene2 = "-";
			int tmpStartPos_gene2 = tmpChrPos_gene2 - anchorSeqLength + 1;;
			int tmpEndPos_gene2 = tmpChrPos_gene2;
			seq_ofs << tmpChrName_gene2 << "\t" << tmpStartPos_gene2 << "\t" << tmpEndPos_gene2 << "\t" << tmpStrandStr_gene2 << endl;
			string tmpSeq_for_gene2 = indexInfo->returnChromStrSubstr(tmpChrNameInt_gene2, tmpStartPos_gene2, anchorSeqLength);
			string tmpSeq_rev_gene2 = convertStringToReverseComplement(tmpSeq_for_gene2);	
			seq_ofs << tmpSeq_for_gene2 << endl << tmpSeq_rev_gene2 << endl;
		}
		seq_ofs << endl;
	}
	seq_ofs.close();
	delete indexInfo;
	return 0;
}