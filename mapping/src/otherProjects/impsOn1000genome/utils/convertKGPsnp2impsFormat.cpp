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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

void getIdVec(string& dataListFile, vector<string>& idVec)
{
	ifstream id_ifs(dataListFile.c_str());
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		if(tabLoc == string::npos)
			idVec.push_back(tmpStr);
		else
			idVec.push_back(tmpStr.substr(0, tabLoc));
	}
	id_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexDir" << endl;
		cout << "#2 inputIdList" << endl;
		cout << "#3 inputDataDir" << endl;
		cout << "#4 outputDir" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputIdList = argv[2];
	string inputDataDir = argv[3];
	string outputDir = argv[4];
	inputDataDir += "/";
	outputDir += "/";

	indexFolderPath += "/";
	cout << "start to initiate indexInfo" << endl;
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
	free(chrom);

	vector<string> idVec;
	getIdVec(inputIdList, idVec);

	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		string input_snp_hap1 = inputDataDir + tmpId + "/" + tmpId + "_1.snp";
		string input_snp_hap2 = inputDataDir + tmpId + "/" + tmpId + "_2.snp";
		string output_snp_hap1 = outputDir + tmpId + "_hap1_SNP.txt";
		string output_snp_hap2 = outputDir + tmpId + "_hap2_SNP.txt";
		// hap1 snp
		ifstream snp_hap1_ifs(input_snp_hap1.c_str());
		ofstream snp_hap1_ofs(output_snp_hap1.c_str());
		while(!snp_hap1_ifs.eof())
		{
			string tmpStr;
			getline(snp_hap1_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
			string tmpChr = tmpStr.substr(0, tabLoc_1);
			int tmpChrInt = indexInfo->convertStringToInt(tmpChr);
			if(tmpChrInt < 0)
				continue;
			string tmpPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			string tmpSNPinfo = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
			if(tmpSNPinfo.length() != 4)
				continue;
			string tmpRefBase = tmpSNPinfo.substr(0, 1);
			string tmpAltBase = tmpSNPinfo.substr(3, 1);
			if(((tmpRefBase != "A")&&(tmpRefBase != "C")&&(tmpRefBase != "G")&&(tmpRefBase != "T")&&(tmpRefBase != "N")
				&&(tmpRefBase != "a")&&(tmpRefBase != "c")&&(tmpRefBase != "g")&&(tmpRefBase != "t")&&(tmpRefBase != "n"))
				||((tmpAltBase != "A")&&(tmpAltBase != "C")&&(tmpAltBase != "G")&&(tmpAltBase != "T")&&(tmpAltBase != "N")
					&&(tmpAltBase != "a")&&(tmpAltBase != "c")&&(tmpAltBase != "g")&&(tmpAltBase != "t")&&(tmpAltBase != "n")))
				continue;
			
			if(tmpRefBase == "a")
				tmpRefBase = "A";
			else if(tmpRefBase == "c")
				tmpRefBase = "C";
			else if(tmpRefBase == "g")
				tmpRefBase = "G";
			else if(tmpRefBase == "t")
				tmpRefBase = "T";
			else
			{}

			if(tmpAltBase == "a")
				tmpAltBase = "A";
			else if(tmpAltBase == "c")
				tmpAltBase = "C";
			else if(tmpAltBase == "g")
				tmpAltBase = "G";
			else if(tmpAltBase == "t")
				tmpAltBase = "T";
			else
			{}
			snp_hap1_ofs << tmpChr << "\t" << tmpPosStr << "\t" << tmpRefBase << "\t" << tmpAltBase << endl;
		}
		snp_hap1_ofs.close();
		snp_hap1_ifs.close();
		// hap2 snp
		ifstream snp_hap2_ifs(input_snp_hap2.c_str());
		ofstream snp_hap2_ofs(output_snp_hap2.c_str());
		while(!snp_hap2_ifs.eof())
		{
			string tmpStr;
			getline(snp_hap2_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
			string tmpChr = tmpStr.substr(0, tabLoc_1);
			int tmpChrInt = indexInfo->convertStringToInt(tmpChr);
			if(tmpChrInt < 0)
				continue;
			string tmpPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			string tmpSNPinfo = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
			if(tmpSNPinfo.length() != 4)
				continue;
			string tmpRefBase = tmpSNPinfo.substr(0, 1);
			string tmpAltBase = tmpSNPinfo.substr(3, 1);
			if(((tmpRefBase != "A")&&(tmpRefBase != "C")&&(tmpRefBase != "G")&&(tmpRefBase != "T")&&(tmpRefBase != "N")
				&&(tmpRefBase != "a")&&(tmpRefBase != "c")&&(tmpRefBase != "g")&&(tmpRefBase != "t")&&(tmpRefBase != "n"))
				||((tmpAltBase != "A")&&(tmpAltBase != "C")&&(tmpAltBase != "G")&&(tmpAltBase != "T")&&(tmpAltBase != "N")
					&&(tmpAltBase != "a")&&(tmpAltBase != "c")&&(tmpAltBase != "g")&&(tmpAltBase != "t")&&(tmpAltBase != "n")))
				continue;
			
			if(tmpRefBase == "a")
				tmpRefBase = "A";
			else if(tmpRefBase == "c")
				tmpRefBase = "C";
			else if(tmpRefBase == "g")
				tmpRefBase = "G";
			else if(tmpRefBase == "t")
				tmpRefBase = "T";
			else
			{}

			if(tmpAltBase == "a")
				tmpAltBase = "A";
			else if(tmpAltBase == "c")
				tmpAltBase = "C";
			else if(tmpAltBase == "g")
				tmpAltBase = "G";
			else if(tmpAltBase == "t")
				tmpAltBase = "T";
			else
			{}
			snp_hap2_ofs << tmpChr << "\t" << tmpPosStr << "\t" << tmpRefBase << "\t" << tmpAltBase << endl;
		}
		snp_hap2_ofs.close();
		snp_hap2_ifs.close();
	}

	return 0;
}