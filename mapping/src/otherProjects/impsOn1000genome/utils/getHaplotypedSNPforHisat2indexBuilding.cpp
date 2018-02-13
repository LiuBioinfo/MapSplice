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
#include <sstream>
#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../../incorporateGenomicVariants/general/SNPhash_info.h"

using namespace std;

void get_startPos_endPos_idVec_hap1(
	int& tmpChr_hap1_startPos, int& tmpChr_hap1_endPos, vector<string>& tmpChr_hap1_snpIdVec,
	vector<int>& tmpSnpPosVec, vector<string>& tmpSnpIdVec, vector<int>& tmpSnpHapInfoVec)
{
	bool startPos_assigned_bool = false;
	for(int tmp = 0; tmp < tmpSnpPosVec.size(); tmp++)
	{
		int tmpSNPpos = tmpSnpPosVec[tmp];
		string tmpSNPid = tmpSnpIdVec[tmp];
		int tmpSnpHapInfo = tmpSnpHapInfoVec[tmp];
		if((tmpSnpHapInfo == 1)||(tmpSnpHapInfo == 0))
		{
			if(!startPos_assigned_bool)
			{	
				tmpChr_hap1_startPos = tmpSNPpos;
				startPos_assigned_bool = true;
			}
			tmpChr_hap1_endPos = tmpSNPpos;
			tmpChr_hap1_snpIdVec.push_back(tmpSNPid);
		}
	}
}

void get_startPos_endPos_idVec_hap2(
	int& tmpChr_hap2_startPos, int& tmpChr_hap2_endPos, vector<string>& tmpChr_hap2_snpIdVec,
	vector<int>& tmpSnpPosVec, vector<string>& tmpSnpIdVec, vector<int>& tmpSnpHapInfoVec)
{
	bool startPos_assigned_bool = false;
	for(int tmp = 0; tmp < tmpSnpPosVec.size(); tmp++)
	{
		int tmpSNPpos = tmpSnpPosVec[tmp];
		string tmpSNPid = tmpSnpIdVec[tmp];
		int tmpSnpHapInfo = tmpSnpHapInfoVec[tmp];
		if((tmpSnpHapInfo == 2)||(tmpSnpHapInfo == 0))
		{
			if(!startPos_assigned_bool)
			{	
				tmpChr_hap2_startPos = tmpSNPpos;
				startPos_assigned_bool = true;
			}
			tmpChr_hap2_endPos = tmpSNPpos;
			tmpChr_hap2_snpIdVec.push_back(tmpSNPid);
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputIndexFolderPath" << endl;
		cout << "#2 inputSNPfile_unique_1" << endl;
		cout << "#3 inputSNPfile_unique_2" << endl;
		cout << "#4 inputSNPfile_shared" << endl; 
		cout << "#5 outputDir" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputSNPfile_unique_1 = argv[2];
	string inputSNPfile_unique_2 = argv[3];
	string inputSNPfile_shared = argv[4];
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputFile_log = outputFolderStr + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	log_ofs << "Command: \n" << argv[0] << endl << argv[1] << endl 
		<< argv[2] << endl << argv[3] << endl << argv[4] << endl;

	cout << "initiate indexInfo ..." << endl;
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
	parameter_ifs.close();

	string total_SNP_file = outputFolderStr + "total_SNP.txt";
	ofstream SNP_total_ofs(total_SNP_file.c_str());
	ifstream SNP_ifs_hap1(inputSNPfile_unique_1.c_str());
	while(!SNP_ifs_hap1.eof())
	{
		string tmpStr;
		getline(SNP_ifs_hap1, tmpStr);
		if(tmpStr == "")
			break;
		SNP_total_ofs << tmpStr << "\tHap1" << endl;
	}
	SNP_ifs_hap1.close();
	ifstream SNP_ifs_hap2(inputSNPfile_unique_2.c_str());
	while(!SNP_ifs_hap2.eof())
	{
		string tmpStr;
		getline(SNP_ifs_hap2, tmpStr);
		if(tmpStr == "")
			break;
		SNP_total_ofs << tmpStr << "\tHap2" << endl;
	}
	SNP_ifs_hap2.close();
	ifstream SNP_ifs_twoHap(inputSNPfile_shared.c_str());
	while(!SNP_ifs_twoHap.eof())
	{
		string tmpStr;
		getline(SNP_ifs_twoHap, tmpStr);
		if(tmpStr == "")
			break;
		SNP_total_ofs << tmpStr << "\tTwoHap" << endl;
	}
	SNP_ifs_twoHap.close();
	SNP_total_ofs.close();

	string invalid_SNP_file = outputFolderStr + "invalid_SNP.txt";
	ofstream invalidSNP_ofs(invalid_SNP_file.c_str());
	string total_SNP_file_sorted = outputFolderStr + "total_SNP.sorted.txt";
	string cmd_sort_snp = "sort -k2 -n " + total_SNP_file 
		+ " > " + total_SNP_file_sorted;
	system(cmd_sort_snp.c_str());

	vector< vector<int> > snpPosVecVec;
	vector< vector<string> > snpBaseVecVec;
	vector< vector<int> > snpHapInfoVecVec; // 1 -- hap1; 2 -- hap2; 0 -- twoHap;
 	
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		vector<int> tmpSNPposVec;
		vector<string> tmpSNPbaseVec;
		vector<int> tmpSNPhapInfoVec;
		snpPosVecVec.push_back(tmpSNPposVec);
		snpBaseVecVec.push_back(tmpSNPbaseVec);
		snpHapInfoVecVec.push_back(tmpSNPhapInfoVec);
	}

	ifstream total_SNP_ifs(total_SNP_file_sorted.c_str());
	while(!total_SNP_ifs.eof())
	{
		string tmpStr;
		getline(total_SNP_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		string tmpChrName = tmpStr.substr(0, tabLoc_1);
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		string tmpChrPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		int tmpChrPos = atoi(tmpChrPosStr.c_str());
		string tmpRefBase = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpChrPos, 1);
		string tmpAltBase = tmpStr.substr(tabLoc_3 + 1, 1);
		if(tmpRefBase == tmpAltBase)
		{
			invalidSNP_ofs << tmpStr << endl;
			continue;
		}
		string tmpHapInfo = tmpStr.substr(tabLoc_4 + 1);
	
		snpPosVecVec[tmpChrNameInt].push_back(tmpChrPos);
		snpBaseVecVec[tmpChrNameInt].push_back(tmpAltBase);		
		if(tmpHapInfo == "TwoHap")
			snpHapInfoVecVec[tmpChrNameInt].push_back(0);
		else if(tmpHapInfo == "Hap1")
			snpHapInfoVecVec[tmpChrNameInt].push_back(1);
		else if(tmpHapInfo == "Hap2")
			snpHapInfoVecVec[tmpChrNameInt].push_back(2);
		else
		{
			cout << "invalid tmpHapInfo: " << tmpHapInfo << endl;
			exit(1); 
		}
	}
	total_SNP_ifs.close();
	invalidSNP_ofs.close();

	vector< vector<string> > snpIdVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		vector<string> tmpSNPidVec;
		snpIdVecVec.push_back(tmpSNPidVec);
	}
	int tmpSNP_NO = 0;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{	
		int tmpChr_snpNum = snpPosVecVec[tmpChr].size();
		for(int tmpSNP = 0; tmpSNP < tmpChr_snpNum; tmpSNP ++)
		{
			tmpSNP_NO ++;
			string tmpSNP_id = "SNP_" + int_to_str(tmpSNP_NO);
			snpIdVecVec[tmpChr].push_back(tmpSNP_id);
		}
	}

	string SNP_info_hisat2_file = outputFolderStr + "SNP_info.hisat2.txt";
	string hap_info_hisat2_file = outputFolderStr + "hap_info.hisat2.txt";
	ofstream snp_ofs(SNP_info_hisat2_file.c_str());
	ofstream hap_ofs(hap_info_hisat2_file.c_str());
	
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		int tmpChr_snpNum = snpPosVecVec[tmpChr].size();
		// snp
		for(int tmpSNP = 0; tmpSNP < tmpChr_snpNum; tmpSNP ++)
		{
			string tmpSNP_id = (snpIdVecVec[tmpChr])[tmpSNP];
			snp_ofs << tmpSNP_id << "\tsingle\t" << tmpChrName << "\t" 
				<< (snpPosVecVec[tmpChr])[tmpSNP] - 1 << "\t" << (snpBaseVecVec[tmpChr])[tmpSNP] << endl;
		}
		// hap1
		int tmpChr_hap1_startPos, tmpChr_hap1_endPos;
		vector<string> tmpChr_hap1_snpIdVec;
		get_startPos_endPos_idVec_hap1(tmpChr_hap1_startPos, tmpChr_hap1_endPos, tmpChr_hap1_snpIdVec,
			snpPosVecVec[tmpChr], snpIdVecVec[tmpChr], snpHapInfoVecVec[tmpChr]);
		if(tmpChr_hap1_snpIdVec.size() > 0)
		{	
			hap_ofs << "hap" << tmpChr * 2 + 1 << "\t" << tmpChrName << "\t" << tmpChr_hap1_startPos - 1 << "\t" << tmpChr_hap1_endPos - 1 << "\t";
			hap_ofs << tmpChr_hap1_snpIdVec[0];
			if(tmpChr_hap1_snpIdVec.size() > 1)
			{
				for(int tmpSNPindex = 1; tmpSNPindex < tmpChr_hap1_snpIdVec.size(); tmpSNPindex++)
					hap_ofs << "," << tmpChr_hap1_snpIdVec[tmpSNPindex];
			}
			hap_ofs << endl;
		}
		// hap2
		int tmpChr_hap2_startPos, tmpChr_hap2_endPos;
		vector<string> tmpChr_hap2_snpIdVec;
		get_startPos_endPos_idVec_hap2(tmpChr_hap2_startPos, tmpChr_hap2_endPos, tmpChr_hap2_snpIdVec,
			snpPosVecVec[tmpChr], snpIdVecVec[tmpChr], snpHapInfoVecVec[tmpChr]);
		if(tmpChr_hap2_snpIdVec.size() > 0)
		{	
			hap_ofs << "hap" << tmpChr * 2 + 2 << "\t" << tmpChrName << "\t" << tmpChr_hap2_startPos - 1 << "\t" << tmpChr_hap2_endPos - 1 << "\t";
			hap_ofs << tmpChr_hap2_snpIdVec[0];
			if(tmpChr_hap2_snpIdVec.size() > 1)
			{
				for(int tmpSNPindex = 1; tmpSNPindex < tmpChr_hap2_snpIdVec.size(); tmpSNPindex++)
					hap_ofs << "," << tmpChr_hap2_snpIdVec[tmpSNPindex];
			}
			hap_ofs << endl;
		}
	}
	snp_ofs.close();
	hap_ofs.close();

	delete indexInfo;
	free(chrom);
	log_ofs.close();
	return 0;
}