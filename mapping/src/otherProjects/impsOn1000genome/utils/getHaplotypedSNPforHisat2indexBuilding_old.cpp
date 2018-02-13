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

	SNPhash_Info tmpSNPhashInfo_unique_hap1;
	SNPhash_Info tmpSNPhashInfo_unique_hap2;
	SNPhash_Info tmpSNPhashInfo_shared;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo_unique_hap1.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_unique_hap2.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_shared.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo_unique_hap1.generateSNPhash_formattedSNPfile(inputSNPfile_unique_1, indexInfo);
	tmpSNPhashInfo_unique_hap2.generateSNPhash_formattedSNPfile(inputSNPfile_unique_2, indexInfo);
	tmpSNPhashInfo_shared.generateSNPhash_formattedSNPfile(inputSNPfile_shared, indexInfo);

	int SNP_num_shared = tmpSNPhashInfo_shared.returnSNPnum();
	int SNP_num_unique_hap1 = tmpSNPhashInfo_unique_hap1.returnSNPnum();
	int SNP_num_unique_hap2 = tmpSNPhashInfo_unique_hap2.returnSNPnum();
	vector< vector<int> > snpIndexVecVec_shared;
	vector< vector<int> > snpIndexVecVec_unique_hap1;
	vector< vector<int> > snpIndexVecVec_unique_hap2;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		vector<int> tmp_snpIndexVec_shared;
		vector<int> tmp_snpIndexVec_unique_hap1;
		vector<int> tmp_snpIndexVec_unique_hap2;
		snpIndexVecVec_shared.push_back(tmp_snpIndexVec_shared);
		snpIndexVecVec_unique_hap1.push_back(tmp_snpIndexVec_unique_hap1);
		snpIndexVecVec_unique_hap2.push_back(tmp_snpIndexVec_unique_hap2);
	}

	string SNP_info_hisat2_file = outputFolderStr + "SNP_info.hisat2.txt";
	ofstream snp_ofs(SNP_info_hisat2_file.c_str());
	
	// shared SNP
	int tmpSNP_id_1st_shared = 1;
	for(int tmp = 0; tmp < SNP_num_shared; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo_shared.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo_shared.returnSNP_chrPos(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo_shared.returnSNP_alterBase(tmp);
		snp_ofs << "SNP_" << tmpSNP_id_1st_shared + tmp << "\tsingle\t" << tmpSNP_chrNameStr 
			<< "\t" << tmpSNP_chrPos - 1 << "\t" << tmpSNP_alterBase << endl;
		snpIndexVecVec_shared[tmpSNP_chrNameInt].push_back(tmp);
	}

	// hap1 SNP
	int tmpSNP_id_1st_unique_hap1 = 1 + SNP_num_shared;
	for(int tmp = 0; tmp < SNP_num_unique_hap1; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo_unique_hap1.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo_unique_hap1.returnSNP_chrPos(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo_unique_hap1.returnSNP_alterBase(tmp);
		snp_ofs << "SNP_" << tmpSNP_id_1st_unique_hap1 + tmp << "\tsingle\t" << tmpSNP_chrNameStr 
			<< "\t" << tmpSNP_chrPos - 1 << "\t" << tmpSNP_alterBase << endl;
		snpIndexVecVec_unique_hap1[tmpSNP_chrNameInt].push_back(tmp);		
	}
	// hap2 SNP
	int tmpSNP_id_1st_unique_hap2 = 1 + SNP_num_shared + SNP_num_unique_hap1;
	for(int tmp = 0; tmp < SNP_num_unique_hap2; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo_unique_hap2.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo_unique_hap2.returnSNP_chrPos(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo_unique_hap2.returnSNP_alterBase(tmp);
		snp_ofs << "SNP_" << tmpSNP_id_1st_unique_hap2 + tmp << "\tsingle\t" << tmpSNP_chrNameStr 
			<< "\t" << tmpSNP_chrPos - 1 << "\t" << tmpSNP_alterBase << endl;
		snpIndexVecVec_unique_hap2[tmpSNP_chrNameInt].push_back(tmp);		
	}
	snp_ofs.close();
	string hap_info_hisat2_file = outputFolderStr + "hap_info.hisat2.txt";
	ofstream hap_ofs(hap_info_hisat2_file.c_str());
	cout << "start to print hapInfo" << endl;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		cout << "tmpChromIndex: " << tmp << endl;
		string tmpChrName = indexInfo->returnChrNameStr(tmp);
		cout << "tmpChrName: " << tmpChrName << endl;
		int tmpChrLength = indexInfo->returnChromLength(tmp);
		if(tmp == chromNum - 1)
			tmpChrLength = tmpChrLength - 1;
		cout << "tmpChrLength: " << tmpChrLength << endl;
		if((snpIndexVecVec_shared[tmp]).size() > 0)
		{	
			// hap1 
			hap_ofs << "hap1\t" << tmpChrName << "\t1\t" << tmpChrLength-1 << "\t";			
			// print the 1st shared SNP from tmpChr 
			hap_ofs << "SNP_" << ((snpIndexVecVec_shared[tmp])[0] + tmpSNP_id_1st_shared);
			// print other shared SNPs from tmpChr
			for(int tmpSNPindex = 1; tmpSNPindex < (snpIndexVecVec_shared[tmp]).size(); tmpSNPindex ++)
				hap_ofs << ",SNP_" << ((snpIndexVecVec_shared[tmp])[tmpSNPindex] + tmpSNP_id_1st_shared);
			// print hap1 unique SNPs from tmpChr
			for(int tmpSNPindex = 0; tmpSNPindex < (snpIndexVecVec_unique_hap1[tmp]).size(); tmpSNPindex ++)
				hap_ofs << ",SNP_" << ((snpIndexVecVec_unique_hap1[tmp])[tmpSNPindex] + tmpSNP_id_1st_unique_hap1);
			hap_ofs << endl;
			// hap2
			hap_ofs << "hap2\t" << tmpChrName << "\t1\t" << tmpChrLength-1 << "\t";
			// print the 1st shared SNP from tmpChr 
			hap_ofs << "SNP_" << ((snpIndexVecVec_shared[tmp])[0] + tmpSNP_id_1st_shared);
			// print other shared SNPs from tmpChr
			for(int tmpSNPindex = 1; tmpSNPindex < (snpIndexVecVec_shared[tmp]).size(); tmpSNPindex ++)
				hap_ofs << ",SNP_" << ((snpIndexVecVec_shared[tmp])[tmpSNPindex] + tmpSNP_id_1st_shared);
			// print hap2 unique SNPs from tmpChr
			for(int tmpSNPindex = 0; tmpSNPindex < (snpIndexVecVec_unique_hap2[tmp]).size(); tmpSNPindex ++)
				hap_ofs << ",SNP_" << ((snpIndexVecVec_unique_hap2[tmp])[tmpSNPindex] + tmpSNP_id_1st_unique_hap2);
			hap_ofs << endl;
		}
		else
		{
			if((snpIndexVecVec_unique_hap1[tmp]).size() > 0)
			{	
				hap_ofs << "hap1\t" << tmpChrName << "\t1\t" << tmpChrLength-1 << "\t";
				hap_ofs << "SNP_" << ((snpIndexVecVec_unique_hap1[tmp])[0] + tmpSNP_id_1st_unique_hap1);
				for(int tmpSNPindex = 1; tmpSNPindex < (snpIndexVecVec_unique_hap1[tmp]).size(); tmpSNPindex ++)
					hap_ofs << ",SNP_" << ((snpIndexVecVec_unique_hap1[tmp])[tmpSNPindex] + tmpSNP_id_1st_unique_hap1);
				hap_ofs << endl;
			}
			if((snpIndexVecVec_unique_hap2[tmp]).size() > 0)
			{
				hap_ofs << "hap2\t" << tmpChrName << "\t1\t" << tmpChrLength-1 << "\t";
				hap_ofs << "SNP_" << ((snpIndexVecVec_unique_hap2[tmp])[0] + tmpSNP_id_1st_unique_hap2);
				for(int tmpSNPindex = 1; tmpSNPindex < (snpIndexVecVec_unique_hap2[tmp]).size(); tmpSNPindex ++)
					hap_ofs << ",SNP_" << ((snpIndexVecVec_unique_hap2[tmp])[tmpSNPindex] + tmpSNP_id_1st_unique_hap2);
				hap_ofs << endl;
			}			
		}		
	}
	hap_ofs.close();
	delete indexInfo;
	free(chrom);
	log_ofs.close();
	return 0;
}