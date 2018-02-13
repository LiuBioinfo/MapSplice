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

using namespace std;

int fileSize(const char *add)
{
    ifstream mySource;
    mySource.open(add, ios_base::binary);
    mySource.seekg(0,ios_base::end);
    int size = mySource.tellg();
    mySource.close();
    return size;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputDataIdListFile" << endl;
		cout << "#2 inputDataDir" << endl;
		cout << "#3 outputDir" << endl;
		exit(1);
	}
	string inputDataIdListFile = argv[1];
	string inputDataDir = argv[2]; inputDataDir += "/";
	string outputDir = argv[3]; outputDir += "/";

	string mkdir_outputDir = "mkdir -p " + outputDir;
	string output_log = outputDir + "log.txt";
	string output_valid = outputDir + "valid.txt";
	string output_invalid = outputDir + "invalid.txt";
	string outputDir_data = outputDir + "data/";
	string mkdir_outputDir_data = "mkdir -p " + outputDir_data;
	system(mkdir_outputDir.c_str());
	system(mkdir_outputDir_data.c_str());

	ofstream log_ofs(output_log.c_str());
	ofstream valid_ofs(output_valid.c_str());
	ofstream invalid_ofs(output_invalid.c_str());

	cout << "start to process each individual's datasets" << endl;
	log_ofs << "start to process each individual's datasets" << endl;
	ifstream idList_ifs(inputDataIdListFile.c_str());
	int id_num = 0;
	while(!idList_ifs.eof())
	{
		string tmpIdStr;
		getline(idList_ifs, tmpIdStr);
		if(tmpIdStr == "")
			break;
		id_num ++;
		//if(id_num >= 10)
		//	exit(1);
		cout << endl << "tmpId:\t" << tmpIdStr << endl;
		log_ofs << endl << "tmpId:\t" << tmpIdStr << endl;
		log_ofs << "tmpId_mum:\t" << id_num << endl;
		string tmpSNPfile_hap1_input = inputDataDir + tmpIdStr + "/SNPFiles/genome1/totalSNP.snp";
		string tmpSNPfile_hap2_input = inputDataDir + tmpIdStr + "/SNPFiles/genome2/totalSNP.snp"; 
		string tmpFq_1_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_1.fastq.gz";
		string tmpFq_2_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_2.fastq.gz";

		string tmpIdDir_output = outputDir_data + tmpIdStr + "/";
		string cmd_mkdir_tmpIdDir_output = "mkdir -p " + tmpIdDir_output;
		string cmd_cp_SNP_hap1 = "cp " + tmpSNPfile_hap1_input + " " + tmpIdDir_output + tmpIdStr + "_1.snp";
		string cmd_cp_SNP_hap2 = "cp " + tmpSNPfile_hap2_input + " " + tmpIdDir_output + tmpIdStr + "_2.snp";
		string cmd_gunzip_fq1 = "gunzip -c " + tmpFq_1_input + " > " + tmpIdDir_output + tmpIdStr + "_1.fastq";
		string cmd_gunzip_fq2 = "gunzip -c " + tmpFq_2_input + " > " + tmpIdDir_output + tmpIdStr + "_2.fastq";
		log_ofs << "cmd_mkdir_tmpIdDir_output: " << cmd_mkdir_tmpIdDir_output << endl;
		log_ofs << "cmd_cp_SNP_hap1: " << cmd_cp_SNP_hap1 << endl;
		log_ofs << "cmd_cp_SNP_hap2: " << cmd_cp_SNP_hap2 << endl;
		log_ofs << "cmd_gunzip_fq1: " << cmd_gunzip_fq1 << endl;
		log_ofs << "cmd_gunzip_fq2: " << cmd_gunzip_fq2 << endl;

		const int tmp_mkdir_err = system(cmd_mkdir_tmpIdDir_output.c_str());
		if(tmp_mkdir_err == -1)
			log_ofs << "error in mkdir" << endl;
		const int tmp_cp_snp_hap1_err = system(cmd_cp_SNP_hap1.c_str());
		if(tmp_cp_snp_hap1_err == -1)
			log_ofs << "error in cp_snp_hap1" << endl;
		const int tmp_cp_snp_hap2_err = system(cmd_cp_SNP_hap2.c_str());
		if(tmp_cp_snp_hap2_err == -1)
			log_ofs << "error in cp_snp_hap2" << endl;
		const int tmp_gunzip_fq1_err = system(cmd_gunzip_fq1.c_str());
		if(tmp_gunzip_fq1_err == -1)
			log_ofs << "error in gunzip_fq1" << endl;
		const int tmp_gunzip_fq2_err = system(cmd_gunzip_fq2.c_str());
		if(tmp_gunzip_fq2_err == -1)
			log_ofs << "error in gunzip_fq2" << endl;
		if((tmp_mkdir_err == -1)||(tmp_cp_snp_hap1_err == -1)||(tmp_cp_snp_hap2_err == -1)
			||(tmp_gunzip_fq1_err == -1)||(tmp_gunzip_fq2_err == -1))
			invalid_ofs << tmpIdStr << endl;
		else
			valid_ofs << tmpIdStr << endl;
	}
	idList_ifs.close();

	log_ofs.close();
	valid_ofs.close();
	invalid_ofs.close();
	return 0;
}