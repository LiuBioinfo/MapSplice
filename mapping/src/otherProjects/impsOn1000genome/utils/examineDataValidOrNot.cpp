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

int fileSize(string& fileName)
{
    ifstream mySource;
    mySource.open(fileName.c_str(), ios_base::binary);
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
	string output_toCheck = outputDir + "toCheck.txt";
	//string outputDir_data = outputDir + "data/";
	//string mkdir_outputDir_data = "mkdir -p " + outputDir_data;
	system(mkdir_outputDir.c_str());
	//system(mkdir_outputDir_data.c_str());

	ofstream log_ofs(output_log.c_str());
	ofstream valid_ofs(output_valid.c_str());
	ofstream invalid_ofs(output_invalid.c_str());
	ofstream toCheck_ofs(output_toCheck.c_str());

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
		string tmpSNPfile_hap1_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_1.snp";
		string tmpSNPfile_hap2_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_2.snp";
		string tmpFq_1_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_1.fastq";
		string tmpFq_2_input = inputDataDir + tmpIdStr + "/" + tmpIdStr + "_2.fastq";
		int tmpSNPfileSize_1 = fileSize(tmpSNPfile_hap1_input);
		int tmpSNPfileSize_2 = fileSize(tmpSNPfile_hap2_input);
		int tmpFqFileSize_1 = fileSize(tmpFq_1_input);
		int tmpFqFileSize_2 = fileSize(tmpFq_2_input);
		log_ofs << "tmpId: " << tmpIdStr << endl;
		log_ofs << "tmpSNPfileSize_1: " << tmpSNPfileSize_1 << endl;
		log_ofs << "tmpSNPfileSize_2: " << tmpSNPfileSize_2 << endl;
		log_ofs << "tmpFqFileSize_1: " << tmpFqFileSize_1 << endl;
		log_ofs << "tmpFqFileSize_2: " << tmpFqFileSize_2 << endl;
		if((tmpSNPfileSize_1 == 0)||(tmpSNPfileSize_2 == 0)
			||(tmpFqFileSize_1 == 0)||(tmpFqFileSize_2 == 0))
			invalid_ofs << tmpIdStr << "\t" << tmpSNPfileSize_1 << "\t" << tmpSNPfileSize_2 
				<< "\t" << tmpFqFileSize_1 << "\t" << tmpFqFileSize_2 << endl;
		else if(tmpFqFileSize_1 == tmpFqFileSize_2)
			valid_ofs << tmpIdStr << "\t" << tmpSNPfileSize_1 << "\t" << tmpSNPfileSize_2 
				<< "\t" << tmpFqFileSize_1 << "\t" << tmpFqFileSize_2 << endl;
		else
			toCheck_ofs << tmpIdStr << "\t" << tmpSNPfileSize_1 << "\t" << tmpSNPfileSize_2 
				<< "\t" << tmpFqFileSize_1 << "\t" << tmpFqFileSize_2 << endl;
	}
	idList_ifs.close();

	log_ofs.close();
	valid_ofs.close();
	invalid_ofs.close();
	toCheck_ofs.close();
	return 0;
}