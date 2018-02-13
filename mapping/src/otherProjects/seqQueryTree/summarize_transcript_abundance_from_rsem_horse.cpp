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
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

void generateRnaIdVec(vector<string>& rnaIdVec, 
	string& rnaId_transcriptId_map_file)
{
	ifstream rnaId_ifs(rnaId_transcriptId_map_file.c_str());
	while(!rnaId_ifs.eof())
	{
		string tmpStr;
		getline(rnaId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpRnaId = tmpStr.substr(0, tabLoc);
		rnaIdVec.push_back(tmpRnaId);
	}
	rnaId_ifs.close();
}

void generateRsemFileVec(vector<string>& rsemFileVec, string& rsemFileList)
{
	ifstream list_ifs(rsemFileList.c_str());
	while(!list_ifs.eof())
	{
		string tmpStr;
		getline(list_ifs, tmpStr);
		if(tmpStr == "")
			break;
		rsemFileVec.push_back(tmpStr);
	}
	list_ifs.close();
}

void getRsemIdVec(vector<string>& rsemIdVec, string& tmpRsemFile)
{
	ifstream rsem_ifs(tmpRsemFile.c_str());
	string tmpHeader;
	getline(rsem_ifs, tmpHeader);
	while(!rsem_ifs.eof())
	{
		string tmpStr;
		getline(rsem_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		string tmpRsemRnaId = tmpStr.substr(0, tabLoc_1);
		rsemIdVec.push_back(tmpRsemRnaId);
	}
	rsem_ifs.close();		
}

void generateRnaAbundanceVec(vector<double>& expectedCountVec,
	vector<double>& tpmVec, vector<double>& fpkmVec, string& tmpRsemFile)
{
	ifstream rsem_ifs(tmpRsemFile.c_str());
	string tmpHeader;
	getline(rsem_ifs, tmpHeader);
	while(!rsem_ifs.eof())
	{
		string tmpStr;
		getline(rsem_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
		int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
		string tmpEffectedCountStr = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string tmpTpmStr = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		string tmpFpkmStr = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
		double tmpEffectedCount = atof(tmpEffectedCountStr.c_str());
		double tmpTpm = atof(tmpTpmStr.c_str());
		double tmpFpkm = atof(tmpFpkmStr.c_str());
		expectedCountVec.push_back(tmpEffectedCount);
		tpmVec.push_back(tmpTpm);
		fpkmVec.push_back(tmpFpkm);
	}
	rsem_ifs.close();
}

bool searchRsemRnaIdInRnaIdVec(string& tmpRsemRnaId, vector<string>& rnaIdVec)
{
	for(int tmp = 0; tmp < rnaIdVec.size(); tmp++)
	{
		if(tmpRsemRnaId == rnaIdVec[tmp])
			return true;
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_rsemFileList" << endl; // /scratch/xli262/seqQueryTree/rsem/batch_2/1.isoforms.results
		cout << "#2 input_rnaId_transcriptId_map" << endl; // /scratch/xli262/seqQueryTree/rnaId_transcriptId.map
		cout << "#3 output_file_prefix" << endl;
		exit(1);
	}

	string input_rsemFileList = argv[1];
	string input_rnaId_transcriptId_map = argv[2];
	string output_file_prefix = argv[3];
	string output_expectedCount = output_file_prefix + ".expectedCount";
	string output_TPM = output_file_prefix + ".TPM";
	string output_FPKM = output_file_prefix + ".FPKM";

	vector<string> rnaIdVec;
	generateRnaIdVec(rnaIdVec, input_rnaId_transcriptId_map);
	vector<string> rsemFileVec;
	generateRsemFileVec(rsemFileVec, input_rsemFileList);
	vector< vector<double> > effectedCountVecVec;
	vector< vector<double> > tpmVecVec;
	vector< vector<double> > fpkmVecVec;
	for(int tmp = 0; tmp < rsemFileVec.size(); tmp++)
	{
		string tmpRsemFile = rsemFileVec[tmp];
		vector<double> expectedCountVec; 
		vector<double> tpmVec; 
		vector<double> fpkmVec;
		generateRnaAbundanceVec(expectedCountVec, tpmVec, fpkmVec, tmpRsemFile);
		effectedCountVecVec.push_back(expectedCountVec);
		tpmVecVec.push_back(tpmVec);
		fpkmVecVec.push_back(fpkmVec);
	}

	vector<string> rsemRnaIdVec;
	getRsemIdVec(rsemRnaIdVec, rsemFileVec[0]);

	ofstream expectedCount_ofs(output_expectedCount.c_str());
	ofstream tpm_ofs(output_TPM.c_str());
	ofstream fpkm_ofs(output_FPKM.c_str());
	expectedCount_ofs << "rna";
	tpm_ofs << "rna";
	fpkm_ofs << "rna";
	for(int tmp = 0; tmp < rsemFileVec.size(); tmp++)
	{
		expectedCount_ofs << "\t" << tmp + 1;
		tpm_ofs << "\t" << tmp + 1;
		fpkm_ofs << "\t" << tmp + 1;
	}
	expectedCount_ofs << endl;
	tpm_ofs << endl;
	fpkm_ofs << endl;	
	for(int tmp = 0; tmp < rsemRnaIdVec.size(); tmp++)
	{
		string tmpRsemRnaId = rsemRnaIdVec[tmp];
		bool exists_in_fa_bool = searchRsemRnaIdInRnaIdVec(tmpRsemRnaId, rnaIdVec);
		if(exists_in_fa_bool)
		{
			expectedCount_ofs << tmpRsemRnaId;
			tpm_ofs << tmpRsemRnaId;
			fpkm_ofs << tmpRsemRnaId;
			for(int tmp2 = 0; tmp2 < rsemFileVec.size(); tmp2 ++)
			{
				expectedCount_ofs << "\t" << (effectedCountVecVec[tmp2])[tmp];
				tpm_ofs << "\t" << (tpmVecVec[tmp2])[tmp];
				fpkm_ofs << "\t" << (fpkmVecVec[tmp2])[tmp];				
			}			
			expectedCount_ofs << endl;
			tpm_ofs << endl;
			fpkm_ofs << endl;
		}
	}

	expectedCount_ofs.close();
	tpm_ofs.close();
	fpkm_ofs.close();
	return 0;
}