// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
// input-1: original wget cmd 
// wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060890/SRR060890.sra
// wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060891/SRR060891.sra
// wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060892/SRR060892.sra
// wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060893/SRR060893.sra
// input-2: already downloaded sra
// SRR067891.sra
// SRR067893.sra
// SRR067895.sra


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

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_originalWgetCmdListFile" << endl;
		cout << "#2 input_alreadyDownloadedSraListFile" << endl;
		cout << "#3 output_updatedWgetCmdListFile" << endl;
		exit(1);
	}

	string input_originalWgetCmdListFile = argv[1];
	string input_alreadyDownloadedSraListFile = argv[2];
	string output_updatedWgetCmdListFile = argv[3];

	vector<string> rawWgetCmdVec;
	vector<string> rawSRRvec;
	ifstream originalWgetCmdList_ifs(input_originalWgetCmdListFile.c_str());
	while(!originalWgetCmdList_ifs.eof())
	{
		string tmpStr;
		getline(originalWgetCmdList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int SRR_loc_1 = tmpStr.find("SRR");
		int SRR_loc_2 = tmpStr.find("SRR", SRR_loc_1 + 1);
		int SRR_loc_3 = tmpStr.find("SRR", SRR_loc_2 + 1); 
		int SRR_loc_4 = tmpStr.find("SRR", SRR_loc_3 + 1); 
		string tmpRawSRR = tmpStr.substr(SRR_loc_4);
		cout << "tmpRawWget: " << tmpStr << endl;
		cout << "tmpRawSRR: " << tmpRawSRR << endl;
		rawWgetCmdVec.push_back(tmpStr);
		rawSRRvec.push_back(tmpRawSRR);
	}
	originalWgetCmdList_ifs.close();

	vector<string> downloadedSRRvec;
	ifstream alreadyDownloadedSraList_ifs(input_alreadyDownloadedSraListFile.c_str());
	while(!alreadyDownloadedSraList_ifs.eof())
	{
		string tmpStr;
		getline(alreadyDownloadedSraList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		cout << "tmpDownloadedSRR: " << tmpStr << endl;
		downloadedSRRvec.push_back(tmpStr);
	}
	alreadyDownloadedSraList_ifs.close();

	ofstream updatedWgetCmdList_ofs(output_updatedWgetCmdListFile.c_str());
	for(int tmp = 0; tmp < rawSRRvec.size(); tmp++)
	{
		string tmpRawSRR = rawSRRvec[tmp];
		bool already_downloaded_bool = false;
		for(int tmp2 = 0; tmp2 < downloadedSRRvec.size(); tmp2++)
		{
			if(tmpRawSRR == downloadedSRRvec[tmp2])
			{
				already_downloaded_bool = true;
				break;
			}
		}
		if(!already_downloaded_bool)
			updatedWgetCmdList_ofs << rawWgetCmdVec[tmp] << endl;
	}
	updatedWgetCmdList_ofs.close();
	return 0;
}