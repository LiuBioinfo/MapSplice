// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
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

void getSRRidFromBothPath(string& tmpSRRid_1, string& tmpSRRid_2, 
	string& tmpFqPath_1, string& tmpFqPath_2)
{
	int SRRloc_start_1 = tmpFqPath_1.rfind("SRR");
	int SRRloc_end_1 = tmpFqPath_1.rfind("_1.fastq");
	int SRRloc_start_2 = tmpFqPath_2.rfind("SRR");
	int SRRloc_end_2 = tmpFqPath_2.rfind("_2.fastq");	
	tmpSRRid_1 = tmpFqPath_1.substr(SRRloc_start_1, SRRloc_end_1 - SRRloc_start_1);
	tmpSRRid_2 = tmpFqPath_2.substr(SRRloc_start_2, SRRloc_end_2 - SRRloc_start_2);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputFastqFileList" << endl;
		cout << "#2 outputMergedFqDir" << endl;
		cout << "#3 outputCatCmdListFile" << endl;
		exit(1);
	}
	string inputFastqFileList = argv[1];
	string outputMergedFqDir = argv[2]; outputMergedFqDir += "/";
	string outputCatCmdListFile = argv[3];

	vector<string> fqPathVec;
	ifstream fq_ifs(inputFastqFileList.c_str());
	while(!fq_ifs.eof())
	{
		string tmpStr;
		getline(fq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		fqPathVec.push_back(tmpStr);
	}
	fq_ifs.close();

	vector<string> SRRidVec;
	int sample_num = fqPathVec.size()/2;
	for(int tmp = 0; tmp < sample_num; tmp++)
	{
		string tmpFqPath_1 = fqPathVec[tmp * 2];
		string tmpFqPath_2 = fqPathVec[tmp * 2 + 1];
		string tmpSRRid_1, tmpSRRid_2;
		getSRRidFromBothPath(tmpSRRid_1, tmpSRRid_2, tmpFqPath_1, tmpFqPath_2);
		if(tmpSRRid_1 != tmpSRRid_2)
		{
			cout << "tmpSRRid_1 != tmpSRRid_2" << endl;
			exit(1);
		}
		else
			SRRidVec.push_back(tmpSRRid_1);
	}

	ofstream cmd_ofs(outputCatCmdListFile.c_str());
	for(int tmp = 0; tmp < sample_num; tmp++)
	{
		string tmpFqPath_1 = fqPathVec[tmp * 2];
		string tmpFqPath_2 = fqPathVec[tmp * 2 + 1];	
		cmd_ofs << "cat \\" << endl << tmpFqPath_1 << " \\" << endl
			<< tmpFqPath_2 << " > \\" << endl << outputMergedFqDir 
			<< SRRidVec[tmp] << ".fastq" << endl << endl;
		cmd_ofs << "rm " << tmpFqPath_1 << endl << endl 
			<< "rm " << tmpFqPath_2 << endl << endl;
	}
	cmd_ofs.close();
	return 0;
}