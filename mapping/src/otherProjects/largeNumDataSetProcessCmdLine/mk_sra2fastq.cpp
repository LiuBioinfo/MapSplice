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
		cout << "#1 fastq_dump_path" << endl;
		cout << "#2 inputSRAList" << endl; 
		cout << "#3 outputCmdListFile" << endl;
		exit(1);
	}
	string fastq_dump_path = argv[1];
	string inputSRAList = argv[2];
	string outputCmdListFile = argv[3];

	vector<string> sraPathVec;
	ifstream sra_ifs(inputSRAList.c_str());
	while(!sra_ifs.eof())
	{
		string tmpStr;
		getline(sra_ifs, tmpStr);
		if(tmpStr == "")
			break;
		sraPathVec.push_back(tmpStr);
	}
	sra_ifs.close();

	ofstream cmd_ofs(outputCmdListFile.c_str());
	for(int tmp = 0; tmp < sraPathVec.size(); tmp++)
	{
		string tmpSraPath = sraPathVec[tmp];
		cmd_ofs << fastq_dump_path << " --split-3 " << tmpSraPath << endl;
	}
	cmd_ofs.close();

	return 0;
}