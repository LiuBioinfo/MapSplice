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

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable STARindex inputDataSetFolderPath inputDataListFile outputFolderPath scriptFile" << endl;
		exit(1);
	}
	string dataFileSuffix_1 = "1_paired.fastq";
	string dataFileSuffix_2 = "2_paired.fastq";
	string STARindexPath = argv[1];
	string inputDataSetFolderPath = argv[2];
	string inputDataListFile = argv[3];
	string outputFolderPath = argv[4];
	string scriptFile = argv[5];

	cout << "start to generate dataSetNameVec" << endl;
	vector<string> dataSetNameVec;
	ifstream dataSetList_ifs(inputDataListFile.c_str());
	while(!dataSetList_ifs.eof())
	{
		string tmpStr;
		getline(dataSetList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		dataSetNameVec.push_back(tmpStr);
	}
	dataSetList_ifs.close();

	cout << "start to generate STAR running scripts" << endl;
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	{
		string tmpDataSetName = dataSetNameVec[tmp];
		string tmpReadFilePath_1 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_1;
		string tmpReadFilePath_2 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_2;
		string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName + "/";
		string tmpCmdLine = "./STAR --runThreadN 16 --genomeDir " + STARindexPath + "/ --readFilesIn \\\n" 
			+ tmpReadFilePath_1 + " \\\n"
			+ tmpReadFilePath_2 + " \\\n" + "--outFileNamePrefix " + tmpOutputFolderPath;
		script_ofs << endl << tmpCmdLine << endl;
	}
	script_ofs.close();
	return 0;
}