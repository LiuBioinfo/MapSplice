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

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable globalIndexPath localIndexPath ";
		cout << "inputDataSetFolderPath inputDataListFile outputFolderPath scriptFile" << endl;
		exit(1);
	}
	string dataFileSuffix_1 = "1_paired.fastq";
	string dataFileSuffix_2 = "2_paired.fastq";
	string globalIndexPath = argv[1];
	string localIndexPath = argv[2];
	string inputDataSetFolderPath = argv[3];
	string inputDataListFile = argv[4];
	string outputFolderPath = argv[5];
	string scriptFile = argv[6];

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

	cout << "start to generate read File Path Vec and output Folder Path Vec" << endl;
	// vector<string> readFilePathVec_1;
	// vector<string> readFilePathVec_2;
	// vector<string> outputFolderPathVec;
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	{
		string tmpDataSetName = dataSetNameVec[tmp];
		string tmpReadFilePath_1 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_1;
		string tmpReadFilePath_2 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_2;
		string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName;
		string tmpCmdLine = "./mps -J -B -T 16 -G " + globalIndexPath + " \\\n" 
			+ "-L" + localIndexPath + " \\\n" + "-1 " + tmpReadFilePath_1 + " \\\n"
			+ "-2 " + tmpReadFilePath_2 + " \\\n" + "-O " + tmpOutputFolderPath;
		script_ofs << endl << tmpCmdLine << endl;
	}
	script_ofs.close();
	return 0;
}