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
	if(argc != 6)
	{
		cout << "Executable globalIndexPath inputDataSetFolderPath inputDataListFile outputFolderPath scriptFile" << endl;
		exit(1);
	}
	string mps3Sam_path_suffix = "output.sam";
	string sam2junc_path_suffix = "sam2junc";
	string globalIndexPath = argv[1];
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

	cout << "start to generate script for sam2junc cmds" << endl;
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	{
		string tmpDataSetName = dataSetNameVec[tmp];
		string tmpInputSam = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + mps3Sam_path_suffix;
		string tmpOutputSam2junc = outputFolderPath + "/" + tmpDataSetName + "/" + sam2junc_path_suffix;
		string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName;
		string tmpCmdLine = "./sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite \\\n" 
			+ globalIndexPath + " 16 \\\n" + tmpOutputSam2junc + " \\\n" + tmpInputSam;
		script_ofs << endl << tmpCmdLine << endl;
	}
	script_ofs.close();
	return 0;
}