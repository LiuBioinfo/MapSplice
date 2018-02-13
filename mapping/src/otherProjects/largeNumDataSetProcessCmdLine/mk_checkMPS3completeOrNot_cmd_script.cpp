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
		cout << "Executable inputDataListFile outputFolderPath scriptFile" << endl;
		exit(1);
	}
	string inputDataListFile = argv[1];
	string outputFolderPath = argv[2];
	string scriptFile = argv[3];

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

	cout << "start to generate checkMPS3completeOrNot_script" << endl;
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	{
		string tmpDataSetName = dataSetNameVec[tmp];
		string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName;
		string tmpFile_stats = tmpOutputFolderPath + "/stats.txt";

		//script_ofs << endl;
		script_ofs << "wc -l " << tmpFile_stats << endl;
	}
	script_ofs.close();
	return 0;
}