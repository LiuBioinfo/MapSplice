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

	cout << "start to generate mkdir_script" << endl;
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	{
		string tmpDataSetName = dataSetNameVec[tmp];
		string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName;
		string tmpDir_snpMerLearned = tmpOutputFolderPath + "/SNPmer_learned";
		string tmpDir_phase1 = tmpOutputFolderPath + "/phase1_output";
		string tmpDir_phase1_repeat_region = tmpOutputFolderPath + "/phase1_output/repeat_region";
		string tmpDir_phase1_oneEndUnmapped = tmpOutputFolderPath + "/phase1_output/oneEndUnmapped";
		string tmpDir_phase1_incomplete = tmpOutputFolderPath + "/phase1_output/incomplete";
		string tmpDir_phase1_completePair = tmpOutputFolderPath + "/phase1_output/completePair";
		string tmpDir_phase1_bothEndsUnmapped = tmpOutputFolderPath + "/phase1_output/bothEndsUnmapped";

		string tmpDir_phase2 = tmpOutputFolderPath + "/phase2_output";
		string tmpDir_logs = tmpOutputFolderPath + "/logs";
		string tmpDir_logs_phase1 = tmpOutputFolderPath + "/logs/phase1_log";
		string tmpDir_logs_phase2 = tmpOutputFolderPath + "/logs/phase2_log";
		string tmpDir_logs_phase2_fixHeadTail = tmpOutputFolderPath + "/logs/phase2_log/fixHeadTail";
		string tmpDir_logs_phase2_fixOneEndUnmapped = tmpOutputFolderPath + "/logs/phase2_log/fixOneEndUnmapped";

		script_ofs << endl;
		script_ofs << "mkdir " << tmpOutputFolderPath << endl;
		script_ofs << "mkdir " << tmpDir_snpMerLearned << endl;
		script_ofs << "mkdir " << tmpDir_phase1 << endl;
		script_ofs << "mkdir " << tmpDir_phase1_repeat_region << endl;
		script_ofs << "mkdir " << tmpDir_phase1_oneEndUnmapped << endl;
		script_ofs << "mkdir " << tmpDir_phase1_incomplete << endl;
		script_ofs << "mkdir " << tmpDir_phase1_completePair << endl;
		script_ofs << "mkdir " << tmpDir_phase1_bothEndsUnmapped << endl;
		script_ofs << "mkdir " << tmpDir_phase2 << endl;
		script_ofs << "mkdir " << tmpDir_logs << endl;
		script_ofs << "mkdir " << tmpDir_logs_phase1 << endl;
		script_ofs << "mkdir " << tmpDir_logs_phase2 << endl;
		script_ofs << "mkdir " << tmpDir_logs_phase2_fixHeadTail << endl;
		script_ofs << "mkdir " << tmpDir_logs_phase2_fixOneEndUnmapped << endl;
	}
	script_ofs.close();
	return 0;
}