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
	if(argc < 4)
	{
		cout << "Executable outputFolder sampleName_1 inputMPS2result_1 (sampleName_2 inputMPS2result_2 ... sampleName_N inputMPS2result_N)" << endl;
		exit(1);
	}
	string outputFolderPath = argv[1];
	cout << "creating results folder ...." << endl;
	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());
   	string logStr = outputFolderPath + "/log.txt";
   	ofstream log_ofs(logStr.c_str());

	int inputMPS2result_num = (argc - 2)/2;
	if(inputMPS2result_num * 2 + 2 != argc)
	{
		cout << "invalid parameters for sampleName_N and inputMPS2result_N " << endl;
		exit(1);
	}
	cout << "inputMPS2result_num: " << inputMPS2result_num << endl;
	log_ofs << "inputMPS2result_num: " << inputMPS2result_num << endl;
	vector<string> inputMPS2resultPathVec;
	vector<string> inputSampeNameVec;
	for(int tmp = 0; tmp < inputMPS2result_num; tmp++)
	{
		string tmpInputSampleName = argv[2+tmp*2];
		inputSampeNameVec.push_back(tmpInputSampleName);		
		string tmpInputMPS2result = argv[3+tmp*2];
		inputMPS2resultPathVec.push_back(tmpInputMPS2result);
	}
	string merge_total_file = outputFolderPath + "/mps2_fusion_merge_total.txt";
	ofstream mps2_fusion_merge_total_ofs(merge_total_file.c_str());
	for(int tmp = 0; tmp < inputMPS2result_num; tmp++)
	{
		string tmpInputMPS2result = inputMPS2resultPathVec[tmp];
		ifstream tmpMPS2fusion_ifs(tmpInputMPS2result.c_str());
		while(!tmpMPS2fusion_ifs.eof())
		{
			string tmpStr;
			getline(tmpMPS2fusion_ifs, tmpStr);
			if(tmpStr == "")
				break;
			mps2_fusion_merge_total_ofs << inputSampeNameVec[tmp] << "\t" << tmpStr << endl;
		}
		tmpMPS2fusion_ifs.close();
	}
	mps2_fusion_merge_total_ofs.close();
	log_ofs.close();
	return 0;
}