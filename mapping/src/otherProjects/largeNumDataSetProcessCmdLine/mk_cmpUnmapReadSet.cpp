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
	if(argc != 5)
	{
		cout << "Executable inputFolder outputFolder dataNum scriptFile" << endl;
		exit(1);
	}

	string inputFolderPath = argv[1];
	string outputFolderPath = argv[2];
	string dataNum_str = argv[3];
	int data_num = atoi(dataNum_str.c_str());
	string scriptFile = argv[4];
	ofstream script_ofs(scriptFile.c_str());
	for(int tmp1 = 1; tmp1 <= data_num; tmp1++)
	{
		for(int tmp2 = 1; tmp2 <= data_num; tmp2++)
		{
			if(tmp1 == tmp2)
				continue;
			string exeName = "./cmpDataKmerAgainstBloomFilter \\\n";
			script_ofs << exeName << inputFolderPath <<  "/" << tmp1 << "/unmappedRead.merged.fq.jf.fa \\\n"
				<< inputFolderPath <<  "/" << tmp2 << "/unmappedRead.merged.fq.jf.fa \\\n" << outputFolderPath 
				<< "/against_" << tmp1 << "/" << tmp2 << endl << endl;
		}
	}
	script_ofs.close();

	// cout << "start to generate dataSetNameVec" << endl;
	// vector<string> dataSetNameVec;
	// ifstream dataSetList_ifs(inputDataListFile.c_str());
	// while(!dataSetList_ifs.eof())
	// {
	// 	string tmpStr;
	// 	getline(dataSetList_ifs, tmpStr);
	// 	if(tmpStr == "")
	// 		break;
	// 	dataSetNameVec.push_back(tmpStr);
	// }
	// dataSetList_ifs.close();

	// cout << "start to generate STAR running scripts" << endl;
	// ofstream script_ofs(scriptFile.c_str());
	// for(int tmp = 0; tmp < dataSetNameVec.size(); tmp++)
	// {
	// 	string tmpDataSetName = dataSetNameVec[tmp];
	// 	string tmpReadFilePath_1 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_1;
	// 	string tmpReadFilePath_2 = inputDataSetFolderPath + "/" + tmpDataSetName + "/" + dataFileSuffix_2;
	// 	string tmpOutputFolderPath = outputFolderPath + "/" + tmpDataSetName + "/";
	// 	string tmpCmdLine = "./STAR --runThreadN 16 --genomeDir " + STARindexPath + "/ --readFilesIn \\\n" 
	// 		+ tmpReadFilePath_1 + " \\\n"
	// 		//+ tmpReadFilePath_2 + " \\\n" 
	// 		+ "--outFileNamePrefix " + tmpOutputFolderPath;
	// 	script_ofs << endl << tmpCmdLine << endl;
	// }
	script_ofs.close();
	return 0;
}