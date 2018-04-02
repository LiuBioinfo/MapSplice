/*
FILE 1: 
inputMergedChrFa:
>id1xxxxx
AAACATAAGG
>id2xxxxx
ATGCACAACA
ATGCACAACA

FILE 2:
inputChrNameList:
id1
id2

FILE 3:
inputChrNameLineNumList:
1
3
6

Note: do not forget the last line
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <vector>
#include <map>
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputMergedChrFa" << endl;
		// /home/xli262/cab3/eca_ref_EquCab3.0_chrUn.fa
		cout << "#2 inputChrNameList" << endl;
		// /home/xli262/cab3/eca_ref_EquCab3.0_chrUn.fa.chrNameList
		cout << "#3 inputChrNameLineNumList" << endl; 
		// chrNameLineNumList.size() = chrNameList.size() + 1
		// the last num val = mergedChrFaLineNum + 1
		// /home/xli262/cab3/eca_ref_EquCab3.0_chrUn.fa.chrNameLineNum
		cout << "#4 outputDir" << endl;
		exit(1);
	}

	string inputMergedChrFa = argv[1];
	string inputChrNameList = argv[2];
	string inputChrNameLineNumList = argv[3];
	string outputDir = argv[4]; outputDir += "/";

	cout << "loading inputMergedChrFa ..." << endl;
	vector<string> faStrVec;
	ifstream fa_ifs(inputMergedChrFa.c_str());
	while(!fa_ifs.eof())
	{
		string tmpStr;
		getline(fa_ifs, tmpStr);
		if(tmpStr == "")
			break;
		faStrVec.push_back(tmpStr);
	}
	fa_ifs.close();

	cout << "loading chrNameList ..." << endl;
	vector<string> chrNameVec;
	ifstream name_ifs(inputChrNameList.c_str());
	while(!name_ifs.eof())
	{
		string tmpStr;
		getline(name_ifs, tmpStr);
		if(tmpStr == "")
			break;
		chrNameVec.push_back(tmpStr);
	}
	name_ifs.close();
	
	cout << "loading chrNameLineNumList ..." << endl;
	vector<int> lineNumVec;
	ifstream lineNum_ifs(inputChrNameLineNumList.c_str());
	while(!lineNum_ifs.eof())
	{
		string tmpStr;
		getline(lineNum_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpLineNum = atoi(tmpStr.c_str());
		lineNumVec.push_back(tmpLineNum);
	}
	lineNum_ifs.close();

	cout << "printing chrFa files ..." << endl;
	for(int tmpChrIndex = 0; tmpChrIndex < chrNameVec.size(); tmpChrIndex++)
	{
		string tmpChrName = chrNameVec[tmpChrIndex];
		string tmpChrFilePath = outputDir + tmpChrName + ".fa";
		ofstream tmpFa_ofs(tmpChrFilePath.c_str());
		tmpFa_ofs << ">" << tmpChrName << endl;
		int lineIndexInMergedFa_start = lineNumVec[tmpChrIndex];
		int lineIndexInMergedFa_end = lineNumVec[tmpChrIndex + 1] - 2;
		for(int tmpLineIndex = lineIndexInMergedFa_start;
			tmpLineIndex <= lineIndexInMergedFa_end; tmpLineIndex++)
			tmpFa_ofs << faStrVec[tmpLineIndex] << endl;
		tmpFa_ofs.close();
	}
	cout << "All jobs done!" << endl;
	return 0;
}