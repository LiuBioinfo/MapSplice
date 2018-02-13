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
//#include <omp.h>
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "../../general/index_info.h"

using namespace std;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

void parseAssignmentInfoStr(string& tmpStr, int& trueId_startPos, int& trueId_endPos, int& assignedId)
{
	vector<string> tmpFieldVec;
	parseStr2fieldVec(tmpFieldVec, tmpStr);
	string trueId_startPos_str = tmpFieldVec[2];
	string trueId_endPos_str = tmpFieldVec[3];
	string assignedId_str = tmpFieldVec[4];
	trueId_startPos = atoi(trueId_startPos_str.c_str());
	trueId_endPos = atoi(trueId_endPos_str.c_str());
	assignedId = atoi(assignedId_str.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputAsssignmentInfoFile outputFolder class_num" << endl;
		exit(1);
	}
	string class_num_str = argv[3];
	int class_num = atoi(class_num_str.c_str());
	int region_num = class_num - 2;
	string inputAsssignmentInfoFile = argv[1];
	string outputFolder = argv[2];
	outputFolder += "/";
	string cmd_mkdir = "mkdir " + outputFolder;
	system(cmd_mkdir.c_str());
	string outputFolder_correct = outputFolder + "assigned_correct.txt";
	string outputFolder_incorrect = outputFolder + "assigned_incorrect.txt";
	string outputFolder_unassigned = outputFolder + "unassigned.txt";
	string outputFolder_stats = outputFolder + "stats.txt";
	ifstream raw_ifs(inputAsssignmentInfoFile.c_str());
	ofstream correct_ofs(outputFolder_correct.c_str());
	ofstream incorrect_ofs(outputFolder_incorrect.c_str());
	ofstream unassigned_ofs(outputFolder_unassigned.c_str());
	ofstream stats_ofs(outputFolder_stats.c_str());
	int correct_num = 0;
	int incorrect_num = 0;
	int unassigned_num = 0;
	int total_num = 0;
	while(!raw_ifs.eof())
	{
		string tmpStr;
		getline(raw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		total_num ++;
		int trueId_startPos, trueId_endPos, assignedId;
		parseAssignmentInfoStr(tmpStr, trueId_startPos, trueId_endPos, assignedId);
		if((assignedId < 1)||(assignedId > region_num))
		{
			unassigned_ofs << tmpStr << endl;
			unassigned_num ++;
		}
		else if((assignedId == trueId_startPos)||(assignedId == trueId_endPos))
		{
			correct_ofs << tmpStr << endl;
			correct_num ++;
		}
		else
		{
			incorrect_ofs << tmpStr << endl;
			incorrect_num ++;
		}
	}
	double correct_perc = (double)correct_num/(double)total_num;
	double incorrect_perc = (double)incorrect_num/(double)total_num;
	double assigned_perc = correct_perc + incorrect_perc;
	double unassigned_perc = (double)unassigned_num/(double)total_num;
	stats_ofs << "Total #: " << total_num << endl << "  Assigned: " << correct_num + incorrect_num << endl
		<< "    Correct #: " << correct_num << endl << "    Incorrect #: " << incorrect_num << endl 
		<< "  Unassigned #: " << unassigned_num << endl << endl;
	stats_ofs << endl << "  Assigned %: " << assigned_perc << endl << "    Correct %: " << correct_perc << endl
		<< "    Incorrect %: " << incorrect_perc << endl << "  Unassigned %: " << unassigned_perc << endl;
	stats_ofs.close();
	unassigned_ofs.close();
	incorrect_ofs.close();
	correct_ofs.close();
	raw_ifs.close();
	return 0;
}