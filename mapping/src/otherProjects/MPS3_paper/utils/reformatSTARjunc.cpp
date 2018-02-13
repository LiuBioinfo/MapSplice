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
	if(argc != 3)
	{
		cout << "Executable inputSTARjunc outputReformattedJunc" << endl;
		exit(1);
	}
	string inputSTARjuncFileStr = argv[1];
	string outputReformattedJuncFileStr = argv[2];
	ifstream STARjunc_ifs(inputSTARjuncFileStr.c_str());
	ofstream reformattedJunc_ofs(outputReformattedJuncFileStr.c_str());
	int tmpJunc_IDnum = 0;
	while(!STARjunc_ifs.eof())
	{
		string tmpStr;
		getline(STARjunc_ifs, tmpStr);
		if((STARjunc_ifs.eof())||(tmpStr == ""))
			break;
		tmpJunc_IDnum ++;
		vector<string> tmpSTARjuncFieldVec;
		int tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			if(tmpTabLoc != string::npos)
			{	
				string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
				tmpSTARjuncFieldVec.push_back(tmpFieldStr);
				tmpStartLoc = tmpTabLoc + 1;
			}
			else
			{	
				string tmpLastFieldStr = tmpStr.substr(tmpStartLoc);
				tmpSTARjuncFieldVec.push_back(tmpLastFieldStr);
				break;
			}
		}
		//cout << "tmpSTARjuncFieldVec.size(): " << tmpSTARjuncFieldVec.size() << endl;
		string tmpJunc_chrName = tmpSTARjuncFieldVec[0];
		string tmpJunc_1stPosInIntronStr = tmpSTARjuncFieldVec[1];
		string tmpJunc_2ndPosInIntronStr = tmpSTARjuncFieldVec[2];
		string tmpJunc_uniMapReadCountStr = tmpSTARjuncFieldVec[6];
		string tmpJunc_mulMapReadCountStr = tmpSTARjuncFieldVec[7];
		int tmpJunc_1stPosInIntron = atoi(tmpJunc_1stPosInIntronStr.c_str());
		int tmpJunc_endPosInIntron = atoi(tmpJunc_2ndPosInIntronStr.c_str());
		int tmpJunc_uniMapReadCount = atoi(tmpJunc_uniMapReadCountStr.c_str());
		int tmpJunc_mulMapReadCount = atoi(tmpJunc_mulMapReadCountStr.c_str()); 
		int tmpJunc_totalMapReadCount = tmpJunc_uniMapReadCount + tmpJunc_mulMapReadCount;
		int tmpJunc_startPos = tmpJunc_1stPosInIntron - 1;
		int tmpJunc_endPos = tmpJunc_endPosInIntron + 1;
		reformattedJunc_ofs << tmpJunc_chrName << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos 
			<< "\tJUNC_" << tmpJunc_IDnum << "\t" << tmpJunc_totalMapReadCount << endl; 
	}
	STARjunc_ifs.close();
	reformattedJunc_ofs.close();
	return 0;
}