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
		cout << "Executable inputTophat2junc outputReformattedJunc" << endl;
		exit(1);
	}
	string inputTophat2juncFileStr = argv[1];
	string outputReformattedJuncFileStr = argv[2];
	ifstream Tophat2junc_ifs(inputTophat2juncFileStr.c_str());
	ofstream reformattedJunc_ofs(outputReformattedJuncFileStr.c_str());
	int tmpJunc_IDnum = 0;
	string tmp1stLine;
	getline(Tophat2junc_ifs, tmp1stLine);
	while(!Tophat2junc_ifs.eof())
	{
		string tmpStr;
		getline(Tophat2junc_ifs, tmpStr);
		if((Tophat2junc_ifs.eof())||(tmpStr == ""))
			break;
		tmpJunc_IDnum ++;
		vector<string> tmpTophat2juncFieldVec;
		int tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			if(tmpTabLoc != string::npos)
			{	
				string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
				tmpTophat2juncFieldVec.push_back(tmpFieldStr);
				tmpStartLoc = tmpTabLoc + 1;
			}
			else
			{	
				string tmpLastFieldStr = tmpStr.substr(tmpStartLoc);
				tmpTophat2juncFieldVec.push_back(tmpLastFieldStr);
				break;
			}
		}
		//cout << "tmpTophat2juncFieldVec.size(): " << tmpTophat2juncFieldVec.size() << endl;
		string tmpJunc_chrName = tmpTophat2juncFieldVec[0];
		string tmpJunc_startPosIn1stBlockStr = tmpTophat2juncFieldVec[1];
		string tmpJunc_endPosIn2ndBlockStr = tmpTophat2juncFieldVec[2];
		string tmpJunc_supNumStr = tmpTophat2juncFieldVec[4];
		string tmpJunc_twoAnchorLenStr = tmpTophat2juncFieldVec[10];
		//cout << "tmpJunc_supNumStr: " << tmpJunc_supNumStr << endl;
		//cout << "tmpJunc_twoAnchorLenStr: " << tmpJunc_twoAnchorLenStr << endl;

		int tmpJunc_startPosIn1stBlock = atoi(tmpJunc_startPosIn1stBlockStr.c_str());
		int tmpJunc_endPosIn2ndBlock = atoi(tmpJunc_endPosIn2ndBlockStr.c_str());
		int tmpJunc_totalMapReadCount = atoi(tmpJunc_supNumStr.c_str());
		int tmpCommaLoc = tmpJunc_twoAnchorLenStr.find(",");
		string tmpJunc_1stAnchorLenStr = tmpJunc_twoAnchorLenStr.substr(0, tmpCommaLoc);
		string tmpJunc_2ndAnchorLenStr = tmpJunc_twoAnchorLenStr.substr(tmpCommaLoc + 1);
		int tmpJunc_1stAnchorLen = atoi(tmpJunc_1stAnchorLenStr.c_str());
		int tmpJunc_2ndAnchorLen = atoi(tmpJunc_2ndAnchorLenStr.c_str());

		int tmpJunc_startPos = tmpJunc_startPosIn1stBlock + tmpJunc_1stAnchorLen;
		int tmpJunc_endPos = tmpJunc_endPosIn2ndBlock - tmpJunc_2ndAnchorLen + 1;

		reformattedJunc_ofs << tmpJunc_chrName << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos 
			<< "\tJUNC_" << tmpJunc_IDnum << "\t" << tmpJunc_totalMapReadCount << endl; 
	}
	Tophat2junc_ifs.close();
	reformattedJunc_ofs.close();
	return 0;
}