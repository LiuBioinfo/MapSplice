// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE

// STAR: chr1    3207050 3216714 2       2       1       1       0       67
// MPS3: chr1	3207049	3216715 	JUNC_139948	1	133	67	0	0
// STAR: chr1    4482750 4483852 2       4       1       4       0       92
// MPS3: chr1	4482749	4483853	JUNC_113477	4	168	94	0	0
// MPS3: [0] [1]-1 [2]+1 JUNC [6]+[7]	

// TOPHAT2: chr1    4009146 4014915 JUNC00000005    95      -       4009146 4014915 255,0,0 2       83,99   0,5670
// MPS3: chr1    4009229 4014817 JUNC_4  95      83      99      0       2
// MPS3: [0] [1]+[10].first [2]-[10].second+1 JUNC [4]

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

vector<string> parse_str_tab(string& tmpStr)
{
	vector<string> tmpFieldVec;
	int start_loc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmp_tab_loc = tmpStr.find("\t", start_loc);
		if(tmp_tab_loc == string::npos)
		{
			tmpFieldVec.push_back(tmpStr.substr(start_loc));
			break;
		}
		tmpFieldVec.push_back(tmpStr.substr(start_loc, tmp_tab_loc - start_loc));
		start_loc = tmp_tab_loc + 1;
		if(start_loc >= tmpStr.length())
			break;
	}
	return tmpFieldVec;
}

vector<string> parse_str_comma(string& tmpStr)
{
	vector<string> tmpFieldVec;
	int start_loc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmp_tab_loc = tmpStr.find(",", start_loc);
		if(tmp_tab_loc == string::npos)
		{
			tmpFieldVec.push_back(tmpStr.substr(start_loc));
			break;
		}
		tmpFieldVec.push_back(tmpStr.substr(start_loc, tmp_tab_loc - start_loc));
		start_loc = tmp_tab_loc + 1;
		if(start_loc >= tmpStr.length())
			break;
	}
	return tmpFieldVec;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputOtherAlignerJunc outputReformattedJunc alignerName" << endl;
		exit(1);
	}
	string inputRawJuncFileStr = argv[1];
	string outputReformattedJuncFileStr = argv[2];
	string alignerName = argv[3];

	if((alignerName != "MAPSPLICE2")&&(alignerName != "STAR")
		&&(alignerName != "STARx2")&&(alignerName != "TOPHAT2"))
	{
		cout << "invalid aligner: " << alignerName << endl;
		exit(1);
	}

	ifstream rawJunc_ifs(inputRawJuncFileStr.c_str());
	ofstream reformattedJunc_ofs(outputReformattedJuncFileStr.c_str());	
	string header;
	if(alignerName == "TOPHAT2")
		getline(rawJunc_ifs, header);
	while(!rawJunc_ifs.eof())
	{
		string tmpStr;
		getline(rawJunc_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tabVec = parse_str_tab(tmpStr);
		string chrName;
		int startPos, endPos, supNum;
		if((alignerName == "STAR")||(alignerName == "STARx2"))
		{
			chrName = tabVec[0];
			startPos = atoi(tabVec[1].c_str()) - 1;
			endPos = atoi(tabVec[2].c_str()) + 1;
			supNum = atoi(tabVec[6].c_str()) + atoi(tabVec[7].c_str()); 
		}
		else if(alignerName == "TOPHAT2")
		{
			chrName = tabVec[0];
			string anchorLengthStr = tabVec[10];
			vector<string> anchorLengthPair = parse_str_comma(anchorLengthStr);
			int anchorLength_left = atoi(anchorLengthPair[0].c_str());
			int anchorLength_right = atoi(anchorLengthPair[1].c_str());
			startPos = atoi(tabVec[1].c_str()) + anchorLength_left;
			endPos = atoi(tabVec[2].c_str()) - anchorLength_right + 1;
			supNum = atoi(tabVec[4].c_str()); 
		}
		else if(alignerName == "MAPSPLICE2")
		{
			chrName = tabVec[0];
			startPos = atoi(tabVec[1].c_str());
			endPos = atoi(tabVec[1].c_str());
			supNum = atoi(tabVec[4].c_str());
		}
		else
		{
			cout << "invalid alignerName: " << alignerName << endl;
			exit(1);
		}
		reformattedJunc_ofs << chrName << "\t" << startPos << "\t" 
			<< endPos << "\tJUNC\t" << supNum << endl; 
	}
	return 0;
}