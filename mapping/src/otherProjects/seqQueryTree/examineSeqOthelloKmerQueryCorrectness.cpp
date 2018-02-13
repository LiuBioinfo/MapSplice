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

void get_resultsFileVec(string& results_file_list_file, vector<string>& resultsFile_vec)
{
	ifstream results_list_ifs(results_file_list_file.c_str());
	while(!results_list_ifs.eof())
	{
		string tmpStr;
		getline(results_list_ifs, tmpStr);
		if(tmpStr == "")
			break;
		resultsFile_vec.push_back(tmpStr);
	}
	results_list_ifs.close();
}

void loadQueryResults2seqVecKmerVecVec(
	string& input_seqOthelloResults_file, vector<string>& seqVec, 
	vector< vector<int> >& queryVecVec)
{
	ifstream seqOtehllo_ifs(input_seqOthelloResults_file.c_str());
	int tmpQueryNum = 0;
	while(!seqOtehllo_ifs.eof())
	{
		string tmpStr;
		getline(seqOtehllo_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpQueryNum ++;
		if((tmpQueryNum / 10000) * 10000 == tmpQueryNum)
			cout << "tmpQueryNum: " << tmpQueryNum << endl;		
		//int tmpSeqLength = tmpStr.length();
		seqVec.push_back(tmpStr);
		//seqLengthVec.push_back(tmpSeqLength);
		vector<int> tmpKmerCountVec;
		for(int tmp = 0; tmp < 2592; tmp++)
		{
			string tmpStr;
			getline(seqOtehllo_ifs, tmpStr);
			int tmpSpaceLoc = tmpStr.find(" ");
			string tmpKmerCountStr = tmpStr.substr(tmpSpaceLoc + 1);
			int tmpKmerCount = atoi(tmpKmerCountStr.c_str());
			tmpKmerCountVec.push_back(tmpKmerCount);
		}
		queryVecVec.push_back(tmpKmerCountVec);
	}
	seqOtehllo_ifs.close();
}


int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 kmerFileList" << endl;
		cout << "#2 outputDir" << endl;
		cout << "#3 idList" << endl;
		exit(1);
	}
	string kmerFileList = argv[1];
	string outputDir = argv[2]; outputDir += "/";
	string idList = argv[3];
	string error_file = outputDir + "err.txt";

	vector<string> kmerFileVec;
	vector<string> idVec;
	get_resultsFileVec(kmerFileList, kmerFileVec);
	get_resultsFileVec(idList, idVec);
	ofstream err_ofs(error_file.c_str());
	for(int tmp = 0; tmp < kmerFileVec.size(); tmp++)
	{
		string tmpQueryResults = outputDir + idVec[tmp] + ".query";	
		// step 1: query all the kmers
		string query_cmd_1 = "/home/xli262/seqQueryTree/build/SeqOthello/Query2592";
		string query_cmd_2 = "/scratch/lji226/projects/xal_seq/kmer/info/map";
		string query_cmd_3 = kmerFileVec[tmp];
		string query_cmd_4 = "20 N Y N >";
		string query_cmd = query_cmd_1 + " " + query_cmd_2 + " " 
			+ query_cmd_3 + " " + query_cmd_4 + " " + tmpQueryResults;
		cout << "tmpKmerFile: " << kmerFileVec[tmp] << endl;
		system(query_cmd.c_str());
		// step 2: start to examine
		vector<string> tmpSeqVec;
		vector< vector<int> > tmpQueryVecVec;
		loadQueryResults2seqVecKmerVecVec(tmpQueryResults, tmpSeqVec, tmpQueryVecVec);
		for(int tmp2 = 0; tmp2 < tmpQueryVecVec.size(); tmp2++)
		{
			if((tmpQueryVecVec[tmp2])[tmp] != 1)
				err_ofs << tmpSeqVec[tmp2] << "\t" << tmp << endl;
		}
		string rm_cmd = "rm " + tmpQueryResults;
		system(rm_cmd.c_str());
	}
	err_ofs.close();
	return 0;
}