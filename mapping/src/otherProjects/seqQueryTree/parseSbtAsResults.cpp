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

string getIdFromIdListFileLineStr(string& tmpStr)
{
	return tmpStr;
}

void get_querySeq_presentSampleNum_visitedNodeNum(string& str_1, string& str_2,
	string& query_seq, int& present_sample_num, int& visited_node_num)
{
	int space_loc_line1 = str_1.find(" ");
	query_seq = str_1.substr(1, space_loc_line1 - 1);
	string present_sample_num_str = str_1.substr(space_loc_line1 + 1);
	int space_loc_line2_1 = str_2.find(" ");
	int space_loc_line2_2 = str_2.find(" ", space_loc_line2_1 + 1);
	string visited_node_num_str = str_2.substr(space_loc_line2_1 + 1, space_loc_line2_2 - space_loc_line2_1 - 1);
	present_sample_num = atoi(present_sample_num_str.c_str());
	visited_node_num = atoi(visited_node_num_str.c_str());
}

int convertId2index(map<string, int>& id2indexMap, string& tmp_id)
{
	map<string, int>::iterator tmpIter = id2indexMap.find(tmp_id);
	if(tmpIter != id2indexMap.end())
		return (tmpIter->second);
	else
	{
		cout << "tmp_id not found in map: " << tmp_id << endl;
		exit(1);
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_SBT_AS_results" << endl;
		cout << "#2 input_file_id_list" << endl;
		cout << "#3 output_formatted_results_dir" << endl;
		//cout << "#4 input_uniform_srr_id_list_file" << endl;
		exit(1);
	}

	string input_SBT_AS_results = argv[1];
	string input_file_id_list = argv[2];
	string output_formatted_results_dir = argv[3]; output_formatted_results_dir += "/";
	string input_uniform_srr_id_list_file = argv[4];
	//string mkdir_output = "mkdir " + output_formatted_results_dir;
	system(mkdir_output.c_str());

	string output_formatted_results_query_seq = output_formatted_results_dir + "query_seq.txt";
	string output_formatted_results_present_sample_num = output_formatted_results_dir + "present_sample_num.txt";
	string output_formatted_results_visited_node_num = output_formatted_results_dir + "visited_node_num.txt";
	string output_formatted_results_query_occurrence_map = output_formatted_results_dir + "query_occurrence_map.txt";

	ofstream query_seq_ofs(output_formatted_results_query_seq.c_str());
	ofstream present_sample_num_ofs(output_formatted_results_present_sample_num.c_str());
	ofstream visited_node_num_ofs(output_formatted_results_visited_node_num.c_str());
	ofstream query_occurrence_map_ofs(output_formatted_results_query_occurrence_map.c_str());

	cout << "start to load id_list ..." << endl;
	vector<string> idVec;
	ifstream id_list_ifs(input_file_id_list.c_str());
	while(!id_list_ifs.eof())
	{
		string tmpStr;
		getline(id_list_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpId = getIdFromIdListFileLineStr(tmpStr);
		idVec.push_back(tmpId);
	}
	id_list_ifs.close();

	cout << "start to build id2indexMap ..." << endl;
	map<string, int> id2indexMap;
	for(int tmp = 0; tmp < idVec.size(); tmp++)
		id2indexMap.insert(pair<string,int>(idVec[tmp], tmp));

	cout << "start to check whether id2indexMap is a 1 to 1 map" << endl;
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		map<string, int>::iterator tmpIter = id2indexMap.find(tmpId);
		if(tmpIter != id2indexMap.end())
		{
			int tmpIndex = tmpIter->second;
			if(tmpIndex != tmp)
			{
				cout << "tmpIndex not consistent" << endl;
				exit(1);
			}
		}
		else
		{
			cout << "Not able to find tmpId" << endl;
			exit(1);
		}
	}
	cout << "id2indexMap building successed" << endl;

	cout << "start to parse SBT-AS results" << endl;
	vector<string> querySeqVec;
	vector<int> presentSampleNumVec;
	vector<int> visitedNodeNumVec;
	vector< vector<bool> > occurrenceMapVecVec;
	ifstream SBT_AS_ifs(input_SBT_AS_results.c_str());
	int tmpQueryNum = 0;
	while(!SBT_AS_ifs.eof())
	{	
		string tmpStr_1, tmpStr_2;
		getline(SBT_AS_ifs, tmpStr_1);
		if(tmpStr_1 == "")
			break;
		getline(SBT_AS_ifs, tmpStr_2);
		if(tmpStr_2 == "")
			break;		
		tmpQueryNum ++;
		if((tmpQueryNum / 10000) * 10000 == tmpQueryNum)
			cout << "tmpQueryNum: " << tmpQueryNum << endl;
		vector<bool> tmpOccurrenceMapVec;
		for(int tmp = 0; tmp < idVec.size(); tmp++)
			tmpOccurrenceMapVec.push_back(false);
		string tmp_query_seq;
		int tmp_present_sample_num, tmp_visited_node_num;
		get_querySeq_presentSampleNum_visitedNodeNum(tmpStr_1, tmpStr_2,
			tmp_query_seq, tmp_present_sample_num, tmp_visited_node_num);
		for(int tmpLine = 0; tmpLine < tmp_present_sample_num; tmpLine++)
		{
			string tmpStr;
			getline(SBT_AS_ifs, tmpStr);
			if(tmpStr == "")
				break;
			//cout << "tmpStr: " << tmpStr << endl;
			string tmp_id = tmpStr;
			int tmp_index = convertId2index(id2indexMap,tmp_id);
			tmpOccurrenceMapVec[tmp_index] = true;
		}
		querySeqVec.push_back(tmp_query_seq);
		presentSampleNumVec.push_back(tmp_present_sample_num);
		visitedNodeNumVec.push_back(tmp_visited_node_num);
		occurrenceMapVecVec.push_back(tmpOccurrenceMapVec);
	}
	SBT_AS_ifs.close();
	cout << "start to print query_seq ..." << endl;
	for(int tmp = 0; tmp < querySeqVec.size(); tmp++)
		query_seq_ofs << querySeqVec[tmp] << endl;
	query_seq_ofs.close();
	cout << "start to print present_sample_num ..." << endl;
	for(int tmp = 0; tmp < presentSampleNumVec.size(); tmp++)
		present_sample_num_ofs << presentSampleNumVec[tmp] << endl;
	present_sample_num_ofs.close();
	cout << "start to print visited_node_num ..." << endl;
	for(int tmp = 0; tmp < visitedNodeNumVec.size(); tmp++)
		visited_node_num_ofs << visitedNodeNumVec[tmp] << endl;
	visited_node_num_ofs.close();
	cout << "start to print query_occurrence_map ..." << endl; 
	for(int tmp = 0; tmp < occurrenceMapVecVec.size(); tmp++)
	{
		for(int tmp2 = 0; tmp2 < occurrenceMapVecVec[tmp].size(); tmp2++)
		{
			if((occurrenceMapVecVec[tmp])[tmp2])
				query_occurrence_map_ofs << "1";
			else
				query_occurrence_map_ofs << "0";
		}
		query_occurrence_map_ofs << endl;
	}
	query_occurrence_map_ofs.close();
	return 0;
}