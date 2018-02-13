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
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 seq_fa_file" << endl;
		cout << "#2 srr_id_list_file" << endl;
		cout << "#3 tool_1_results" << endl;
		cout << "#4 tool_2_results" << endl;
		cout << "#5 outputDir" << endl;
		exit(1);
	}
	string seq_fa_file = argv[1];
	string srr_id_list_file = argv[2];
	string tool_1_results = argv[3];
	string tool_2_results = argv[4];
	string outputDir = argv[5]; outputDir += "/";

	string cmd_mkdir = "mkdir " + outputDir;
	system(cmd_mkdir.c_str());

	cout << "start to load fa file" << endl;
	vector<string> fa_id_vec;
	vector<string> fa_seq_vec;
	ifstream fa_ifs(seq_fa_file.c_str());
	while(!fa_ifs.eof())
	{
		string tmpId;
		getline(fa_ifs, tmpId);
		if(tmpId == "")
			break;
		string tmpSeq;
		getline(fa_ifs, tmpSeq);
		fa_id_vec.push_back(tmpId);
		fa_seq_vec.push_back(tmpSeq);
	}
	fa_ifs.close();

	cout << "start to load srr_id file" << endl;
	vector<string> srr_id_vec;
	ifstream srr_ifs(srr_id_list_file.c_str());
	while(!srr_ifs.eof())
	{
		string tmpId;
		getline(srr_ifs, tmpId);
		if(tmpId == "")
			break;
		srr_id_vec.push_back(tmpId);
	}
	srr_ifs.close();

	cout << "start to load tool_1 file" << endl;
	vector<string> results_1_vec;
	ifstream tool_1_ifs(tool_1_results.c_str());
	while(!tool_1_ifs.eof())
	{
		string tmpStr;
		getline(tool_1_ifs, tmpStr);
		if(tmpStr == "")
			break;
		results_1_vec.push_back(tmpStr);	
	}
	tool_1_ifs.close();

	cout << "start to load tool_2 file" << endl;
	vector<string> results_2_vec;
	ifstream tool_2_ifs(tool_2_results.c_str());
	while(!tool_2_ifs.eof())
	{
		string tmpStr;
		getline(tool_2_ifs, tmpStr);
		if(tmpStr == "")
			break;
		results_2_vec.push_back(tmpStr);	
	}
	tool_2_ifs.close();	

	cout << "start to print diff" << endl;
	string output_seq_fa = outputDir + "seq.fa";
	string output_diff = outputDir + "diff.txt";
	ofstream diff_seq_ofs(output_seq_fa.c_str());
	ofstream diff_id_ofs(output_diff.c_str());
	for(int tmp = 0; tmp < results_1_vec.size(); tmp++)
	{
		if(results_1_vec[tmp] != results_2_vec[tmp])
		{
			diff_seq_ofs << fa_id_vec[tmp] << endl
				<< fa_seq_vec[tmp] << endl;
			diff_id_ofs << fa_id_vec[tmp] << endl;
			for(int tmp2 = 0; tmp2 < srr_id_vec.size(); tmp2++)
			{
				if((results_1_vec[tmp])[tmp2] != (results_2_vec[tmp])[tmp2])
					diff_id_ofs << srr_id_vec[tmp2] << "\t"
						<< (results_1_vec[tmp])[tmp2] << "\t"
						<< (results_2_vec[tmp])[tmp2] << "\t" << endl;
			}
		}
	}
	diff_id_ofs.close();
	diff_seq_ofs.close();
	return 0;
}