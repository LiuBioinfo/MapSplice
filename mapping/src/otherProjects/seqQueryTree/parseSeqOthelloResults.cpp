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

void get_idList_valid_and_invalid(
	string& input_idList_file_total, string& input_idList_file_valid, 
	vector<string>& idVec_valid, vector<string>& idVec_invalid)
{
	ifstream valid_ifs(input_idList_file_valid.c_str());
	while(!valid_ifs.eof())
	{
		string tmpStr;
		getline(valid_ifs, tmpStr);
		if(tmpStr == "")
			break;
		idVec_valid.push_back(tmpStr);
	}
	valid_ifs.close();
	ifstream total_ifs(input_idList_file_total.c_str());
	while(!total_ifs.eof())
	{
		string tmpStr;
		getline(total_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//idList_valid.push_back(idList_valid);
		bool invalid_bool = true;
		for(int tmp = 0; tmp < idVec_valid.size(); tmp++)
		{
			if(idVec_valid[tmp] == tmpStr)
			{
				invalid_bool = false;
				break;
			}
		}
		if(invalid_bool)
			idVec_invalid.push_back(tmpStr);
	}
	total_ifs.close();	
}

void loadQueryResults2seqVecKmerVecVec(
	string& input_seqOthelloResults_file, vector<string>& seqVec, 
	vector<int>& seqLengthVec, vector< vector<int> >& queryVecVec, 
	int resultsForEachSeq_line_num_total, int resultsForEachSeq_line_num_valid,
	int invalid_id_num)
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
		int tmpSeqLength = tmpStr.length();
		seqVec.push_back(tmpStr);
		seqLengthVec.push_back(tmpSeqLength);
		vector<int> tmpKmerCountVec;
		for(int tmp = 0; tmp < resultsForEachSeq_line_num_valid - 1; tmp++)
		{
			string tmpStr;
			getline(seqOtehllo_ifs, tmpStr);
			int tmpSpaceLoc = tmpStr.find(" ");
			string tmpKmerCountStr = tmpStr.substr(tmpSpaceLoc + 1);
			int tmpKmerCount = atoi(tmpKmerCountStr.c_str());
			//cout << "tmpKmerCount: " << tmpKmerCount << endl;
			tmpKmerCountVec.push_back(tmpKmerCount);
		}
		for(int tmp = 0; tmp < resultsForEachSeq_line_num_total - resultsForEachSeq_line_num_valid; tmp++)
		{
			string tmpStr;
			getline(seqOtehllo_ifs, tmpStr);
		}
		for(int tmp = 0; tmp < invalid_id_num; tmp++)
			tmpKmerCountVec.push_back(0);
		queryVecVec.push_back(tmpKmerCountVec);
	}
	seqOtehllo_ifs.close();
}

void print_queryResultsVec(
	string& output_query_occurrence_map_tmpThres, string& output_present_sample_num_tmpThres,
	double thres_tmp, vector<int>& seqLengthVec, vector< vector<int> >& queryVecVec)
{
	ofstream query_occurrence_map_ofs_tmpThres(output_query_occurrence_map_tmpThres.c_str());
	ofstream present_sample_num_ofs_tmpThres(output_present_sample_num_tmpThres.c_str());
	for(int tmp = 0; tmp < queryVecVec.size(); tmp++) // each query seq
	{
		int present_sample_num_tmpThres = 0;
		for(int tmp2 = 0; tmp2 < queryVecVec[tmp].size(); tmp2++)
		{
			//cout << "tmpQueryNum: " << (queryVecVec[tmp])[tmp2] << endl;
			//cout << "tmpSeqLength: " << (seqLengthVec[tmp]-19) << endl;
			double tmp_presentKmer_ratio = ((double)((queryVecVec[tmp])[tmp2]))/((double)(seqLengthVec[tmp]-19));
			//cout << "tmp_presentKmer_ratio: " << tmp_presentKmer_ratio << endl;
			if(tmp_presentKmer_ratio >= thres_tmp)
			{	
				query_occurrence_map_ofs_tmpThres << "1";
				present_sample_num_tmpThres ++;
			}
			else
				query_occurrence_map_ofs_tmpThres << "0";
		}
		query_occurrence_map_ofs_tmpThres << endl;
		present_sample_num_ofs_tmpThres << present_sample_num_tmpThres << endl;
	}
	query_occurrence_map_ofs_tmpThres.close();
	present_sample_num_ofs_tmpThres.close();		
}

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_seqOthelloResults_file" << endl;
		cout << "#2 input_idList_file_total" << endl;
		// seq + valid sample #
		cout << "#3 input_idList_file_valid (long enough reads,";
		cout << "have kmers, consistent order with that in results file)" << endl;
		// including the seq line and all the results lines
		cout << "#4 resultsForEachSeq_line_num_total" << endl; 
		// including the seq line and all the "valid" results lines
		cout << "#5 resultsForEachSeq_line_num_valid" << endl;
		cout << "#6 input_thres_list_file" << endl;
		cout << "#7 output_formatted_results_dir" << endl;
		exit(1);
	}
	string input_seqOthelloResults_file = argv[1];
	string input_idList_file_total = argv[2];
	string input_idList_file_valid = argv[3];
	string resultsForEachSeq_line_num_total_str = argv[4];
	string resultsForEachSeq_line_num_valid_str = argv[5];
	string input_thres_list_file = argv[6];
	string output_formatted_results_dir = argv[7];

	//////////////////////////////////////
	//////////////////////////////////////
	vector<double> thresVec;
	vector<string> thresVec_str;
	vector<string> idList_valid;
	vector<string> idList_invalid;
	vector<string> seqVec;
	vector<int> seqLengthVec;
	vector< vector<int> > queryVecVec;
	//////////////////////////////////////
	//////////////////////////////////////

	cout << "start to make output dir" << endl;
	output_formatted_results_dir += "/";
	string cmd_mkdir_output_formatted_results_dir = "mkdir " + output_formatted_results_dir;
	system(cmd_mkdir_output_formatted_results_dir.c_str());

	cout << "start to load thres list" << endl;
	ifstream thres_ifs(input_thres_list_file.c_str());
	while(!thres_ifs.eof())
	{
		string tmpStr;
		getline(thres_ifs, tmpStr);
		if(tmpStr == "")
			break;
		double tmpThres = atof(tmpStr.c_str());
		thresVec_str.push_back(tmpStr);
		thresVec.push_back(tmpThres);
	}
	thres_ifs.close();

	cout << "start to print the id list" << endl;
	string output_formatted_results_id_list = output_formatted_results_dir + "id_list.txt";
	ofstream id_list_ofs(output_formatted_results_id_list.c_str());		
	int resultsForEachSeq_line_num_total = atoi(resultsForEachSeq_line_num_total_str.c_str());
	int resultsForEachSeq_line_num_valid = atoi(resultsForEachSeq_line_num_valid_str.c_str());
	if(resultsForEachSeq_line_num_valid > resultsForEachSeq_line_num_total)
	{
		cout << "Error! resultsForEachSeq_line_num_valid > resultsForEachSeq_line_num_total!" << endl;
		exit(1);
	}
	get_idList_valid_and_invalid(input_idList_file_total, input_idList_file_valid, idList_valid, idList_invalid);
	for(int tmp = 0; tmp < idList_valid.size(); tmp++)
		id_list_ofs << idList_valid[tmp] << "\tvalid" << endl;
	for(int tmp = 0; tmp < idList_invalid.size(); tmp++)
		id_list_ofs << idList_invalid[tmp] << "\tinvalid" << endl;	
	id_list_ofs.close();

	cout << "start to load seq query results to seqVec and KmerVecVec" << endl;
	loadQueryResults2seqVecKmerVecVec(input_seqOthelloResults_file, 
		seqVec, seqLengthVec, queryVecVec, resultsForEachSeq_line_num_total, 
		resultsForEachSeq_line_num_valid, idList_invalid.size());

	cout << "start to print seq" << endl;
	string output_formatted_results_query_seq = output_formatted_results_dir + "query_seq.txt";
	ofstream query_seq_ofs(output_formatted_results_query_seq.c_str());
	for(int tmp = 0; tmp < seqVec.size(); tmp++)
		query_seq_ofs << seqVec[tmp] << endl;
	query_seq_ofs.close();

	cout << "start to print each thres results" << endl;
	for(int tmp = 0; tmp < thresVec.size(); tmp++)
	{	
		cout << "tmpThres: " << thresVec[tmp] << endl;
		string output_present_sample_num_tmpThres = output_formatted_results_dir 
			+ "present_sample_num." + thresVec_str[tmp] + ".txt";
		string output_query_occurrence_map_tmpThres = output_formatted_results_dir 
			+ "query_occurrence_map." + thresVec_str[tmp] + ".txt";
		print_queryResultsVec(output_query_occurrence_map_tmpThres, 
			output_present_sample_num_tmpThres, thresVec[tmp], seqLengthVec, queryVecVec);
	}
	cout << "All jobs done!" << endl;
	return 0;
}