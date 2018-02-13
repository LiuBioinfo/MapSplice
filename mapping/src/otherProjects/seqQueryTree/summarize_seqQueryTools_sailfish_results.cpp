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

void parseStr2FieldVec(string& tmpStr, vector<string>& fieldVec)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			fieldVec.push_back(tmpStr.substr(startLoc));
			return;
		}
		else
		{
			string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
			fieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
	}
	return;
}

void parseMergedSailfishResults(string& inputMergedSailfishResults, 
	vector<string>& sampleIdVec, vector<string>& transcriptIdVec,
	vector< vector<double> >& expVecVec)
{
	ifstream sailfish_ifs(inputMergedSailfishResults.c_str());
	string tmpStr_1st;
	getline(sailfish_ifs, tmpStr_1st);
	vector<string> fieldVec_1stLine;
	parseStr2FieldVec(tmpStr_1st, fieldVec_1stLine);
	for(int tmp = 1; tmp < fieldVec_1stLine.size(); tmp++)
	{
		sampleIdVec.push_back(fieldVec_1stLine[tmp]);
		vector<double> tmpExpVec;
		expVecVec.push_back(tmpExpVec);
	}
	int tmpLineNum = 0;
	while(!sailfish_ifs.eof())
	{
		string tmpStr;
		tmpLineNum ++;
		getline(sailfish_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if((tmpLineNum/10000) * 10000 == tmpLineNum)
			cout << tmpLineNum << " transcripts processed" << endl;
		vector<string> fieldVec_tmp;
		parseStr2FieldVec(tmpStr, fieldVec_tmp);
		transcriptIdVec.push_back(fieldVec_tmp[0]);
		for(int tmp = 0; tmp < sampleIdVec.size(); tmp++)
			expVecVec[tmp].push_back(atof(fieldVec_tmp[tmp+1].c_str()));
	}

	sailfish_ifs.close();
}

double queryExp_transcriptId_sampleId(vector<string>& sampleIdVec, vector<string>& transcriptIdVec,
	vector< vector<double> >& expVecVec, string& tmp_sampleId, string& tmp_transcriptId)
	// return value < 0 if sample_id or transcript_id invalid
{
	int index_sample = -1;
	int index_transcript = -1;
	for(int tmp = 0; tmp < sampleIdVec.size(); tmp++)
	{
		if(sampleIdVec[tmp] == tmp_sampleId)
		{
			index_sample = tmp;
			break;
		}
	}
	for(int tmp = 0; tmp < transcriptIdVec.size(); tmp++)
	{
		if(transcriptIdVec[tmp] == tmp_transcriptId)
		{
			index_transcript = tmp;
			break;
		}
	}
	if((index_sample >= 0)&&(index_transcript >= 0))
		return (expVecVec[index_sample])[index_transcript];
	else
	{
		cout << "Error!" << endl;
		if(index_sample < 0)
			cout << "sample_id not found: " << tmp_sampleId << endl;
		if(index_transcript < 0)
			cout << "transcript_id not found: " << tmp_transcriptId << endl;
		return -1.0;
	}
}

void get_tpm_maxAsAbsent_minAsPresent(string& tpm_range_file, 
	double& tpm_max_as_absent, double& tpm_min_as_present)
{
	ifstream tpm_range_ifs(tpm_range_file.c_str());
	string tpmStr_1, tpmStr_2;
	getline(tpm_range_ifs, tpmStr_1);
	getline(tpm_range_ifs, tpmStr_2);
	tpm_max_as_absent = atof(tpmStr_1.c_str());
	tpm_min_as_present = atof(tpmStr_2.c_str());
	tpm_range_ifs.close();
}

void get_thresVec(string& thres_theta_list_file, vector<double>& thres_vec)
{
	ifstream theta_ifs(thres_theta_list_file.c_str());
	while(!theta_ifs.eof())
	{
		string tmpStr;
		getline(theta_ifs, tmpStr);
		if(tmpStr == "")
			break;
		thres_vec.push_back(atof(tmpStr.c_str()));
	}
	theta_ifs.close();
}

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

void load01formatResults(string& tmpFile, vector< vector<bool> >& tmpResultsVecVec)
{
	//parseStr2FieldVec(string& tmpStr, vector<string>& fieldVec)
	ifstream tmp_ifs(tmpFile.c_str());
	while(!tmp_ifs.eof())
	{
		string tmpStr;
		getline(tmp_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//vector<string> tmpFieldVec;
		//parseStr2FieldVec(tmpStr, tmpFieldVec);
		vector<bool> tmpResultsVec;
		for(int tmp = 0; tmp < tmpStr.length(); tmp++)
		{
			if(tmpStr[tmp] == '0')
				tmpResultsVec.push_back(false);
			else if(tmpStr[tmp] == '1')
				tmpResultsVec.push_back(true);
			else
			{
				cout << "invalid tmpField in resutls file: " << tmpStr[tmp] << endl;
				exit(1); 
			}
		}
		tmpResultsVecVec.push_back(tmpResultsVec);
	}
	tmp_ifs.close();
}

void get_occurrence_situation_from_occurrence_map(string& results_file,
	int minSeqLength, vector<int>& seqLengthVec, vector< vector<double> >& expVecVec, double& tpm_max_as_absent, double& tpm_min_as_present,
	int& tmp_sailfish_exist_seqQueryTool_exist_num, int& tmp_sailfish_nonexist_seqQueryTool_nonexist_num,
	int& tmp_sailfish_exist_seqQueryTool_nonexist_num, int& tmp_sailfish_nonexist_seqQueryTool_exist_num)
{
	vector< vector<bool> > tmpResultsVecVec;
	load01formatResults(results_file, tmpResultsVecVec);
	for(int tmp1 = 0; tmp1 < expVecVec.size(); tmp1 ++)
	{
		int tmpLength = seqLengthVec[tmp1];
		if(tmpLength < minSeqLength)
			continue;
		for(int tmp2 = 0; tmp2 < expVecVec[tmp1].size(); tmp2 ++)
		{	
			bool tmp_seqQueryTool_exist_bool = (tmpResultsVecVec[tmp1])[tmp2];
			double tmpExp = (expVecVec[tmp1])[tmp2];
			if(tmp_seqQueryTool_exist_bool) // present
			{
				if(tmpExp <= tpm_max_as_absent)
					tmp_sailfish_nonexist_seqQueryTool_exist_num ++;
				else if(tmpExp > tpm_min_as_present)
					tmp_sailfish_exist_seqQueryTool_exist_num ++;
				else
				{}
			}
			else // absent
			{
				if(tmpExp <= tpm_max_as_absent)
					tmp_sailfish_nonexist_seqQueryTool_nonexist_num ++;
				else if(tmpExp > tpm_min_as_present)
					tmp_sailfish_exist_seqQueryTool_nonexist_num ++;
				else
				{}
			}
		}
	}
}

int main(int argc, char** argv)
{
	if((argc != 9)&&(argc != 11)&&(argc != 13)&&(argc != 15))
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputQuerySeqFile" << endl;
		cout << "#2 Sailfish_results" << endl;
		cout << "#3 tpm_max_as_absent" << endl; // if tpm <= x, absent
		cout << "#4 tpm_min_as_present" << endl; // if tpm > x, present
		cout << "#5 thres_theta_list" << endl;
		cout << "#6 minSeqLength" << endl;
		cout << "#7 output_file" << endl;
		cout << "#8 seqOthello_results_list" << endl;
		cout << "(#9 otherTool_name_1" << endl;
		cout << " #10 otherTool_results_1_list)" << endl;
		cout << "(#11 otherTool_name_2" << endl;
		cout << " #12 otherTool_results_2_list)" << endl;
		cout << "(#13 otherTool_name_3" << endl;
		cout << " #14 otherTool_results_3_list)" << endl;
		exit(1);
	}
	string inputQuerySeqFile = argv[1];
	string Sailfish_results_file = argv[2];
	string tpm_max_as_absent_str = argv[3];
	string tpm_min_as_present_str = argv[4];
	string thres_theta_list_file = argv[5];
	string minSeqLengthStr = argv[6];
	string output_file = argv[7];
	string seqOthello_results_list = argv[8];

	int minSeqLength = atoi(minSeqLengthStr.c_str());
	double tpm_max_as_absent = atof(tpm_max_as_absent_str.c_str());
	double tpm_min_as_present = atof(tpm_min_as_present_str.c_str());

	cout << "start to get seqLengthVec " << endl;
	ifstream querySeq_ifs(inputQuerySeqFile.c_str());
	vector<int> seqLengthVec;
	while(!querySeq_ifs.eof())
	{
		string tmpStr;
		getline(querySeq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpLength = tmpStr.length();
		seqLengthVec.push_back(tmpLength);
	}
	querySeq_ifs.close();

	cout << "Start to parse merged sailfish results ..." << endl;
	vector<string> sampleIdVec; 
	vector<string> transcriptIdVec;
	vector< vector<double> > expVecVec;
	parseMergedSailfishResults(Sailfish_results_file, 
		sampleIdVec, transcriptIdVec, expVecVec);

	cout << "start to generate results file vec" << endl;
	vector<string> seqOtehllo_results_vec;
	get_resultsFileVec(seqOthello_results_list, seqOtehllo_results_vec);
	int otherToolNum = 0;
	string otherTool_name_1, otherTool_results_list_1,
		   otherTool_name_2, otherTool_results_list_2,
		   otherTool_name_3, otherTool_results_list_3;
	vector<string> otherTool_results_vec_1, otherTool_results_vec_2, otherTool_results_vec_3;
	
	if(argc > 9)
	{
		otherToolNum ++;
		otherTool_name_1 = argv[9];
		otherTool_results_list_1 = argv[10];
		get_resultsFileVec(otherTool_results_list_1, otherTool_results_vec_1);
		if(argc > 11)
		{
			otherToolNum ++;
			otherTool_name_2 = argv[11];
			otherTool_results_list_2 = argv[12];
			get_resultsFileVec(otherTool_results_list_2, otherTool_results_vec_2);
			if(argc > 13)
			{
				otherToolNum ++;
				otherTool_name_3 = argv[13];
				otherTool_results_list_3 = argv[14];
				get_resultsFileVec(otherTool_results_list_3, otherTool_results_vec_3);
			}
		}
	}
	cout << "start parse thres file" << endl;
	vector<double> thres_vec;
	get_thresVec(thres_theta_list_file, thres_vec);
	cout << "start to generate # for each occurrence case" << endl;
	vector< vector<int> > sailfish_exist_seqQueryTool_exist_num_vecVec;
	vector< vector<int> > sailfish_nonexist_seqQueryTool_nonexist_num_vecVec;
	vector< vector<int> > sailfish_exist_seqQueryTool_nonexist_num_vecVec;
	vector< vector<int> > sailfish_nonexist_seqQueryTool_exist_num_vecVec;
	for(int tmpThres = 0; tmpThres < thres_vec.size(); tmpThres ++)
	{
		cout << "tmpThres: " << thres_vec[tmpThres] << endl;
		vector<int> tmp_sailfish_exist_seqQueryTool_exist_num_vec;
		vector<int> tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec;
		vector<int> tmp_sailfish_exist_seqQueryTool_nonexist_num_vec;
		vector<int> tmp_sailfish_nonexist_seqQueryTool_exist_num_vec;
		cout << "start to get occurrence case num for seqOtehllo results" << endl;
		int tmp_sailfish_exist_seqOthello_exist_num = 0;
		int tmp_sailfish_nonexist_seqOthello_nonexist_num = 0;
		int tmp_sailfish_exist_seqOthello_nonexist_num = 0;
		int tmp_sailfish_nonexist_seqOthello_exist_num = 0;
		get_occurrence_situation_from_occurrence_map(seqOtehllo_results_vec[tmpThres], minSeqLength,
			seqLengthVec, expVecVec, tpm_max_as_absent, tpm_min_as_present,
			tmp_sailfish_exist_seqOthello_exist_num, tmp_sailfish_nonexist_seqOthello_nonexist_num,
			tmp_sailfish_exist_seqOthello_nonexist_num, tmp_sailfish_nonexist_seqOthello_exist_num);
		tmp_sailfish_exist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_exist_seqOthello_exist_num);
		tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_nonexist_seqOthello_nonexist_num);
		tmp_sailfish_exist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_exist_seqOthello_nonexist_num);
		tmp_sailfish_nonexist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_nonexist_seqOthello_exist_num);
		if(otherToolNum > 0)
		{
			cout << "start to get occurrence case num for othertool_1 results" << endl;
			int tmp_sailfish_exist_otherTool_1_exist_num = 0;
			int tmp_sailfish_nonexist_otherTool_1_nonexist_num = 0;
			int tmp_sailfish_exist_otherTool_1_nonexist_num = 0;
			int tmp_sailfish_nonexist_otherTool_1_exist_num = 0;
			get_occurrence_situation_from_occurrence_map(otherTool_results_vec_1[tmpThres], minSeqLength,
				seqLengthVec, expVecVec, tpm_max_as_absent, tpm_min_as_present,
				tmp_sailfish_exist_otherTool_1_exist_num, tmp_sailfish_nonexist_otherTool_1_nonexist_num,
				tmp_sailfish_exist_otherTool_1_nonexist_num, tmp_sailfish_nonexist_otherTool_1_exist_num);
			tmp_sailfish_exist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_exist_otherTool_1_exist_num);
			tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_1_nonexist_num);
			tmp_sailfish_exist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_exist_otherTool_1_nonexist_num);
			tmp_sailfish_nonexist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_1_exist_num);
			if(otherToolNum > 1)
			{
				cout << "start to get occurrence case num for othertool_2 results" << endl;
				int tmp_sailfish_exist_otherTool_2_exist_num = 0;
				int tmp_sailfish_nonexist_otherTool_2_nonexist_num = 0;
				int tmp_sailfish_exist_otherTool_2_nonexist_num = 0;
				int tmp_sailfish_nonexist_otherTool_2_exist_num = 0;
				get_occurrence_situation_from_occurrence_map(otherTool_results_vec_2[tmpThres], minSeqLength,
					seqLengthVec, expVecVec, tpm_max_as_absent, tpm_min_as_present,
					tmp_sailfish_exist_otherTool_2_exist_num, tmp_sailfish_nonexist_otherTool_2_nonexist_num,
					tmp_sailfish_exist_otherTool_2_nonexist_num, tmp_sailfish_nonexist_otherTool_2_exist_num);
				tmp_sailfish_exist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_exist_otherTool_2_exist_num);
				tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_2_nonexist_num);
				tmp_sailfish_exist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_exist_otherTool_2_nonexist_num);
				tmp_sailfish_nonexist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_2_exist_num);
				if(otherToolNum > 2)
				{
					cout << "start to get occurrence case num for othertool_3 results" << endl;
					int tmp_sailfish_exist_otherTool_3_exist_num = 0;
					int tmp_sailfish_nonexist_otherTool_3_nonexist_num = 0;
					int tmp_sailfish_exist_otherTool_3_nonexist_num = 0;
					int tmp_sailfish_nonexist_otherTool_3_exist_num = 0;
					get_occurrence_situation_from_occurrence_map(otherTool_results_vec_3[tmpThres], minSeqLength,
						seqLengthVec, expVecVec, tpm_max_as_absent, tpm_min_as_present,
						tmp_sailfish_exist_otherTool_3_exist_num, tmp_sailfish_nonexist_otherTool_3_nonexist_num,
						tmp_sailfish_exist_otherTool_3_nonexist_num, tmp_sailfish_nonexist_otherTool_3_exist_num);
					tmp_sailfish_exist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_exist_otherTool_3_exist_num);
					tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_3_nonexist_num);
					tmp_sailfish_exist_seqQueryTool_nonexist_num_vec.push_back(tmp_sailfish_exist_otherTool_3_nonexist_num);
					tmp_sailfish_nonexist_seqQueryTool_exist_num_vec.push_back(tmp_sailfish_nonexist_otherTool_3_exist_num);
				}
			}
		}
		sailfish_exist_seqQueryTool_exist_num_vecVec.push_back(tmp_sailfish_exist_seqQueryTool_exist_num_vec);
		sailfish_nonexist_seqQueryTool_nonexist_num_vecVec.push_back(tmp_sailfish_nonexist_seqQueryTool_nonexist_num_vec);
		sailfish_exist_seqQueryTool_nonexist_num_vecVec.push_back(tmp_sailfish_exist_seqQueryTool_nonexist_num_vec);
		sailfish_nonexist_seqQueryTool_exist_num_vecVec.push_back(tmp_sailfish_nonexist_seqQueryTool_exist_num_vec);
	}

	cout << "start to get tpr and fpr for each tool under each threshold" << endl;
	ofstream roc_ofs(output_file.c_str());
	for(int tmp1 = 0; tmp1 < thres_vec.size(); tmp1 ++)
	{
		double tmpThres = thres_vec[tmp1];		
		roc_ofs << endl << "tmpThres:\t" << tmpThres << endl;
		int tmp_tp, tmp_fp, tmp_tn, tmp_fn;
		double tmp_tpr, tmp_fpr;

		roc_ofs << "\tSeqOthello:" << endl;
		tmp_tp = (sailfish_exist_seqQueryTool_exist_num_vecVec[tmp1])[0];
		tmp_fp = (sailfish_nonexist_seqQueryTool_exist_num_vecVec[tmp1])[0];
		tmp_tn = (sailfish_nonexist_seqQueryTool_nonexist_num_vecVec[tmp1])[0];
		tmp_fn = (sailfish_exist_seqQueryTool_nonexist_num_vecVec[tmp1])[0];
		tmp_tpr = (double)tmp_tp/(double)(tmp_tp + tmp_fn); 
		tmp_fpr = (double)tmp_fp/(double)(tmp_tn + tmp_fp);
		roc_ofs << "\t\tTP:\t" << tmp_tp << endl << "\t\tFP:\t" << tmp_fp << endl
			<< "\t\tTN:\t" << tmp_tn << endl << "\t\tFN:\t" << tmp_fn << endl
			<< "\t\tTPR:\t" << tmp_tpr << endl << "\t\tFPR:\t" << tmp_tpr << endl;
		if(otherToolNum > 0)
		{
			roc_ofs << "\t" << otherTool_name_1 << ":" << endl;
			tmp_tp = (sailfish_exist_seqQueryTool_exist_num_vecVec[tmp1])[1];
			tmp_fp = (sailfish_nonexist_seqQueryTool_exist_num_vecVec[tmp1])[1];
			tmp_tn = (sailfish_nonexist_seqQueryTool_nonexist_num_vecVec[tmp1])[1];
			tmp_fn = (sailfish_exist_seqQueryTool_nonexist_num_vecVec[tmp1])[1];
			tmp_tpr = (double)tmp_tp/(double)(tmp_tp + tmp_fn); 
			tmp_fpr = (double)tmp_fp/(double)(tmp_tn + tmp_fp);
			roc_ofs << "\t\tTP:\t" << tmp_tp << endl << "\t\tFP:\t" << tmp_fp << endl
				<< "\t\tTN:\t" << tmp_tn << endl << "\t\tFN:\t" << tmp_fn << endl
				<< "\t\tTPR:\t" << tmp_tpr << endl << "\t\tFPR:\t" << tmp_tpr << endl;

			if(otherToolNum > 1)
			{
				roc_ofs << "\t" << otherTool_name_2 << ":" << endl;
				tmp_tp = (sailfish_exist_seqQueryTool_exist_num_vecVec[tmp1])[2];
				tmp_fp = (sailfish_nonexist_seqQueryTool_exist_num_vecVec[tmp1])[2];
				tmp_tn = (sailfish_nonexist_seqQueryTool_nonexist_num_vecVec[tmp1])[2];
				tmp_fn = (sailfish_exist_seqQueryTool_nonexist_num_vecVec[tmp1])[2];
				tmp_tpr = (double)tmp_tp/(double)(tmp_tp + tmp_fn); 
				tmp_fpr = (double)tmp_fp/(double)(tmp_tn + tmp_fp);
				roc_ofs << "\t\tTP:\t" << tmp_tp << endl << "\t\tFP:\t" << tmp_fp << endl
					<< "\t\tTN:\t" << tmp_tn << endl << "\t\tFN:\t" << tmp_fn << endl
					<< "\t\tTPR:\t" << tmp_tpr << endl << "\t\tFPR:\t" << tmp_tpr << endl;				
				if(otherToolNum > 2)
				{
					roc_ofs << "\t" << otherTool_name_3 << ":" << endl;
					tmp_tp = (sailfish_exist_seqQueryTool_exist_num_vecVec[tmp1])[3];
					tmp_fp = (sailfish_nonexist_seqQueryTool_exist_num_vecVec[tmp1])[3];
					tmp_tn = (sailfish_nonexist_seqQueryTool_nonexist_num_vecVec[tmp1])[3];
					tmp_fn = (sailfish_exist_seqQueryTool_nonexist_num_vecVec[tmp1])[3];
					tmp_tpr = (double)tmp_tp/(double)(tmp_tp + tmp_fn); 
					tmp_fpr = (double)tmp_fp/(double)(tmp_tn + tmp_fp);
					roc_ofs << "\t\tTP:\t" << tmp_tp << endl << "\t\tFP:\t" << tmp_fp << endl
						<< "\t\tTN:\t" << tmp_tn << endl << "\t\tFN:\t" << tmp_fn << endl
						<< "\t\tTPR:\t" << tmp_tpr << endl << "\t\tFPR:\t" << tmp_tpr << endl;
				}
			}
		}
	}
	roc_ofs.close();
	return 0;
}