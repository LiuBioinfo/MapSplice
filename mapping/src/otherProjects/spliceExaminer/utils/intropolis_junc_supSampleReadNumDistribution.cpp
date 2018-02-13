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

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIntropolisJuncFile outputFolder supSampleNum_max supReadNum_max" << endl;
		exit(1);
	}
	string inputIntropolisJuncFile = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	string supSampleNum_max_str = argv[3];
	string supReadNum_max_str = argv[4];
	int supSampleNum_max = atoi(supSampleNum_max_str.c_str());
	int supReadNum_max = atoi(supReadNum_max_str.c_str());
	int supSampleNum_max_log2_int = log2((double)supSampleNum_max);
	int supReadNum_max_log2_int = log2((double)supReadNum_max);
	log_ofs << "supSampleNum_max: " << supSampleNum_max << endl;
	log_ofs << "supSampleNum_max_log2_int: " << supSampleNum_max_log2_int << endl;
	log_ofs << "supReadNum_max: " << supReadNum_max << endl;
	log_ofs << "supReadNum_max_log2_int: " << supReadNum_max_log2_int << endl;
	cout << "supSampleNum_max: " << supSampleNum_max << endl;
	cout << "supSampleNum_max_log2_int: " << supSampleNum_max_log2_int << endl;
	cout << "supReadNum_max: " << supReadNum_max << endl;
	cout << "supReadNum_max_log2_int: " << supReadNum_max_log2_int << endl;	

	vector<int> supSample_juncNum_log2_vec;
	for(int tmp = 0; tmp < supSampleNum_max_log2_int; tmp++)
		supSample_juncNum_log2_vec.push_back(0);
	vector<int> supRead_juncNum_log2_vec;
	for(int tmp = 0; tmp < supReadNum_max_log2_int; tmp++)
		supRead_juncNum_log2_vec.push_back(0);
	
	int tmpJuncNum = 0;
	ifstream intropolis_ifs(inputIntropolisJuncFile.c_str());
	while(!intropolis_ifs.eof())
	{
		string tmpStr;
		getline(intropolis_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpJuncNum ++;
		int tmpThousandIndex = tmpJuncNum / 100000;
		if(tmpJuncNum == tmpThousandIndex * 100000)
			cout << "Processed Junc #: " << tmpJuncNum << endl;
		vector<string> tmpJuncFieldVec;
		int tmpStartLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			if(tmpTabLoc != string::npos)
			{	
				string tmpFieldStr = tmpStr.substr(tmpStartLoc, tmpTabLoc-tmpStartLoc);
				tmpJuncFieldVec.push_back(tmpFieldStr);
				tmpStartLoc = tmpTabLoc + 1;
			}
			else
			{	
				string tmpLastFieldStr = tmpStr.substr(tmpStartLoc);
				tmpJuncFieldVec.push_back(tmpLastFieldStr);
				break;
			}
		}
		string tmpJunc_sampleNumStr = tmpJuncFieldVec[5];
		string tmpJunc_readNumStr = tmpJuncFieldVec[6];
		int tmpJunc_sampleNum = atoi(tmpJunc_sampleNumStr.c_str());
		int tmpJunc_readNum = atoi(tmpJunc_readNumStr.c_str());
		int tmpJunc_sampleNum_log2 = log2((double)tmpJunc_sampleNum);
		int tmpJunc_readNum_log2 = log2((double)tmpJunc_readNum);
		supSample_juncNum_log2_vec[tmpJunc_sampleNum_log2] ++;
		supRead_juncNum_log2_vec[tmpJunc_readNum_log2] ++;
	}
	int totalJuncNum = tmpJuncNum;
	string output_supSampleNum_distribution_file = outputFolderStr + "/supSampleNum_distribution.txt";
	string output_supReadNum_distribution_file = outputFolderStr + "/supReadNum_distribution.txt";
	ofstream supSampleNum_distribution_ofs(output_supSampleNum_distribution_file.c_str());
	ofstream supReadNum_distribution_ofs(output_supReadNum_distribution_file.c_str());
	int tmpSupSample_juncNum_cumulative = 0;
	double tmpSupSample_juncPerc_cumulative = 0.0;
	for(int tmp = 0; tmp < supSampleNum_max_log2_int; tmp++)
	{
		int tmpSupSample_min = pow(2, tmp);
		int tmpSupSample_max = pow(2, tmp + 1) - 1;
		tmpSupSample_juncNum_cumulative += supSample_juncNum_log2_vec[tmp];
		double tmpSupSample_juncPerc = ((double)supSample_juncNum_log2_vec[tmp]/(double)totalJuncNum) * 100;
		tmpSupSample_juncPerc_cumulative += tmpSupSample_juncPerc;
		supSampleNum_distribution_ofs << tmp + 1 << "\t" << tmpSupSample_min << "~" 
			<< tmpSupSample_max << "\t" << supSample_juncNum_log2_vec[tmp] << "\t" << tmpSupSample_juncPerc << "%\t"
			<< tmpSupSample_juncNum_cumulative << "\t" << tmpSupSample_juncPerc_cumulative << "%" << endl;
	}
	int tmpSupRead_juncNum_cumulative = 0;
	double tmpSupRead_juncPerc_cumulative = 0.0;
	for(int tmp = 0; tmp < supReadNum_max_log2_int; tmp++)
	{
		int tmpSupRead_min = pow(2, tmp);
		int tmpSupRead_max = pow(2, tmp + 1) - 1;
		tmpSupRead_juncNum_cumulative += supRead_juncNum_log2_vec[tmp];
		double tmpSupRead_juncPerc = ((double)supRead_juncNum_log2_vec[tmp]/(double)totalJuncNum) * 100;
		tmpSupRead_juncPerc_cumulative += tmpSupRead_juncPerc;
		supReadNum_distribution_ofs << tmp + 1 << "\t" << tmpSupRead_min << "~" 
			<< tmpSupRead_max << "\t" << supRead_juncNum_log2_vec[tmp] << "\t" << tmpSupRead_juncPerc << "%\t"
			<< tmpSupRead_juncNum_cumulative << "\t" << tmpSupRead_juncPerc_cumulative << "%" << endl;
	}
	supSampleNum_distribution_ofs.close();
	supReadNum_distribution_ofs.close();
	log_ofs << endl << "totalJuncNum: " << totalJuncNum << endl;
	cout << endl << "totalJuncNum: " << totalJuncNum << endl;
	intropolis_ifs.close();
	log_ofs.close();
	return 0;
}