// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
// sailfish results format:
//Name	Length	EffectiveLength	TPM	NumReads
//ENST00000335137.3|ENSG00000186092.4|OTTHUMG00000001094.1|OTTHUMT00000003223.1|OR4F5-001|OR4F5|918|CDS:1-918|	918	717.611	0	0
//ENST00000423372.3|ENSG00000237683.5|-|-|AL627309.1-201|AL627309.1|2661|UTR5:1-70|CDS:71-850|UTR3:851-2661|	2661	2460.61	30.7966	4997
//ENST00000426406.1|ENSG00000235249.1|OTTHUMG00000002860.1|OTTHUMT00000007999.1|OR4F29-001|OR4F29|995|UTR5:1-19|CDS:20-958|UTR3:959-995|	995	794.611	0.0190846	1
//ENST00000332831.2|ENSG00000185097.2|OTTHUMG00000002581.1|OTTHUMT00000007334.1|OR4F16-001|OR4F16|995|UTR5:1-19|CDS:20-958|UTR3:959-995|	995	794.611	0.0190846	1
//ENST00000599533.1|ENSG00000269831.1|-|-|AL669831.1-201|AL669831.1|129|CDS:1-129|	129	40.7436	0	0

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

void parseSailfishResultsStr(string& tmpStr, string& tmpId, int& tmpLength, 
	double& tmpEffectiveLength, double& tmpTpm, double& tmpReadNum)
{
	int tabLoc_1 = tmpStr.find("\t");
	int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);

	tmpId = tmpStr.substr(0, tabLoc_1);
	string tmpLengthStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEffectiveLengthStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	string tmpTpmStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	string tmpReadNumStr = tmpStr.substr(tabLoc_4 + 1);
	tmpLength = atoi(tmpLengthStr.c_str());
	tmpEffectiveLength = atof(tmpEffectiveLengthStr.c_str());
	tmpTpm = atof(tmpTpmStr.c_str());
	tmpReadNum = atof(tmpReadNumStr.c_str());
}

void parseSailfishResultsStr_tpm_readNum_only(string& tmpStr, 
	double& tmpTpm, double& tmpReadNum)
{
	int tabLoc_1 = tmpStr.find("\t");
	int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
	string tmpTpmStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	string tmpReadNumStr = tmpStr.substr(tabLoc_4 + 1);
	tmpTpm = atof(tmpTpmStr.c_str());
	tmpReadNum = atof(tmpReadNumStr.c_str());	
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_sailfish_results_list_file" << endl;
		cout << "#2 input_id_list_file" << endl;
		cout << "#3 output_merged_sailfish_results_file_prefix" << endl;
		exit(1);
	}

	string input_sailfish_results_list_file = argv[1];
	string input_sampleId_list_file = argv[2];
	string output_merged_sailfish_results_file_prefix = argv[3];

	/*******************************************/
	/*******************************************/
	vector<string> sampleId_vec;
	vector<string> transcriptIdVec;
	vector< pair<int, double> > lengthPairVec; // <length, effective length>
	vector< vector<double> > tpmVecVec;
	vector< vector<double> > readNumVecVec;	
	/*******************************************/
	/*******************************************/

	cout << "start to get sampleId list" << endl;
	ifstream sampleIdList_ifs(input_sampleId_list_file.c_str());
	while(!sampleIdList_ifs.eof())
	{
		string tmpStr;
		getline(sampleIdList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		sampleId_vec.push_back(tmpStr);
	}
	sampleIdList_ifs.close();

	cout << "start to get sailfish results file list" << endl;
	ifstream sailfishList_ifs(input_sailfish_results_list_file.c_str());
	vector<string> sailfish_results_file_vec;
	while(!sailfishList_ifs.eof())
	{
		string tmpStr;
		getline(sailfishList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		sailfish_results_file_vec.push_back(tmpStr);
	}
	sailfishList_ifs.close();

	cout << "start to get transcriptId and lengthVec (and results) from the 1st sailfish results file, " << endl;
	//cout << "sailfish_results_file_vec.size(): " << sailfish_results_file_vec.size() << endl;
	string sailfishResults_file_1st = sailfish_results_file_vec[0];
	//cout << "sailfishResults_file_1st: " << sailfishResults_file_1st << endl;
	ifstream sailfish_1st_ifs(sailfishResults_file_1st.c_str());
	string tmp1stLine;
	getline(sailfish_1st_ifs, tmp1stLine);
	vector<double> tpmVec_1st;
	vector<double> readNumVec_1st;
	while(!sailfish_1st_ifs.eof())
	{
		string tmpLine;
		getline(sailfish_1st_ifs, tmpLine);
		if(tmpLine == "")
			break;
		//cout << "tmpLine: " << tmpLine << endl;
		string id_tmp;
		int length_tmp;
		double effectiveLength_tmp, tpm_tmp, readNum_tmp;
		parseSailfishResultsStr(tmpLine, id_tmp, length_tmp, effectiveLength_tmp, tpm_tmp, readNum_tmp);
		transcriptIdVec.push_back(id_tmp);
		lengthPairVec.push_back(pair<int,double>(length_tmp, effectiveLength_tmp));
		tpmVec_1st.push_back(tpm_tmp);
		readNumVec_1st.push_back(readNum_tmp);
	}
	sailfish_1st_ifs.close();
	tpmVecVec.push_back(tpmVec_1st);
	readNumVecVec.push_back(readNumVec_1st);

	cout << "lengthPairVec.size(): " << lengthPairVec.size() << endl << endl;

	cout << "transcriptIdVec.size(): " << transcriptIdVec.size() << endl;
	cout << "sampleId_vec.size(): " << sampleId_vec.size() << endl;
	cout << "tpmVecVec.size(): " << tpmVecVec.size() << endl;
	cout << "tpmVecVec[0].size(): " << tpmVecVec[0].size() << endl;
	cout << "tpmVecVec[tpmVecVec.size()-1].size(): " << tpmVecVec[tpmVecVec.size()-1].size() << endl;
	cout << "readNumVecVec.size(): " << readNumVecVec.size() << endl;
	cout << "readNumVecVec[0].size(): " << readNumVecVec[0].size() << endl;
	cout << "readNumVecVec[readNumVecVec.size()-1].size(): " << readNumVecVec[readNumVecVec.size()-1].size() << endl;

	cout << "start to parse each sailfish results file" << endl;
	for(int tmp = 1; tmp < sailfish_results_file_vec.size(); tmp++)
	{
		cout << "file_index: " << tmp << endl;
		string tmp_sailfish_results_file = sailfish_results_file_vec[tmp];
		ifstream sailfish_tmp_ifs(tmp_sailfish_results_file.c_str());
		string tmpStr_1st;
		getline(sailfish_tmp_ifs, tmpStr_1st);
		vector<double> tpmVec_tmp;
		vector<double> readNumVec_tmp;		
		while(!sailfish_tmp_ifs.eof())
		{
			string tmpLine;
			getline(sailfish_tmp_ifs, tmpLine);
			if(tmpLine == "")
				break;
			double tpm_tmp, readNum_tmp;
			parseSailfishResultsStr_tpm_readNum_only(tmpLine, tpm_tmp, readNum_tmp);
			tpmVec_tmp.push_back(tpm_tmp);
			readNumVec_tmp.push_back(readNum_tmp);
		}
		cout << "tpmVec_tmp.size(): " << tpmVec_tmp.size() << endl;
		cout << "readNumVec_tmp.size(): " << readNumVec_tmp.size() << endl;
		tpmVecVec.push_back(tpmVec_tmp);
		readNumVecVec.push_back(readNumVec_tmp);
		sailfish_tmp_ifs.close();
	}


	cout << "start to merge multi sailfish results" << endl;
	
	string output_merged_sailfish_results_file_tpm = output_merged_sailfish_results_file_prefix + ".tpm";
	string output_merged_sailfish_results_file_readNum = output_merged_sailfish_results_file_prefix + ".readNum";
	ofstream merge_ofs_tpm(output_merged_sailfish_results_file_tpm.c_str());
	ofstream merge_ofs_readNum(output_merged_sailfish_results_file_readNum.c_str());
	//vector<string> sampleId_vec;
	//vector<string> transcriptIdVec;
	//vector< pair<int, double> > lengthPairVec; // <length, effective length>
	//vector< vector<double> > tpmVecVec;
	//vector< vector<double> > readNumVecVec;
	merge_ofs_tpm << "Transcript_id\\Sample_id";
	for(int tmp = 0; tmp < sampleId_vec.size(); tmp++)
		merge_ofs_tpm << "\t" << sampleId_vec[tmp];
	merge_ofs_tpm << endl;
	merge_ofs_readNum << "Transcript_id\\Sample_id";
	for(int tmp = 0; tmp < sampleId_vec.size(); tmp++)
		merge_ofs_readNum << "\t" << sampleId_vec[tmp];
	merge_ofs_readNum << endl;

	cout << "transcriptIdVec.size(): " << transcriptIdVec.size() << endl;
	cout << "sampleId_vec.size(): " << sampleId_vec.size() << endl;
	cout << "tpmVecVec.size(): " << tpmVecVec.size() << endl;
	cout << "tpmVecVec[0].size(): " << tpmVecVec[0].size() << endl;
	cout << "tpmVecVec[tpmVecVec.size()-1].size(): " << tpmVecVec[tpmVecVec.size()-1].size() << endl;
	cout << "readNumVecVec.size(): " << readNumVecVec.size() << endl;
	cout << "readNumVecVec[0].size(): " << readNumVecVec[0].size() << endl;
	cout << "readNumVecVec[readNumVecVec.size()-1].size(): " << readNumVecVec[readNumVecVec.size()-1].size() << endl;
	for(int tmpTranscript = 0; tmpTranscript < transcriptIdVec.size(); tmpTranscript++)
	{
		merge_ofs_tpm << transcriptIdVec[tmpTranscript] << "\t" << lengthPairVec[tmpTranscript].first 
			<< "\t" << lengthPairVec[tmpTranscript].second;
		merge_ofs_readNum << transcriptIdVec[tmpTranscript] << "\t" << lengthPairVec[tmpTranscript].first 
			<< "\t" << lengthPairVec[tmpTranscript].second;			
		for(int tmpSample = 0; tmpSample < sampleId_vec.size(); tmpSample++)
		{
			merge_ofs_tpm << "\t" << (tpmVecVec[tmpSample])[tmpTranscript];
			merge_ofs_readNum << "\t" << (readNumVecVec[tmpSample])[tmpTranscript];
		}
		merge_ofs_tpm << endl;
		merge_ofs_readNum << endl;
	}

	merge_ofs_tpm.close();
	merge_ofs_readNum.close();
	return 0;
}
