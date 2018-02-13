// merged sailfish results
// Transcript_id\sample_id sample_1 sample_2 ... sample_N
// Transcript_1 tpm(readNum)_1_1 tpm(readNum)_1_2 ... tpm(readNum)_1_N  
// Transcript_2 tpm(readNum)_2_1 tpm(readNum)_2_2 ... tpm(readNum)_2_N
// ....
// Transcript_N tpm(readNum)_N_1 tpm(readNum)_N_2 ... tpm(readNum)_N_N 
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

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputMergedSailfishResults" << endl;
		cout << "#2 sample_id" << endl;
		cout << "#3 transcript_id" << endl;
		exit(1);
	}
	string inputMergedSailfishResults = argv[1];
	string sample_id = argv[2];
	string transcript_id = argv[3];

	vector<string> sampleIdVec; 
	vector<string> transcriptIdVec;
	vector< vector<double> > expVecVec;	

	cout << "Start to parse merged sailfish results ..." << endl;
	parseMergedSailfishResults(inputMergedSailfishResults, 
		sampleIdVec, transcriptIdVec, expVecVec);

	cout << "Start to query ..." << endl;
	double tmpExp = queryExp_transcriptId_sampleId(sampleIdVec, 
		transcriptIdVec, expVecVec, sample_id, transcript_id);

	cout << "Sample_id: " << sample_id << endl;
	cout << "Transcript_id: " << transcript_id << endl;
	cout << "Exp: " << tmpExp << endl;

	return 0;
}