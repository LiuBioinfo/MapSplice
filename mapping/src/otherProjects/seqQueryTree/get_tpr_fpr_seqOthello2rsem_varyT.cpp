#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

void generateRnaIdAndAbundanceVec_from_rsem(vector<string>& rnaIdVec, 
	vector< vector<double> >& abundanceVecVec, string& summarizedRsemResults)
{
	ifstream rsem_ifs(summarizedRsemResults.c_str());
	string tmpHeader_rsem;
	getline(rsem_ifs, tmpHeader_rsem);
	while(!rsem_ifs.eof())
	{
		string tmpStr;
		getline(rsem_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpStartLoc = 0;
		vector<string> tmpFieldVec;
		for(int tmp = 0; ; tmp++)
		{
			int tmpTabLoc = tmpStr.find("\t", tmpStartLoc);
			if(tmpTabLoc == string::npos)
				break;
			string tmpField = tmpStr.substr(tmpStartLoc, tmpTabLoc - tmpStartLoc);
			tmpFieldVec.push_back(tmpField);
			tmpStartLoc = tmpTabLoc + 1;
		}
		tmpFieldVec.push_back(tmpStr.substr(tmpStartLoc));
		rnaIdVec.push_back(tmpFieldVec[0]);
		vector<double> tmpAbundanceVec;
		for(int tmp = 1; tmp < tmpFieldVec.size(); tmp++)
		{
			double tmpAbundance = atof(tmpFieldVec[tmp].c_str());
			tmpAbundanceVec.push_back(tmpAbundance);
		}
		abundanceVecVec.push_back(tmpAbundanceVec);
	}
	rsem_ifs.close();
}

void generateSeqOthelloResults_rnaIdVec(vector<string>& seqOthelloResults_rnaIdVec, string& rnaIdFile)
{
	ifstream rnaId_ifs(rnaIdFile.c_str());
	while(!rnaId_ifs.eof())
	{
		string tmpStr;
		getline(rnaId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpRnaId = tmpStr.substr(0, tabLoc);
		seqOthelloResults_rnaIdVec.push_back(tmpRnaId);
	}
	rnaId_ifs.close();
}

void generateSeqOthelloResults_sampleExistKmerNumVecVec(int sample_num, vector<int>& seqLengthVec, 
	vector< vector<int> >& seqOthelloResults_sampleExistKmerNumVecVec, string& seqOthelloResults)
{
	ifstream seqOthello_ifs(seqOthelloResults.c_str());
	while(!seqOthello_ifs.eof())
	{
		string tmpQuery_seq;
		getline(seqOthello_ifs, tmpQuery_seq);
		if(tmpQuery_seq == "")
			break;
		int tmpQuery_seq_length = tmpQuery_seq.length();
		seqLengthVec.push_back(tmpQuery_seq_length);
		vector<int> tmpSampleExistKmerNumVec;
		for(int tmp = 0; tmp < sample_num; tmp++)
		{
			string tmpKmerNumInTmpSampleStr;
			getline(seqOthello_ifs, tmpKmerNumInTmpSampleStr);
			int tmpKmerNumInTmpSampleInt = atoi(tmpKmerNumInTmpSampleStr.c_str());
			tmpSampleExistKmerNumVec.push_back(tmpKmerNumInTmpSampleInt);
		}
		string tmpStr_1, tmpStr_2, tmpStr_3, tmpStr_4;
		getline(seqOthello_ifs, tmpStr_1);//149
		getline(seqOthello_ifs, tmpStr_2);//150
		getline(seqOthello_ifs, tmpStr_3);//151
		getline(seqOthello_ifs, tmpStr_4);//152
		seqOthelloResults_sampleExistKmerNumVecVec.push_back(tmpSampleExistKmerNumVec);
	}
	seqOthello_ifs.close();
}

void regenerateSeqOthelloExistOrNotVecVecVec(vector< vector< vector<bool> > >& seqOthelloTranscriptExistOrNotBoolVecVecVec,
	vector<int>& seqLengthVec, vector< vector<int> >& seqOthelloResults_sampleExistKmerNumVecVec, 
	int rna_num, int sample_num, vector<double>& cutoff_existKmerPercMin_vec)
{
	for(int tmp = 0; tmp < cutoff_existKmerPercMin_vec.size(); tmp ++)
	{
		double tmp_cutoff_existKmerPercMin = cutoff_existKmerPercMin_vec[tmp];
		vector< vector<bool> > tmpSeqOthelloTranscriptExistOrNotBoolVecVec;
		for(int tmp1 = 0; tmp1 < rna_num; tmp1++)
		{
			//cout << "tmp1: " << tmp1 << endl;
			int tmpRnaLength = seqLengthVec[tmp1];
			//cout << "tmpRnaLength: " << tmpRnaLength << endl;
			vector<bool> tmpSeqOthelloTranscriptExistOrNotBoolVec;
			for(int tmp2 = 0; tmp2 < sample_num; tmp2++)
			{
				//cout << "tmp2: " << tmp2 << endl;
				double tmpExistKmerPerc = (double)(seqOthelloResults_sampleExistKmerNumVecVec[tmp1])[tmp2]/(double)(tmpRnaLength-19);
				//cout << "tmpExistKmerPerc: " << tmpExistKmerPerc << endl;
				if(tmpExistKmerPerc >= tmp_cutoff_existKmerPercMin)	
					tmpSeqOthelloTranscriptExistOrNotBoolVec.push_back(true);
				else
					tmpSeqOthelloTranscriptExistOrNotBoolVec.push_back(false);
			}
			tmpSeqOthelloTranscriptExistOrNotBoolVecVec.push_back(tmpSeqOthelloTranscriptExistOrNotBoolVec);
		}
		seqOthelloTranscriptExistOrNotBoolVecVecVec.push_back(tmpSeqOthelloTranscriptExistOrNotBoolVecVec);
	}
}

void regenerateRsemExistOrNotVecVec(vector< vector<bool> >& rsemTranscriptExistOrNotBoolVecVec, 
	vector< vector<double> >& rsemAbundanceVecVec, double cutoff)
{
	int rnaNum = rsemAbundanceVecVec.size();
	int sampleNum = rsemAbundanceVecVec[0].size();
	for(int tmp1 = 0; tmp1 < rnaNum; tmp1++)
	{
		vector<bool> tmpRsemTranscriptExistOrNotBoolVec;
		for(int tmp2 = 0; tmp2 < sampleNum; tmp2++)
			tmpRsemTranscriptExistOrNotBoolVec.push_back(false);
		rsemTranscriptExistOrNotBoolVecVec.push_back(tmpRsemTranscriptExistOrNotBoolVec);
	}

	for(int tmp1 = 0; tmp1 < rnaNum; tmp1++)
	{
		for(int tmp2 = 0; tmp2 < sampleNum; tmp2++)
		{
			double tmpAbundance = (rsemAbundanceVecVec[tmp1])[tmp2];
			if(tmpAbundance >= cutoff)
				(rsemTranscriptExistOrNotBoolVecVec[tmp1])[tmp2] = true;
			else
				(rsemTranscriptExistOrNotBoolVecVec[tmp1])[tmp2] = false;
		}
	}
}

int main(int argc, char** argv)
{
	if(argc < 7)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 summarizedRsemResults" << endl; // /scratch/xli262/seqQueryTree/rsem/summarize/twoBatch.[TPM/FPKM/expectedCount]
		cout << "#2 rnaIdFile" << endl; // /scratch/xli262/seqQueryTree/transcriptSeq/transcript.rnaId2transcriptId.map
		cout << "#3 seqOthelloResults" << endl; // example: /scratch/xli262/seqQueryTree/SBT_horse_query/transcript.query.t50
		cout << "#4 outputDir" << endl;
		cout << "#5 cutoff_rsemTpm" << endl;
		cout << "#6 cutoffListFile_existKmerPercMin" << endl;
		cout << "#7 toolName" << endl;
		exit(1);
	}
	string summarizedRsemResults = argv[1];
	string rnaIdFile = argv[2];
	string seqOthelloResults = argv[3];
	string outputDir = argv[4]; outputDir += "/";
	string cutoffStr_rsemTpm = argv[5];
	string cutoffListFile_existKmerPercMin = argv[6];
	string toolName = argv[7];
	vector<double> cutoff_existKmerPercMin_vec;
	ifstream existKmerPercMinCutoffList_ifs(cutoffListFile_existKmerPercMin.c_str());
	while(!existKmerPercMinCutoffList_ifs.eof())
	{
		string tmpStr;
		getline(existKmerPercMinCutoffList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		double tmpExistKmerPercMinCutoff = atof(tmpStr.c_str());
		cutoff_existKmerPercMin_vec.push_back(tmpExistKmerPercMinCutoff);
	}
	existKmerPercMinCutoffList_ifs.close();

	string mkdir_cmd = "mkdir " + outputDir;
	system(mkdir_cmd.c_str());

	// start to summarize RSEM TPM/FPKM/expectedCount results;
	cout << "start to summarize RSEM TPM/FPKM/expectedCount results" << endl;
	vector<string> rsemRnaIdVec;
	vector< vector<double> > rsemAbundanceVecVec;
	generateRnaIdAndAbundanceVec_from_rsem(rsemRnaIdVec, rsemAbundanceVecVec, summarizedRsemResults);

	int rna_num = rsemAbundanceVecVec.size();
	int sample_num = rsemAbundanceVecVec[0].size();

	// start to generate seqOthelloResults_rnaIdVec and seqOthelloResults_foundSampleIdVecVec
	cout << "start to start to generate seqOthelloResults_rnaIdVec and seqOthelloResults_foundSampleIdVecVec" << endl;
	vector<string> seqOthelloResults_rnaIdVec;
	generateSeqOthelloResults_rnaIdVec(seqOthelloResults_rnaIdVec, rnaIdFile);

	vector<int> seqLengthVec;
	vector< vector<int> > seqOthelloResults_sampleExistKmerNumVecVec;
	cout << "start to start to generateSeqOthelloResults_sampleExistKmerNumVecVec" << endl;
	generateSeqOthelloResults_sampleExistKmerNumVecVec(sample_num, seqLengthVec, 
		seqOthelloResults_sampleExistKmerNumVecVec, seqOthelloResults);

	vector< vector< vector<bool> > > seqOthelloTranscriptExistOrNotBoolVecVecVec;
	//double cutoff_existKmerPercMin = atof(cutoffStr_existKmerPercMin.c_str());
	cout << "start to regenerateSeqOthelloExistOrNotVecVecVec" << endl;
	regenerateSeqOthelloExistOrNotVecVecVec(seqOthelloTranscriptExistOrNotBoolVecVecVec, seqLengthVec, 
		seqOthelloResults_sampleExistKmerNumVecVec, rna_num, sample_num, cutoff_existKmerPercMin_vec);

	// reorder rsem results according to rnaIds in seqOthelloResults
	// vector<string> rsemRnaIdVec_reordered;
	double cutoff_rsemTpm = atof(cutoffStr_rsemTpm.c_str());
	vector< vector<bool> > rsemTranscriptExistOrNotBoolVecVec;
	cout << "start to regenerateRsemExistOrNotVecVec" << endl;
	regenerateRsemExistOrNotVecVec(rsemTranscriptExistOrNotBoolVecVec, rsemAbundanceVecVec, cutoff_rsemTpm);

	cout << "start to print results" << endl;
	vector<int> rsem_exist_seqOthello_exist_num_vec;
	vector<int> rsem_exist_seqOthello_nonexist_num_vec;
	vector<int> rsem_nonexist_seqOthello_exist_num_vec;
	vector<int> rsem_nonexist_seqOthello_nonexist_num_vec;
	for(int tmp = 0; tmp < cutoff_existKmerPercMin_vec.size(); tmp++)
	{	
		int tmp_rsem_exist_seqOthello_exist_num = 0;
		int tmp_rsem_nonexist_seqOthello_nonexist_num = 0;
		int tmp_rsem_exist_seqOthello_nonexist_num = 0;
		int tmp_rsem_nonexist_seqOthello_exist_num = 0;
		for(int tmp1 = 0; tmp1 < rna_num; tmp1++)
		{
			for(int tmp2 = 0; tmp2 < sample_num; tmp2++)
			{
				bool tmp_rsem_exist_or_not = (rsemTranscriptExistOrNotBoolVecVec[tmp1])[tmp2];
				bool tmp_seqOthello_exist_or_not = ((seqOthelloTranscriptExistOrNotBoolVecVecVec[tmp])[tmp1])[tmp2];
				if(tmp_rsem_exist_or_not && tmp_seqOthello_exist_or_not)
					tmp_rsem_exist_seqOthello_exist_num ++;
				else if(tmp_rsem_exist_or_not && (!tmp_seqOthello_exist_or_not))
					tmp_rsem_exist_seqOthello_nonexist_num ++;
				else if((!tmp_rsem_exist_or_not) && tmp_seqOthello_exist_or_not)
					tmp_rsem_nonexist_seqOthello_exist_num ++;
				else
					tmp_rsem_nonexist_seqOthello_nonexist_num ++;
			}
		}
		rsem_exist_seqOthello_exist_num_vec.push_back(tmp_rsem_exist_seqOthello_exist_num);
		rsem_exist_seqOthello_nonexist_num_vec.push_back(tmp_rsem_exist_seqOthello_nonexist_num);
		rsem_nonexist_seqOthello_exist_num_vec.push_back(tmp_rsem_nonexist_seqOthello_exist_num);
		rsem_nonexist_seqOthello_nonexist_num_vec.push_back(tmp_rsem_nonexist_seqOthello_nonexist_num);
	}

	string cmp_file = outputDir + "cmp.txt";
	string cmp_file_tpr_fpr = outputDir + "cmp_roc.txt";
	ofstream cmp_ofs(cmp_file.c_str());
	ofstream cmp_tpr_fpr_ofs(cmp_file_tpr_fpr.c_str());
	for(int tmp = 0; tmp < cutoff_existKmerPercMin_vec.size(); tmp++)
	{	
		double tmp_cutoff_existKmerPercMin = cutoff_existKmerPercMin_vec[tmp];
		int rsem_exist_seqOthello_exist_num = rsem_exist_seqOthello_exist_num_vec[tmp];
		int rsem_nonexist_seqOthello_nonexist_num = rsem_nonexist_seqOthello_nonexist_num_vec[tmp];
		int rsem_exist_seqOthello_nonexist_num = rsem_exist_seqOthello_nonexist_num_vec[tmp];
		int rsem_nonexist_seqOthello_exist_num = rsem_nonexist_seqOthello_exist_num_vec[tmp];
		int total_num = rsem_exist_seqOthello_exist_num + rsem_nonexist_seqOthello_nonexist_num 
			+ rsem_exist_seqOthello_nonexist_num + rsem_nonexist_seqOthello_exist_num;
		cmp_ofs << "Cutoff_existKmerPercMin:\t" << tmp_cutoff_existKmerPercMin << endl;
		cmp_ofs << "total_num:\t" << total_num << endl << endl;
		cmp_ofs << "rsem_exist_seqOthello_exist_num:\t" << rsem_exist_seqOthello_exist_num << endl;
		cmp_ofs << "rsem_exist_seqOthello_nonexist_num:\t" << rsem_exist_seqOthello_nonexist_num << endl;
		cmp_ofs << "rsem_nonexist_seqOthello_exist_num:\t" << rsem_nonexist_seqOthello_exist_num << endl;	
		cmp_ofs << "rsem_nonexist_seqOthello_nonexist_num:\t" << rsem_nonexist_seqOthello_nonexist_num << endl;
		cmp_ofs << endl;
		cmp_ofs << "rsem_exist_seqOthello_exist_perc:\t" << ((double)rsem_exist_seqOthello_exist_num/(double)total_num)*100 << "%" << endl;
		cmp_ofs << "rsem_exist_seqOthello_nonexist_perc:\t" << ((double)rsem_exist_seqOthello_nonexist_num/(double)total_num)*100 << "%" << endl;
		cmp_ofs << "rsem_nonexist_seqOthello_exist_perc:\t" << ((double)rsem_nonexist_seqOthello_exist_num/(double)total_num)*100 << "%" << endl;	
		cmp_ofs << "rsem_nonexist_seqOthello_nonexist_perc:\t" << ((double)rsem_nonexist_seqOthello_nonexist_num/(double)total_num)*100 << "%" << endl;
		cmp_ofs << endl;
		cmp_ofs << endl;
		int true_num = rsem_exist_seqOthello_exist_num + rsem_nonexist_seqOthello_nonexist_num;
		int false_num = rsem_exist_seqOthello_nonexist_num + rsem_nonexist_seqOthello_exist_num;
		cmp_ofs << "true_num:\t" << true_num << endl;
		cmp_ofs << "false_num:\t" << false_num << endl;
		cmp_ofs << endl;
		cmp_ofs << "true_perc:\t" << ((double)true_num/(double)total_num)*100 << "%" << endl;
		cmp_ofs << "false_perc:\t" << ((double)false_num/(double)total_num)*100 << "%" << endl;
		int tp_num = rsem_exist_seqOthello_exist_num;
		int fn_num = rsem_exist_seqOthello_nonexist_num;
		int fp_num = rsem_nonexist_seqOthello_exist_num;
		int tn_num = rsem_nonexist_seqOthello_nonexist_num;
		cmp_ofs << endl;
		cmp_ofs << "true_positive_rate:\t" << ((double)tp_num/(double)(tp_num + fn_num))*100 << "%" << endl;
		cmp_ofs << "false_positive_rate:\t" << ((double)fp_num/(double)(fp_num + tn_num))*100 << "%" << endl << endl;
		cmp_tpr_fpr_ofs << tmp_cutoff_existKmerPercMin << "\t" << ((double)tp_num/(double)(tp_num + fn_num)) 
			<< "\t" << ((double)fp_num/(double)(fp_num + tn_num)) << "\ttpm" << cutoff_rsemTpm << "\t" << toolName << endl;
	}
	cmp_tpr_fpr_ofs.close();
	cmp_ofs.close();
	return 0;
}