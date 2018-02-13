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

void generateSbtResults_rnaIdVec(vector<string>& sbtResults_rnaIdVec, string& rnaIdFile)
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
		sbtResults_rnaIdVec.push_back(tmpRnaId);
	}
	rnaId_ifs.close();
}

void generateSbtResults_foundSampleIdVecVec(
	vector< vector<int> >& sbtResults_foundSampleIdVecVec, string& sbtResults)
{
	ifstream sbt_ifs(sbtResults.c_str());
	while(!sbt_ifs.eof())
	{
		string tmpQuery_seq_hitNum;
		getline(sbt_ifs, tmpQuery_seq_hitNum);
		if(tmpQuery_seq_hitNum == "")
			break;
		int tabLoc = tmpQuery_seq_hitNum.find(" ");
		string tmpQuery_seq = tmpQuery_seq_hitNum.substr(0, tabLoc);
		string tmpQuery_hitNum_str = tmpQuery_seq_hitNum.substr(tabLoc + 1);
		int tmpHitNum = atoi(tmpQuery_hitNum_str.c_str());
		//cout << "tmpHitNum: " << tmpHitNum << endl;
		vector<int> tmpFoundSampleIdVec;
		if(tmpHitNum > 0)
		{
			for(int tmp = 0; tmp < tmpHitNum; tmp++)
			{
				string tmpHitSampleFilePath;
				getline(sbt_ifs, tmpHitSampleFilePath);
				int dot_bf_dot_bv_loc = tmpHitSampleFilePath.find(".bf.bv");
				string tmpHitSampleIdStr = tmpHitSampleFilePath.substr(44, dot_bf_dot_bv_loc - 44);
				//cout << "tmpHitSampleIdStr: " << tmpHitSampleIdStr << endl;
				int tmpHitSampleId = atoi(tmpHitSampleIdStr.c_str());
				tmpFoundSampleIdVec.push_back(tmpHitSampleId);
			}
		}
		sbtResults_foundSampleIdVecVec.push_back(tmpFoundSampleIdVec);
	}
	sbt_ifs.close();
}

void regenerateSbtExistOrNotVecVec(vector< vector<bool> >& sbtTranscriptExistOrNotBoolVecVec, 
	vector< vector<int> >& sbtResults_foundSampleIdVecVec, int rna_num, int sample_num)
{
	for(int tmp1 = 0; tmp1 < rna_num; tmp1++)
	{
		vector<bool> tmpSbtTranscriptExistOrNotBoolVec;
		for(int tmp2 = 0; tmp2 < sample_num; tmp2++)
			tmpSbtTranscriptExistOrNotBoolVec.push_back(false);
		sbtTranscriptExistOrNotBoolVecVec.push_back(tmpSbtTranscriptExistOrNotBoolVec);
	}
	for(int tmp1 = 0; tmp1 < rna_num; tmp1++)
	{
		int existedSampleNum = sbtResults_foundSampleIdVecVec[tmp1].size();
		for(int tmp2 = 0; tmp2 < existedSampleNum; tmp2++)
		{
			int tmpExistedSampleId = (sbtResults_foundSampleIdVecVec[tmp1])[tmp2];
			(sbtTranscriptExistOrNotBoolVecVec[tmp1])[tmpExistedSampleId - 1] = true;
		}
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
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 summarizedRsemResults" << endl; // /scratch/xli262/seqQueryTree/rsem/summarize/twoBatch.[TPM/FPKM/expectedCount]
		cout << "#2 rnaIdFile" << endl; // /scratch/xli262/seqQueryTree/transcriptSeq/transcript.rnaId2transcriptId.map
		cout << "#3 sbtResults" << endl; // example: /scratch/xli262/seqQueryTree/SBT_horse_query/transcript.query.t50
		cout << "#4 outputDir" << endl;
		cout << "#5 cutoff" << endl;
		exit(1);
	}
	string summarizedRsemResults = argv[1];
	string rnaIdFile = argv[2];
	string sbtResults = argv[3];
	string outputDir = argv[4]; outputDir += "/";
	string cutoffStr = argv[5];
	string mkdir_cmd = "mkdir " + outputDir;
	system(mkdir_cmd.c_str());

	// start to summarize RSEM TPM/FPKM/expectedCount results;
	vector<string> rsemRnaIdVec;
	vector< vector<double> > rsemAbundanceVecVec;
	generateRnaIdAndAbundanceVec_from_rsem(rsemRnaIdVec, rsemAbundanceVecVec, summarizedRsemResults);

	int rna_num = rsemAbundanceVecVec.size();
	int sample_num = rsemAbundanceVecVec[0].size();

	// start to generate sbtResults_rnaIdVec and sbtResults_foundSampleIdVecVec
	vector<string> sbtResults_rnaIdVec;
	generateSbtResults_rnaIdVec(sbtResults_rnaIdVec, rnaIdFile);
	vector< vector<int> > sbtResults_foundSampleIdVecVec;
	generateSbtResults_foundSampleIdVecVec(sbtResults_foundSampleIdVecVec, sbtResults);
	vector< vector<bool> > sbtTranscriptExistOrNotBoolVecVec;
	regenerateSbtExistOrNotVecVec(sbtTranscriptExistOrNotBoolVecVec, 
		sbtResults_foundSampleIdVecVec, rna_num, sample_num);

	// reorder rsem results according to rnaIds in sbtResults
	//vector<string> rsemRnaIdVec_reordered;
	double cutoff = atof(cutoffStr.c_str());
	vector< vector<bool> > rsemTranscriptExistOrNotBoolVecVec;
	regenerateRsemExistOrNotVecVec(rsemTranscriptExistOrNotBoolVecVec, rsemAbundanceVecVec, cutoff);

	int rsem_exist_sbt_exist_num = 0;
	int rsem_nonexist_sbt_nonexist_num = 0;
	int rsem_exist_sbt_nonexist_num = 0;
	int rsem_nonexist_sbt_exist_num = 0;

	string rsem_exist_sbt_exist_file = outputDir + "rsem_exist_sbt_exist_rna_sample.txt";
	string rsem_exist_sbt_nonexist_file = outputDir + "rsem_exist_sbt_nonexist_rna_sample.txt";
	string rsem_nonexist_sbt_exist_file = outputDir + "rsem_nonexist_sbt_exist_rna_sample.txt";
	string rsem_nonexist_sbt_nonexist_file = outputDir + "rsem_nonexist_sbt_nonexist_rna_sample.txt";
	ofstream rsem_exist_sbt_exist_ofs(rsem_exist_sbt_exist_file.c_str());
	ofstream rsem_exist_sbt_nonexist_ofs(rsem_exist_sbt_nonexist_file.c_str());
	ofstream rsem_nonexist_sbt_exist_ofs(rsem_nonexist_sbt_exist_file.c_str());
	ofstream rsem_nonexist_sbt_nonexist_ofs(rsem_nonexist_sbt_nonexist_file.c_str());
	for(int tmp1 = 0; tmp1 < rna_num; tmp1++)
	{
		for(int tmp2 = 0; tmp2 < sample_num; tmp2++)
		{
			bool tmp_rsem_exist_or_not = (rsemTranscriptExistOrNotBoolVecVec[tmp1])[tmp2];
			bool tmp_sbt_exist_or_not = (sbtTranscriptExistOrNotBoolVecVec[tmp1])[tmp2];
			if(tmp_rsem_exist_or_not && tmp_sbt_exist_or_not)
			{
				rsem_exist_sbt_exist_num ++;
				rsem_exist_sbt_exist_ofs << sbtResults_rnaIdVec[tmp1] << "\t" << tmp2 + 1 << endl;
			}
			else if(tmp_rsem_exist_or_not && (!tmp_sbt_exist_or_not))
			{
				rsem_exist_sbt_nonexist_num ++;
				rsem_exist_sbt_nonexist_ofs << sbtResults_rnaIdVec[tmp1] << "\t" << tmp2 + 1 << endl;
			}
			else if((!tmp_rsem_exist_or_not) && tmp_sbt_exist_or_not)
			{
				rsem_nonexist_sbt_exist_num ++;
				rsem_nonexist_sbt_exist_ofs << sbtResults_rnaIdVec[tmp1] << "\t" << tmp2 + 1 << endl;
			}
			else
			{
				rsem_nonexist_sbt_nonexist_num ++;
				rsem_nonexist_sbt_nonexist_ofs << sbtResults_rnaIdVec[tmp1] << "\t" << tmp2 + 1 << endl;
			}
		}
	}
	rsem_exist_sbt_exist_ofs.close();
	rsem_exist_sbt_nonexist_ofs.close();
	rsem_nonexist_sbt_exist_ofs.close();
	rsem_nonexist_sbt_nonexist_ofs.close();

	int total_num = rsem_exist_sbt_exist_num + rsem_nonexist_sbt_nonexist_num 
		+ rsem_exist_sbt_nonexist_num + rsem_nonexist_sbt_exist_num;
	string cmp_file = outputDir + "cmp.txt";
	ofstream cmp_ofs(cmp_file.c_str());
	cmp_ofs << "total_num:\t" << total_num << endl << endl;
	cmp_ofs << "rsem_exist_sbt_exist_num:\t" << rsem_exist_sbt_exist_num << endl;
	cmp_ofs << "rsem_exist_sbt_nonexist_num:\t" << rsem_exist_sbt_nonexist_num << endl;
	cmp_ofs << "rsem_nonexist_sbt_exist_num:\t" << rsem_nonexist_sbt_exist_num << endl;	
	cmp_ofs << "rsem_nonexist_sbt_nonexist_num:\t" << rsem_nonexist_sbt_nonexist_num << endl;
	cmp_ofs << endl;
	cmp_ofs << "rsem_exist_sbt_exist_perc:\t" << ((double)rsem_exist_sbt_exist_num/(double)total_num)*100 << "%" << endl;
	cmp_ofs << "rsem_exist_sbt_nonexist_perc:\t" << ((double)rsem_exist_sbt_nonexist_num/(double)total_num)*100 << "%" << endl;
	cmp_ofs << "rsem_nonexist_sbt_exist_perc:\t" << ((double)rsem_nonexist_sbt_exist_num/(double)total_num)*100 << "%" << endl;	
	cmp_ofs << "rsem_nonexist_sbt_nonexist_perc:\t" << ((double)rsem_nonexist_sbt_nonexist_num/(double)total_num)*100 << "%" << endl;
	cmp_ofs << endl;
	cmp_ofs << endl;
	int true_num = rsem_exist_sbt_exist_num + rsem_nonexist_sbt_nonexist_num;
	int false_num = rsem_exist_sbt_nonexist_num + rsem_nonexist_sbt_exist_num;
	cmp_ofs << "true_num:\t" << true_num << endl;
	cmp_ofs << "false_num:\t" << false_num << endl;
	cmp_ofs << endl;
	cmp_ofs << "true_perc:\t" << ((double)true_num/(double)total_num)*100 << "%" << endl;
	cmp_ofs << "false_perc:\t" << ((double)false_num/(double)total_num)*100 << "%" << endl;
	int tp_num = rsem_exist_sbt_exist_num;
	int fn_num = rsem_exist_sbt_nonexist_num;
	int fp_num = rsem_nonexist_sbt_exist_num;
	int tn_num = rsem_nonexist_sbt_nonexist_num;

	cmp_ofs << endl;
	cmp_ofs << "true_positive_rate:\t" << ((double)tp_num/(double)(tp_num + fn_num))*100 << "%" << endl;
	cmp_ofs << "false_positive_rate:\t" << ((double)fp_num/(double)(fp_num + tn_num))*100 << "%" << endl;
	cmp_ofs.close();
	return 0;
}