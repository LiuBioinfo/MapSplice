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

void generateRnaId2transcriptIdVec(
	vector< pair<string, string> >& rnaId2transcriptIdVec,
	string& input_gff3transcriptIdFile)
{
	ifstream gff3transcriptId_ifs(input_gff3transcriptIdFile.c_str());
	while(!gff3transcriptId_ifs.eof())
	{
		string tmpStr;
		getline(gff3transcriptId_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpLoc_rnaId = tmpStr.find("ID=");
		int tmpLoc_rnaId_nextComma = tmpStr.find(";", tmpLoc_rnaId + 1);
		int tmpLoc_transcriptId = tmpStr.find("transcript_id=");
		string tmpRNAid = tmpStr.substr(tmpLoc_rnaId + 3, 
			tmpLoc_rnaId_nextComma - 1 - (tmpLoc_rnaId + 3) + 1);
		string tmpTranscriptId = tmpStr.substr(tmpLoc_transcriptId + 14);
		rnaId2transcriptIdVec.push_back(pair<string,string>(tmpRNAid, tmpTranscriptId));
	}
	gff3transcriptId_ifs.close();
}

void generateFaTranscriptIdSeqVec(
	vector< pair<string, string> >& faTranscriptIdSeqVec, string& input_faFile)
{
	vector<string> faLineVec;
	ifstream fa_ifs(input_faFile.c_str());
	while(!fa_ifs.eof())
	{
		string tmpStr;
		getline(fa_ifs, tmpStr);
		if(tmpStr == "")
			break;
		faLineVec.push_back(tmpStr);
	}
	fa_ifs.close();
	int tmpVecIndex = -1;
	for(int tmp = 0; tmp < faLineVec.size(); tmp++)
	{
		string tmpStr = faLineVec[tmp];
		if(tmpStr.at(0) == '>')
		{
			//string tmpStr;
			int tmpLoc_ref = tmpStr.find("ref|");
			int tmpLoc_nextStraightLine = tmpStr.find("|", tmpLoc_ref + 4);
			string tmpId = tmpStr.substr(tmpLoc_ref + 4, 
				tmpLoc_nextStraightLine - 1 - (tmpLoc_ref + 4) + 1);
			string tmpSeq = "";
			faTranscriptIdSeqVec.push_back(pair<string,string>(tmpId, tmpSeq));
			tmpVecIndex ++;
		}
		else
			(faTranscriptIdSeqVec[tmpVecIndex].second) += tmpStr;
	}
}

void generateRsemRnaIdVec(vector<string>& rsemRnaIdVec, string& input_rsemFile)
{
	ifstream rsem_ifs(input_rsemFile.c_str());
	string tmpHeader;
	getline(rsem_ifs, tmpHeader);
	while(!rsem_ifs.eof())
	{
		string tmpStr;
		getline(rsem_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpLoc = tmpStr.find("\t");
		string tmpRnaId = tmpStr.substr(0, tmpLoc);
		rsemRnaIdVec.push_back(tmpRnaId);
	}
	rsem_ifs.close();
}

bool searchRnaIdAndReturnTranscriptId(string& tmpRsemRnaId, string& tmpTranscriptIdInGff3, 
	vector< pair<string, string> >& rnaId2transcriptIdVec)
{
	for(int tmp = 0; tmp < rnaId2transcriptIdVec.size(); tmp++)
	{
		string tmpRnaId = rnaId2transcriptIdVec[tmp].first;
		if(tmpRsemRnaId == tmpRnaId)
		{
			tmpTranscriptIdInGff3 = rnaId2transcriptIdVec[tmp].second;
			return true;
		}
	}
	return false;
}

bool searchTranscriptIdAndReturnSeq(string& tmpTranscriptIdInGff3, string& tmpSeq, 
	vector< pair<string, string> >& faTranscriptIdSeqVec)
{
	for(int tmp = 0; tmp < faTranscriptIdSeqVec.size(); tmp++)
	{
		string tmpTranscriptId = faTranscriptIdSeqVec[tmp].first;
		if(tmpTranscriptIdInGff3 == tmpTranscriptId)
		{
			tmpSeq = faTranscriptIdSeqVec[tmp].second;
			return true;
		}
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 input_gff3transcriptIdFile" << endl; // /scratch/xli262/seqQueryTree/gff3.rnaId.transcriptId.idField
		cout << "#2 input_faFile" << endl; // /scratch/xli262/seqQueryTree/rna.fa
		cout << "#3 input_rsemFile" << endl; // /scratch/xli262/seqQueryTree/rsem/batch_2/1.isoforms.results
		cout << "#4 output_file_prefix" << endl;
		exit(1);
	}
	string input_gff3transcriptIdFile = argv[1];
	string input_faFile = argv[2];
	string input_rsemFile = argv[3];
	string output_file_prefix = argv[4];
	string output_fa_file = output_file_prefix + ".fa";
	string output_seq_file = output_file_prefix + ".seq";
	string output_info_file = output_file_prefix + ".rnaId2transcriptId.map";

	vector< pair<string, string> > rnaId2transcriptIdVec;
	generateRnaId2transcriptIdVec(rnaId2transcriptIdVec, input_gff3transcriptIdFile);
	cout << "rnaId2transcriptIdVec.size(): " << rnaId2transcriptIdVec.size() << endl;
	vector< pair<string, string> > faTranscriptIdSeqVec;
	generateFaTranscriptIdSeqVec(faTranscriptIdSeqVec, input_faFile);
	cout << "faTranscriptIdSeqVec.size(): " << faTranscriptIdSeqVec.size() << endl;
	vector<string> rsemRnaIdVec;
	generateRsemRnaIdVec(rsemRnaIdVec, input_rsemFile);
	cout << "rsemRnaIdVec.size(): " << rsemRnaIdVec.size() << endl;

	ofstream fa_ofs(output_fa_file.c_str());
	ofstream seq_ofs(output_seq_file.c_str());
	ofstream info_ofs(output_info_file.c_str());
	for(int tmp = 0; tmp < rsemRnaIdVec.size(); tmp++)
	{
		string tmpRsemRnaId = rsemRnaIdVec[tmp];
		string tmpTranscriptIdInGff3;
		bool foundInRnaId2transcriptIdVec_bool = searchRnaIdAndReturnTranscriptId(
			tmpRsemRnaId, tmpTranscriptIdInGff3, rnaId2transcriptIdVec);
		if(!foundInRnaId2transcriptIdVec_bool)
			continue;
		string tmpSeq;
		bool foundInFaVec_bool = searchTranscriptIdAndReturnSeq(
			tmpTranscriptIdInGff3, tmpSeq, faTranscriptIdSeqVec);
		if(!foundInFaVec_bool)
			continue;
		fa_ofs << ">" << tmpRsemRnaId << "_" << tmpTranscriptIdInGff3 << endl << tmpSeq << endl;
		seq_ofs << tmpSeq << endl;
		info_ofs << tmpRsemRnaId << "\t" << tmpTranscriptIdInGff3 << endl;
	}
	fa_ofs.close();
	seq_ofs.close();
	info_ofs.close();
	return 0;
}