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
//#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;


// inline char getCharRevComp(char ch)
// {
// 	int chInt = ch - 'A';
// 	static const char alphatChar[26] = {'T', 'N', 'G', 'N', 'N', 'N', 'C',
// 		'N', 'N', 'N', 'N', 'N', 'N', 'N',
// 		'N', 'N', 'N', 'N', 'N', 'A',
// 		'N', 'N', 'N', 'N', 'N', 'N'};
// 	return alphatChar[chInt];
// }

// string getRcmSeq(const string& readSeq)
// {
// 	int readSeqLength = readSeq.length();

// 	char readRcmSeqChar[readSeqLength];

// 	readRcmSeqChar[0] = getCharRevComp(readSeq.at(readSeqLength-1));

// 	for(int tmp = 1; tmp < readSeqLength; tmp ++)
// 	{
// 		readRcmSeqChar[tmp] = getCharRevComp((readSeq.at(readSeqLength - tmp - 1)));
// 	}
// 	string rcmSeq = readRcmSeqChar;
// 	return rcmSeq.substr(0, readSeqLength);
// }

bool sortExonPosOrderBool(int i, int j) { return (i < j);}

int returnLeftMostExonStartSite(vector< int >& donerVec)
{
	int leftMost = donerVec[0];
	if(donerVec.size() > 1)
	{
		for(int tmp = 1; tmp < donerVec.size(); tmp++)
		{
			int tmpDonerSite = donerVec[tmp];
			if(tmpDonerSite < leftMost)
				leftMost = tmpDonerSite;
		}
	}
	return leftMost;
}

int returnRightMostExonEndSite(vector< int >& acceptorVec)
{
	int rightMost = acceptorVec[0];
	if(acceptorVec.size() > 1)
	{
		for(int tmp = 1; tmp < acceptorVec.size(); tmp++)
		{
			int tmpAcceptorSite = acceptorVec[tmp];
			if(tmpAcceptorSite > rightMost)
				rightMost = tmpAcceptorSite;
		}
	}
	return rightMost;
}

bool twoTranscriptSegTheSameOrNot(vector< pair<int,int> >& exonVec_1, vector< pair<int,int> >& exonVec_2)
{
	vector<int> startPosVec_1;
	vector<int> startPosVec_2;
	vector<int> endPosVec_1;
	vector<int> endPosVec_2;	
	for(int tmp = 0; tmp < exonVec_1.size(); tmp++)
	{
		startPosVec_1.push_back(exonVec_1[tmp].first);
		endPosVec_1.push_back(exonVec_1[tmp].second);
	}
	for(int tmp = 0; tmp < exonVec_2.size(); tmp++)
	{
		startPosVec_2.push_back(exonVec_2[tmp].first);
		endPosVec_2.push_back(exonVec_2[tmp].second);
	}
	sort(startPosVec_1.begin(), startPosVec_1.end(), sortExonPosOrderBool);
	sort(endPosVec_1.begin(), endPosVec_1.end(), sortExonPosOrderBool);
	sort(startPosVec_2.begin(), startPosVec_2.end(), sortExonPosOrderBool);
	sort(endPosVec_2.begin(), endPosVec_2.end(), sortExonPosOrderBool);		
	if(!((startPosVec_1.size() == startPosVec_2.size())&&(endPosVec_1.size() == endPosVec_2.size())
		&&(startPosVec_1.size() == endPosVec_1.size())))
		return false;
	int exonNum = startPosVec_1.size();
	for(int tmp = 0; tmp < exonNum; tmp++)
	{
		if(startPosVec_1[tmp] != startPosVec_2[tmp])
			return false;
		if(endPosVec_1[tmp] != endPosVec_2[tmp])
			return false;		
	}
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolder inputExonGTFfile inputBackSpliceJuncFile outputPrefix" << endl;
		exit(1);
	}
	int seqLength = 30;

	cout << "loading indexInfo parameters ......" << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	cout << "start to initiate annotated exonVec" << endl;
	string inputExonGTFfile = argv[2];
	vector< vector< pair<int,int> > > exonPosPairVecVec;
	vector< vector< string > > exonTranscriptIdVecVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		vector< pair<int,int> > tmpExonPosPairVec;
		vector< string > tmpExonTranscriptIdVec;
		exonPosPairVecVec.push_back(tmpExonPosPairVec);
		exonTranscriptIdVecVec.push_back(tmpExonTranscriptIdVec);
	}
	cout << "start to load annotated exons" << endl;
	ifstream exonGTF_ifs(inputExonGTFfile.c_str());
	while(!exonGTF_ifs.eof())
	{
		string tmpStr;
		getline(exonGTF_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		string tmpExon_chrName = tmpStr.substr(0, tabLoc_1);
		int tmpExon_chrNameInt = indexInfo->convertStringToInt(tmpExon_chrName);
		string tmpExon_startPosStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		int tmpExon_startPos = atoi(tmpExon_startPosStr.c_str());
		string tmpExon_endPosStr = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		int tmpExon_endPos = atoi(tmpExon_endPosStr.c_str());		
		int tmpTranscriptIdLoc = tmpStr.find("transcript_id");
		int tmpTranscriptIdLoc_end = tmpStr.find(";", tmpTranscriptIdLoc + 1);
		string tmpExon_transcriptId = tmpStr.substr(tmpTranscriptIdLoc + 15, tmpTranscriptIdLoc_end - 2 - tmpTranscriptIdLoc - 15 + 1);
		exonPosPairVecVec[tmpExon_chrNameInt].push_back(pair<int,int>(tmpExon_startPos, tmpExon_endPos));
		exonTranscriptIdVecVec[tmpExon_chrNameInt].push_back(tmpExon_transcriptId);
	}
	cout << "start to search for exons within back splice junctions" << endl;
	string outputPrefix = argv[4];
	string outputExonWithinBackSplice_log = outputPrefix + "_log.txt";
	string outputExonWithinBackSplice = outputPrefix + "_exonInfo.txt";
	string outputExonWithinBackSplice_seq = outputPrefix + "_sequence.txt";
	ofstream exonWithinBackSplice_log_ofs(outputExonWithinBackSplice_log.c_str());
	ofstream exonWithinBackSplice_ofs(outputExonWithinBackSplice.c_str());
	ofstream exonWithinBackSplice_seq_ofs(outputExonWithinBackSplice_seq.c_str());

	exonWithinBackSplice_seq_ofs << "chrName\tdonor_end\tacceptor_start\tflank_string\tgene\t";
	exonWithinBackSplice_seq_ofs << "genomicSeq_1_startPos\tgenomicSeq_1_endPos\tgenomicSeq_1\tgenomicSeq_1_reverseComplement";
	exonWithinBackSplice_seq_ofs << "\tgenomicSeq_2_startPos\tgenomicSeq_2_endPos\tgenomicSeq_2\tgenomicSeq_2_reverseComplement" << endl;

	string inputBackSpliceJuncFile = argv[3];
	ifstream backSplice_ifs(inputBackSpliceJuncFile.c_str());
	string firstLine;
	getline(backSplice_ifs, firstLine);
	while(!backSplice_ifs.eof())
	{
		string tmpStr;
		getline(backSplice_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		string tmpBackSplice_chrName = tmpStr.substr(0, tabLoc_1);
		string tmpBackSplice_donerStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpBackSplice_acceptorStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpBackSplice_flankString = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		string tmpBackSplice_gene = tmpStr.substr(tabLoc_4 + 1);

		int tmpBackSplice_chrNameInt = indexInfo->convertStringToInt(tmpBackSplice_chrName);
		int tmpBackSplice_doner = atoi(tmpBackSplice_donerStr.c_str());
		int tmpBackSplice_acceptor = atoi(tmpBackSplice_acceptorStr.c_str());

		int tmpBackSplice_doner_startPos = tmpBackSplice_doner - seqLength + 1;
		int tmpBackSplice_doner_endPos = tmpBackSplice_doner;
		int tmpBackSplice_acceptor_startPos = tmpBackSplice_acceptor;
		int tmpBackSplice_acceptor_endPos = tmpBackSplice_acceptor + seqLength - 1;

		vector<int> tmpWithinExonVec_start_doner;
		vector<int> tmpWithinExonVec_end_doner;
		vector<string> tmpWithinExonVec_transcriptId_doner;
		vector<int> tmpWithinExonVec_start_acceptor;
		vector<int> tmpWithinExonVec_end_acceptor;
		vector<string> tmpWithinExonVec_transcriptId_acceptor;

		int tmpChrExonNum = exonPosPairVecVec[tmpBackSplice_chrNameInt].size();
		for(int tmp = 0; tmp < tmpChrExonNum; tmp++)
		{
			int tmpChrExon_start = (exonPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].first;
			int tmpChrExon_end = (exonPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].second;
			if((tmpChrExon_start <= tmpBackSplice_doner_startPos)&&(tmpBackSplice_doner_endPos <= tmpChrExon_end))
			{
				tmpWithinExonVec_start_doner.push_back(tmpChrExon_start);
				tmpWithinExonVec_end_doner.push_back(tmpChrExon_end);
				tmpWithinExonVec_transcriptId_doner.push_back((exonTranscriptIdVecVec[tmpBackSplice_chrNameInt])[tmp]);
			}
			if((tmpChrExon_start <= tmpBackSplice_acceptor_startPos)&&(tmpBackSplice_acceptor_endPos <= tmpChrExon_end))
			{
				tmpWithinExonVec_start_acceptor.push_back(tmpChrExon_start);
				tmpWithinExonVec_end_acceptor.push_back(tmpChrExon_end);
				tmpWithinExonVec_transcriptId_acceptor.push_back((exonTranscriptIdVecVec[tmpBackSplice_chrNameInt])[tmp]);
			}
		}
		exonWithinBackSplice_ofs << "--------------------------------------------------------" << endl;
		exonWithinBackSplice_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr << "--" << tmpBackSplice_acceptorStr << ":" << endl;
		int tmpWithinExonNum_doner = tmpWithinExonVec_start_doner.size();
		exonWithinBackSplice_ofs << "Doner site: " << endl;
		for(int tmp = 0; tmp < tmpWithinExonNum_doner; tmp++)
			exonWithinBackSplice_ofs << tmpWithinExonVec_start_doner[tmp] << "\t" << tmpWithinExonVec_end_doner[tmp] << "\t" << tmpWithinExonVec_transcriptId_doner[tmp] << endl;
		int tmpWithinExonNum_acceptor = tmpWithinExonVec_start_acceptor.size();
		exonWithinBackSplice_ofs << "Acceptor site: " << endl;
		for(int tmp = 0; tmp < tmpWithinExonNum_acceptor; tmp++)
			exonWithinBackSplice_ofs << tmpWithinExonVec_start_acceptor[tmp] << "\t" << tmpWithinExonVec_end_acceptor[tmp] << "\t" << tmpWithinExonVec_transcriptId_acceptor[tmp] << endl;
		//exonWithinBackSplice_seq_ofs << "--------------------------------------------------------" << endl;
		exonWithinBackSplice_seq_ofs << tmpStr;// << endl;
		//if(tmpWithinExonNum_doner > 0)
		//{
			exonWithinBackSplice_seq_ofs << "\t" << tmpBackSplice_doner_startPos << "\t" << tmpBackSplice_doner_endPos;
			string tmpSeq_doner = indexInfo->returnChromStrSubstr(tmpBackSplice_chrNameInt, 
				tmpBackSplice_doner_startPos, seqLength);
			exonWithinBackSplice_seq_ofs << "\t" << tmpSeq_doner << "\t" << getRcmSeq(tmpSeq_doner);
		//}
		//if(tmpWithinExonNum_acceptor > 0)
		//{
			exonWithinBackSplice_seq_ofs << "\t" << tmpBackSplice_acceptor_startPos << "\t" << tmpBackSplice_acceptor_endPos;
			string tmpSeq_acceptor = indexInfo->returnChromStrSubstr(tmpBackSplice_chrNameInt,
				tmpBackSplice_acceptor_startPos, seqLength);
			exonWithinBackSplice_seq_ofs << "\t" << tmpSeq_acceptor << "\t" << getRcmSeq(tmpSeq_acceptor);
		//}
		exonWithinBackSplice_seq_ofs << endl;
		if((tmpWithinExonNum_doner == 0)||(tmpWithinExonNum_acceptor == 0))
		{
			cout << "----------------------------------------------------" << endl;
			cout << "error !" << endl;
			cout << "tmpWithinExonNum_doner: " << tmpWithinExonNum_doner << endl;
			cout << "tmpWithinExonNum_acceptor: " << tmpWithinExonNum_acceptor << endl;
			cout << tmpBackSplice_chrName << "\t" << tmpBackSplice_doner_startPos << "\t" << tmpBackSplice_doner_endPos << endl;
			exonWithinBackSplice_log_ofs << "----------------------------------------------------" << endl;
			exonWithinBackSplice_log_ofs << "error !" << endl;
			exonWithinBackSplice_log_ofs << "tmpWithinExonNum_doner: " << tmpWithinExonNum_doner << endl;
			exonWithinBackSplice_log_ofs << "tmpWithinExonNum_acceptor: " << tmpWithinExonNum_acceptor << endl;
			exonWithinBackSplice_log_ofs << tmpBackSplice_chrName << "\t" << tmpBackSplice_doner_startPos << "\t" << tmpBackSplice_doner_endPos << endl;

		}
	}
	backSplice_ifs.close();
	exonWithinBackSplice_ofs.close();
	exonGTF_ifs.close();

	return 0;
}