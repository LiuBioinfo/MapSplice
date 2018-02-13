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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

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
		cout << "Executable inputIndexFolder inputExonGTFfile inputBackSpliceJuncFile outputFilePrefix" << endl;
		exit(1);
	}
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
		//cout << "tmpExon_chrNameInt: " << tmpExon_chrNameInt << endl;
		//cout << "tmpExon_startPos: " << tmpExon_startPos << endl;
		//cout << "tmpExon_endPos: " << tmpExon_endPos << endl;
		exonPosPairVecVec[tmpExon_chrNameInt].push_back(pair<int,int>(tmpExon_startPos, tmpExon_endPos));
		exonTranscriptIdVecVec[tmpExon_chrNameInt].push_back(tmpExon_transcriptId);
	}
	cout << "start to search for exons within back splice junctions" << endl;
	string outputFilePrefix = argv[4];
	string outputExonWithinBackSplice = outputFilePrefix + "_exon.txt";
	ofstream exonWithinBackSplice_ofs(outputExonWithinBackSplice.c_str());
	string outputTranscriptWithinBackSplice = outputFilePrefix + "_transcript.txt";
	ofstream transcriptWithinBackSplice_ofs(outputTranscriptWithinBackSplice.c_str());
	string outputTranscriptWithinBackSplice_grouped = outputFilePrefix + "_transcript_grouped.txt";
	ofstream transcriptWithinBackSplice_grouped_ofs(outputTranscriptWithinBackSplice_grouped.c_str());	
	string outputTranscriptWithinBackSplice_grouped_nonDistance = outputFilePrefix + "_transcript_grouped_nonDistance.txt";
	ofstream transcriptWithinBackSplice_grouped_nonDistance_ofs(outputTranscriptWithinBackSplice_grouped_nonDistance.c_str());	

	string inputBackSpliceJuncFile = argv[3];
	ifstream backSplice_ifs(inputBackSpliceJuncFile.c_str());
	while(!backSplice_ifs.eof())
	{
		string tmpStr;
		getline(backSplice_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpBackSplice_chrName = tmpStr.substr(0, tabLoc_1);
		string tmpBackSplice_donerStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpBackSplice_acceptorStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpBackSplice_chrNameInt = indexInfo->convertStringToInt(tmpBackSplice_chrName);
		int tmpBackSplice_doner = atoi(tmpBackSplice_donerStr.c_str());
		int tmpBackSplice_acceptor = atoi(tmpBackSplice_acceptorStr.c_str());
		int tmpRegion_start = tmpBackSplice_acceptor;
		int tmpRegion_end = tmpBackSplice_doner;
		vector<int> tmpWithinExonVec_start;
		vector<int> tmpWithinExonVec_end;
		vector<string> tmpWithinExonVec_transcriptId;
		int tmpChrExonNum = exonPosPairVecVec[tmpBackSplice_chrNameInt].size();
		for(int tmp = 0; tmp < tmpChrExonNum; tmp++)
		{
			int tmpChrExon_start = (exonPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].first;
			int tmpChrExon_end = (exonPosPairVecVec[tmpBackSplice_chrNameInt])[tmp].second;
			if((tmpChrExon_start >= tmpRegion_start)&&(tmpChrExon_end <= tmpRegion_end))
			{
				tmpWithinExonVec_start.push_back(tmpChrExon_start);
				tmpWithinExonVec_end.push_back(tmpChrExon_end);
				tmpWithinExonVec_transcriptId.push_back((exonTranscriptIdVecVec[tmpBackSplice_chrNameInt])[tmp]);
			}
		}
		exonWithinBackSplice_ofs << "--------------------------------------------------------" << endl;
		exonWithinBackSplice_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr << "--" << tmpBackSplice_acceptorStr << ":" << endl;
		int tmpWithinExonNum = tmpWithinExonVec_start.size();
		for(int tmp = 0; tmp < tmpWithinExonNum; tmp++)
			exonWithinBackSplice_ofs << tmpWithinExonVec_start[tmp] << "\t" << tmpWithinExonVec_end[tmp] << "\t"
				<< tmpWithinExonVec_transcriptId[tmp] << endl;
		
		vector<string> tmpWithinTranscript_idVec;
		for(int tmp = 0; tmp < tmpWithinExonNum; tmp++)
		{
			string tmpTranscriptId = tmpWithinExonVec_transcriptId[tmp];
			int currentTranscriptNum = tmpWithinTranscript_idVec.size();
			bool tmpTranscriptExists = false;
			for(int tmp2 = 0; tmp2 < currentTranscriptNum; tmp2++)
			{
				if(tmpTranscriptId == tmpWithinTranscript_idVec[tmp2])
				{
					tmpTranscriptExists = true;
					break;
				}
			}
			if(!tmpTranscriptExists)
				tmpWithinTranscript_idVec.push_back(tmpTranscriptId);
		}

		vector< vector< pair<int,int> > > tmpWithinTranscript_exonPosPairVecVec_raw;
		int tmpWithinTranscript_idVec_size = tmpWithinTranscript_idVec.size();
		for(int tmpTran = 0; tmpTran < tmpWithinTranscript_idVec_size; tmpTran ++)
		{
			vector< pair<int,int> > tmpWithinTranscript_exonPosPairVec_raw;
			tmpWithinTranscript_exonPosPairVecVec_raw.push_back(tmpWithinTranscript_exonPosPairVec_raw);
		}

		for(int tmp = 0; tmp < tmpWithinExonNum; tmp++)
		{
			string tmpTranscriptId = tmpWithinExonVec_transcriptId[tmp];
			for(int tmpTran = 0; tmpTran < tmpWithinTranscript_idVec_size; tmpTran ++)
			{
				if(tmpTranscriptId == tmpWithinTranscript_idVec[tmpTran])
				{
					int tmpExon_startPos = tmpWithinExonVec_start[tmp];
					int tmpExon_endPos = tmpWithinExonVec_end[tmp];
					tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran].push_back(
						pair<int,int>(tmpExon_startPos, tmpExon_endPos));
					break;
				}
			}
		}

		vector< vector<int> > tmpTranscriptGroupedVec;
		for(int tmpTran = 0; tmpTran < tmpWithinTranscript_idVec_size; tmpTran ++)
		{
			int currentGroupedVecSize = tmpTranscriptGroupedVec.size();
			bool theSameWithExistingGroup = false;
			for(int tmpGroup = 0; tmpGroup < currentGroupedVecSize; tmpGroup ++)
			{
				int tmpTranscriptIndex_in_tmpWithinTranscript_idVec = (tmpTranscriptGroupedVec[tmpGroup])[0];
				bool twoTranscriptSegTheSameBool = twoTranscriptSegTheSameOrNot(
					tmpWithinTranscript_exonPosPairVecVec_raw[tmpTranscriptIndex_in_tmpWithinTranscript_idVec],
					tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran]);
				if(twoTranscriptSegTheSameBool)
				{	
					theSameWithExistingGroup = true;
					tmpTranscriptGroupedVec[tmpGroup].push_back(tmpTran);
					break;
				}
			}
			if(!theSameWithExistingGroup)
			{
				vector<int> tmpGroup;
				tmpGroup.push_back(tmpTran);
				tmpTranscriptGroupedVec.push_back(tmpGroup);
			}
		}
		transcriptWithinBackSplice_grouped_ofs << "--------------------------------------------------------" << endl;
		transcriptWithinBackSplice_grouped_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr 
			<< "--" << tmpBackSplice_acceptorStr << ":" << endl;
		transcriptWithinBackSplice_grouped_nonDistance_ofs << "--------------------------------------------------------" << endl;
		transcriptWithinBackSplice_grouped_nonDistance_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr 
			<< "--" << tmpBackSplice_acceptorStr << ":" << endl;

		for(int tmpGroup = 0; tmpGroup < tmpTranscriptGroupedVec.size(); tmpGroup ++)
		{
			int tmpTranscriptNumInGroup = (tmpTranscriptGroupedVec[tmpGroup]).size();
			for(int tmpTranInGroup = 0; tmpTranInGroup < tmpTranscriptNumInGroup; tmpTranInGroup ++)
			{
				int tmpTran = (tmpTranscriptGroupedVec[tmpGroup])[tmpTranInGroup];
				string tmpTran_id = tmpWithinTranscript_idVec[tmpTran];
				transcriptWithinBackSplice_grouped_ofs << tmpTran_id << ";";
			}
			int tmpTran_1st = (tmpTranscriptGroupedVec[tmpGroup])[0];
			int tmpTran_1st_exonNum = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran_1st]).size();
			vector<int> tmpTran_1st_exonStartPosVec;
			vector<int> tmpTran_1st_exonEndPosVec;
			for(int tmpExon = 0; tmpExon < tmpTran_1st_exonNum; tmpExon ++)
			{
				int tmpExon_startPos = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran_1st])[tmpExon].first;
				int tmpExon_endPos = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran_1st])[tmpExon].second;
				tmpTran_1st_exonStartPosVec.push_back(tmpExon_startPos);
				tmpTran_1st_exonEndPosVec.push_back(tmpExon_endPos);
			}
			sort(tmpTran_1st_exonStartPosVec.begin(), tmpTran_1st_exonStartPosVec.end(), sortExonPosOrderBool);
			sort(tmpTran_1st_exonEndPosVec.begin(), tmpTran_1st_exonEndPosVec.end(), sortExonPosOrderBool);			
			int tmpLeftDistance = tmpTran_1st_exonStartPosVec[0] - tmpBackSplice_acceptor;
			int tmpRightDistance = tmpBackSplice_doner - tmpTran_1st_exonEndPosVec[tmpTran_1st_exonNum - 1];
			transcriptWithinBackSplice_grouped_ofs << " LEFT_DIST:" << tmpLeftDistance << " RIGHT_DIST:" 
				<< tmpRightDistance << endl;
			transcriptWithinBackSplice_grouped_ofs << tmpTran_1st_exonStartPosVec[0] << ":" << tmpTran_1st_exonEndPosVec[0];
			if(tmpTran_1st_exonNum > 1)
			{
				for(int tmpExon = 1; tmpExon < tmpTran_1st_exonNum; tmpExon ++)
					transcriptWithinBackSplice_grouped_ofs << "-" << tmpTran_1st_exonStartPosVec[tmpExon] 
						<< ":" << tmpTran_1st_exonEndPosVec[tmpExon];
			}
			transcriptWithinBackSplice_grouped_ofs << endl;
		
			if((tmpLeftDistance == 0)&&(tmpRightDistance == 0))
			{
				for(int tmpTranInGroup = 0; tmpTranInGroup < tmpTranscriptNumInGroup; tmpTranInGroup ++)
				{
					int tmpTran = (tmpTranscriptGroupedVec[tmpGroup])[tmpTranInGroup];
					string tmpTran_id = tmpWithinTranscript_idVec[tmpTran];
					transcriptWithinBackSplice_grouped_nonDistance_ofs << tmpTran_id << ";";
				}				
				transcriptWithinBackSplice_grouped_nonDistance_ofs << " LEFT_DIST:0 RIGHT_DIST:0" << endl;
				transcriptWithinBackSplice_grouped_nonDistance_ofs << tmpTran_1st_exonStartPosVec[0] << ":" << tmpTran_1st_exonEndPosVec[0];
				if(tmpTran_1st_exonNum > 1)
				{
					for(int tmpExon = 1; tmpExon < tmpTran_1st_exonNum; tmpExon ++)
						transcriptWithinBackSplice_grouped_nonDistance_ofs << "-" << tmpTran_1st_exonStartPosVec[tmpExon] 
							<< ":" << tmpTran_1st_exonEndPosVec[tmpExon];
				}
				transcriptWithinBackSplice_grouped_nonDistance_ofs << endl;
			}
		}

		transcriptWithinBackSplice_ofs << "--------------------------------------------------------" << endl;
		transcriptWithinBackSplice_ofs << tmpBackSplice_chrName << ":" << tmpBackSplice_donerStr 
			<< "--" << tmpBackSplice_acceptorStr << ":" << endl;
		for(int tmpTran = 0; tmpTran < tmpWithinTranscript_idVec_size; tmpTran ++)
		{
			int tmpTran_exonNum = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran]).size();
			string tmpTran_id = tmpWithinTranscript_idVec[tmpTran];
			vector<int> tmpTran_exonStartPosVec;
			vector<int> tmpTran_exonEndPosVec;
			for(int tmpExon = 0; tmpExon < tmpTran_exonNum; tmpExon ++)
			{
				int tmpExon_startPos = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran])[tmpExon].first;
				int tmpExon_endPos = (tmpWithinTranscript_exonPosPairVecVec_raw[tmpTran])[tmpExon].second;
				tmpTran_exonStartPosVec.push_back(tmpExon_startPos);
				tmpTran_exonEndPosVec.push_back(tmpExon_endPos);
			}
			sort(tmpTran_exonStartPosVec.begin(), tmpTran_exonStartPosVec.end(), sortExonPosOrderBool);
			sort(tmpTran_exonEndPosVec.begin(), tmpTran_exonEndPosVec.end(), sortExonPosOrderBool);

			int tmpLeftDistance = tmpTran_exonStartPosVec[0] - tmpBackSplice_acceptor;
			int tmpRightDistance = tmpBackSplice_doner - tmpTran_exonEndPosVec[tmpTran_exonNum - 1];

			transcriptWithinBackSplice_ofs << tmpTran_id << " LEFT_DIST:" << tmpLeftDistance << " RIGHT_DIST:" << tmpRightDistance << endl;
			transcriptWithinBackSplice_ofs << tmpTran_exonStartPosVec[0] << ":" << tmpTran_exonEndPosVec[0];
			if(tmpTran_exonNum > 1)
			{
				for(int tmpExon = 1; tmpExon < tmpTran_exonNum; tmpExon ++)
					transcriptWithinBackSplice_ofs << "-" << tmpTran_exonStartPosVec[tmpExon] << ":" << tmpTran_exonEndPosVec[tmpExon];
			}
			transcriptWithinBackSplice_ofs << endl;
		}

		int withinExonNum = tmpWithinExonVec_start.size();
		if(withinExonNum >= 1)
		{
			int tmpLeftMost = returnLeftMostExonStartSite(tmpWithinExonVec_start);
			int tmpRightMost = returnRightMostExonEndSite(tmpWithinExonVec_end);
			int tmpLeftDistance = tmpLeftMost - tmpBackSplice_acceptor;
			int tmpRightDistance = tmpBackSplice_doner - tmpRightMost;
			exonWithinBackSplice_ofs << "leftDistance:\t" << tmpLeftDistance << "\trightDistance:\t" << tmpRightDistance << endl;
		}
	}
	backSplice_ifs.close();
	exonWithinBackSplice_ofs.close();
	exonGTF_ifs.close();

	return 0;
}