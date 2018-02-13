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

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputGTF outputFolder" << endl;
		exit(1);
	}
	string inputGTF = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	string annotatedSJ_outputFile = outputFolderStr + "/annotatedSJ.txt";
	string novelSJ_outputFile = outputFolderStr + "/novelSJ_bothEndsAnnotated.txt";
	ofstream annSJ_ofs(annotatedSJ_outputFile.c_str());
	ofstream novelSJ_ofs(novelSJ_outputFile.c_str());

	cout << "start to read GTF file" << endl;
	log_ofs << "start to read GTF file" << endl;
	vector< vector<string> > geneGTFstrVecVec;
	vector<string> geneNameVec;
	vector<string> geneStrandVec;
	vector<string> geneChrNameVec;
	ifstream gtf_ifs(inputGTF.c_str());
	string tmpCurrentGeneName = "NULL_GENE";
	int tmpLineNO = 0;
	while(!gtf_ifs.eof())
	{
		string tmpStr;
		getline(gtf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		tmpLineNO ++;
		int tmpThousandIndex = tmpLineNO / 10000;
		if(tmpLineNO == tmpThousandIndex * 10000)
			cout << "Processed Line #: " << tmpLineNO << endl;			
		int tmp_1stTab_loc = tmpStr.find("\t");
		int tmp_2ndTab_loc = tmpStr.find("\t", tmp_1stTab_loc + 1);
		string tmpFeatureStr = tmpStr.substr(tmp_1stTab_loc + 1, tmp_2ndTab_loc - tmp_1stTab_loc - 1);
		int tmp_3rdTab_loc = tmpStr.find("\t", tmp_2ndTab_loc + 1);
		string tmpTypeStr = tmpStr.substr(tmp_2ndTab_loc + 1, tmp_3rdTab_loc - tmp_2ndTab_loc - 1);
		int tmp_gene_name_loc = tmpStr.find("gene_name");
		int tmp_next_1stQuote_loc = tmpStr.find("\"", tmp_gene_name_loc + 1);
		int tmp_next_2ndQuote_loc = tmpStr.find("\"", tmp_next_1stQuote_loc + 1);
		if(
			(tmpFeatureStr != "protein_coding")||
			(tmpTypeStr != "exon")||(tmp_next_1stQuote_loc == string::npos)||(tmp_next_2ndQuote_loc == string::npos))
			continue;
		string tmp_gene_name_str = tmpStr.substr(tmp_next_1stQuote_loc + 1, tmp_next_2ndQuote_loc - tmp_next_1stQuote_loc - 1);
		if(tmp_gene_name_str != tmpCurrentGeneName)
		{
			vector<string> tmpGeneGtfStrVec;
			tmpGeneGtfStrVec.push_back(tmpStr);
			geneGTFstrVecVec.push_back(tmpGeneGtfStrVec);
			geneNameVec.push_back(tmp_gene_name_str);
			int tmp_4thTab_loc = tmpStr.find("\t", tmp_3rdTab_loc + 1);
			int tmp_5thTab_loc = tmpStr.find("\t", tmp_4thTab_loc + 1);
			int tmp_6thTab_loc = tmpStr.find("\t", tmp_5thTab_loc + 1);
			int tmp_7thTab_loc = tmpStr.find("\t", tmp_6thTab_loc + 1);					
			string tmp_strand_str = tmpStr.substr(tmp_6thTab_loc + 1, tmp_7thTab_loc - tmp_6thTab_loc - 1);
			if((tmp_strand_str != "+")&&(tmp_strand_str != "-"))
			{
				cout << "error in tmp_strand_str: " << tmp_strand_str << endl;
				exit(1);
			}
			geneStrandVec.push_back(tmp_strand_str);
			string tmp_chrName_str = tmpStr.substr(0, tmp_1stTab_loc);
			geneChrNameVec.push_back(tmp_chrName_str);
			tmpCurrentGeneName = tmp_gene_name_str;
			//cout << "new gene !" << endl;
			//cout << "tmpStr: " << tmpStr << endl;
		}
		else
		{
			//cout << "existing gene !" << endl;
			//cout << "tmpStr: " << tmpStr << endl;
			geneGTFstrVecVec[geneGTFstrVecVec.size() - 1].push_back(tmpStr);
		}
	}
	gtf_ifs.close();
	cout << "geneNameVec.size(): " << geneNameVec.size() << endl;
	log_ofs << "geneNameVec.size(): " << geneNameVec.size() << endl;
	int gene_transcriptNumSum = 0;
	int gene_annotatedSJnumSum = 0;
	int gene_novelSJnumSum= 0;
	for(int tmp = 0; tmp < geneNameVec.size(); tmp++)
	{
		//cout << "tmpGene_NO: " << tmp + 1 << endl;
		string tmpGene_name = geneNameVec[tmp];
		string tmpGene_strand = geneStrandVec[tmp];
		string tmpGene_chrName = geneChrNameVec[tmp];
		//cout << "tmpGene_name: " << tmpGene_name << endl;
		// cout << "tmpGene_strand: " << tmpGene_strand << endl;
		//cout << "tmpGene_chrName: " << tmpGene_chrName << endl;
		int tmpGeneGtfLineNum = geneGTFstrVecVec[tmp].size();
		//cout << "tmpGeneGtfLineNum: " << tmpGeneGtfLineNum << endl;
		vector< vector<string> > tmpGeneTranscriptStrVecVec;
		string tmpCurrentTranscriptId = "NULL_transcript";
		for(int tmpLine = 0; tmpLine < tmpGeneGtfLineNum; tmpLine ++)
		{
			string tmpTranscriptGtfStr = (geneGTFstrVecVec[tmp])[tmpLine];
			//cout << "tmpTranscriptGtfStr: " << tmpTranscriptGtfStr << endl;
			int tmp_transcript_id_loc = tmpTranscriptGtfStr.find("transcript_id");
			int tmp_next_1stQuote_loc = tmpTranscriptGtfStr.find("\"", tmp_transcript_id_loc + 1);
			int tmp_next_2ndQuote_loc = tmpTranscriptGtfStr.find("\"", tmp_next_1stQuote_loc + 1);
			string tmp_transcript_id_str = tmpTranscriptGtfStr.substr(tmp_next_1stQuote_loc + 1, tmp_next_2ndQuote_loc - tmp_next_1stQuote_loc - 1);
			if(tmp_transcript_id_str != tmpCurrentTranscriptId)
			{
				vector<string> tmpTranscriptGtfStrVec;
				tmpTranscriptGtfStrVec.push_back(tmpTranscriptGtfStr);
				tmpGeneTranscriptStrVecVec.push_back(tmpTranscriptGtfStrVec);
				tmpCurrentTranscriptId = tmp_transcript_id_str;
			}
			else
				tmpGeneTranscriptStrVecVec[tmpGeneTranscriptStrVecVec.size() - 1].push_back(tmpTranscriptGtfStr);
		}
		//cout << "tmpTranscriptNum: " << tmpGeneTranscriptStrVecVec.size() << endl;
		vector< vector< pair<int,int> > > tmpGeneTranscriptExonVecVec;
		int tmpGene_transcriptNum = tmpGeneTranscriptStrVecVec.size();
		gene_transcriptNumSum += tmpGene_transcriptNum;
		for(int tmpTranscript = 0; tmpTranscript < tmpGene_transcriptNum; tmpTranscript ++)
		{
			int tmpTranscript_exonNum = tmpGeneTranscriptStrVecVec[tmpTranscript].size();
			//gene_exonNum += tmpTranscript_exonNum;
			//cout << "tmpExonNum: " << tmpTranscript_exonNum << endl;
			vector< pair<int,int> > tmpExonPairVec;
			for(int tmpExon = 0; tmpExon < tmpTranscript_exonNum; tmpExon ++)
			{
				string tmpExonStr = (tmpGeneTranscriptStrVecVec[tmpTranscript])[tmpExon];
				//cout << "tmpExonStr: " << tmpExonStr << endl;
				int tmp_1stTab_loc = tmpExonStr.find("\t");
				int tmp_2ndTab_loc = tmpExonStr.find("\t", tmp_1stTab_loc + 1);
				int tmp_3rdTab_loc = tmpExonStr.find("\t", tmp_2ndTab_loc + 1);
				int tmp_4thTab_loc = tmpExonStr.find("\t", tmp_3rdTab_loc + 1);
				int tmp_5thTab_loc = tmpExonStr.find("\t", tmp_4thTab_loc + 1);				
				string tmpExon_startPosStr = tmpExonStr.substr(tmp_3rdTab_loc + 1, tmp_4thTab_loc - tmp_3rdTab_loc - 1);
				string tmpExon_endPosStr = tmpExonStr.substr(tmp_4thTab_loc + 1, tmp_5thTab_loc - tmp_4thTab_loc - 1);
				int tmpExon_startPos = atoi(tmpExon_startPosStr.c_str());
				int tmpExon_endPos = atoi(tmpExon_endPosStr.c_str());
				//cout << "tmpExon_startPos: " << tmpExon_startPos << endl;
				//cout << "tmpExon_endPos: " << tmpExon_endPos << endl;
				tmpExonPairVec.push_back(pair<int,int> (tmpExon_startPos, tmpExon_endPos));
			}
			tmpGeneTranscriptExonVecVec.push_back(tmpExonPairVec);
		}
		//cout << "tmpGene_exonNum: " << tmpGene_exonNum << endl;
		vector< pair<int,int> > tmpSpliceSitePosPairVec;
		set<int> tmpSpliceStartPosSet;
		set<int> tmpSpliceEndPosSet;
		for(int tmpTranscript = 0; tmpTranscript < tmpGene_transcriptNum; tmpTranscript ++)
		{
			int tmpTranscript_exonNum = tmpGeneTranscriptExonVecVec[tmpTranscript].size();
			if(tmpTranscript_exonNum <= 1)
				continue;
			for(int tmpExon = 0; tmpExon < tmpTranscript_exonNum - 1; tmpExon ++)
			{
				int tmpExonIndex_1st, tmpExonIndex_2nd;
				if(tmpGene_strand == "+")
				{
					tmpExonIndex_1st = tmpExon;
					tmpExonIndex_2nd = tmpExon + 1;					
				}
				else if(tmpGene_strand == "-")
				{
					tmpExonIndex_1st = tmpExon + 1;
					tmpExonIndex_2nd = tmpExon;
				}
				else
				{
					cout << "error in strand: " << tmpGene_strand << endl;
					exit(1);					
				}
				int tmpSplice_startPos = (tmpGeneTranscriptExonVecVec[tmpTranscript])[tmpExonIndex_1st].second;
				int tmpSplice_endPos = (tmpGeneTranscriptExonVecVec[tmpTranscript])[tmpExonIndex_2nd].first;								
				tmpSpliceStartPosSet.insert(tmpSplice_startPos);
				tmpSpliceEndPosSet.insert(tmpSplice_endPos);
				int tmpExistingSJnum = tmpSpliceSitePosPairVec.size();
				bool tmpSJ_novel_bool = true;
				for(int tmpSJ = 0; tmpSJ < tmpExistingSJnum; tmpSJ ++)
				{
					int tmpExistingSJ_startPos = tmpSpliceSitePosPairVec[tmpSJ].first;
					int tmpExistingSJ_endPos = tmpSpliceSitePosPairVec[tmpSJ].second;
					if((tmpSplice_startPos == tmpExistingSJ_startPos)&&(tmpSplice_endPos == tmpExistingSJ_endPos))
						tmpSJ_novel_bool = false;
				}
				if(tmpSJ_novel_bool)
					tmpSpliceSitePosPairVec.push_back(pair<int,int>(tmpSplice_startPos, tmpSplice_endPos));
			}
		}
		vector< pair<int,int> > tmpSpliceSitePosPairVec_novelNotInGeneAnn;
		set<int>::iterator tmpIntSetIter_startPos;
		set<int>::iterator tmpIntSetIter_endPos;
		for(tmpIntSetIter_startPos = tmpSpliceStartPosSet.begin(); tmpIntSetIter_startPos != tmpSpliceStartPosSet.end();
			tmpIntSetIter_startPos ++)
		{
			int tmpSpliceSite_startPos = (*tmpIntSetIter_startPos);
			for(tmpIntSetIter_endPos = tmpSpliceEndPosSet.begin(); tmpIntSetIter_endPos != tmpSpliceEndPosSet.end();
				tmpIntSetIter_endPos ++)
			{
				int tmpSpliceSite_endPos = (*tmpIntSetIter_endPos);
				//cout << "tmpSJ: " << tmpSpliceSite_startPos << " " << tmpSpliceSite_endPos << endl;
				if(tmpSpliceSite_startPos >= tmpSpliceSite_endPos)
					continue;
				int tmpExistingSJnum = tmpSpliceSitePosPairVec.size();
				bool tmpSJ_novel_bool = true;
				for(int tmpSJ = 0; tmpSJ < tmpExistingSJnum; tmpSJ ++)
				{
					int tmpExistingSJ_startPos = tmpSpliceSitePosPairVec[tmpSJ].first;
					int tmpExistingSJ_endPos = tmpSpliceSitePosPairVec[tmpSJ].second;
					if((tmpSpliceSite_startPos == tmpExistingSJ_startPos)&&(tmpSpliceSite_endPos == tmpExistingSJ_endPos))
						tmpSJ_novel_bool = false;
				}
				//cout << "tmpSJ_novel_bool: " << tmpSJ_novel_bool << endl;
				if(tmpSJ_novel_bool)
					tmpSpliceSitePosPairVec_novelNotInGeneAnn.push_back(pair<int,int>(
						tmpSpliceSite_startPos, tmpSpliceSite_endPos));
			}
		}

		int tmpSJnum_annotated = tmpSpliceSitePosPairVec.size();
		gene_annotatedSJnumSum += tmpSJnum_annotated;
		int tmpSJnum_novel = tmpSpliceSitePosPairVec_novelNotInGeneAnn.size();
		gene_novelSJnumSum += tmpSJnum_novel;
		//cout << "tmpSJnum_annotated: " << tmpSJnum_annotated << endl;
		//cout << "tmpSJnum_novel: " << tmpSJnum_novel << endl;
		for(int tmp = 0; tmp < tmpSJnum_annotated; tmp++)
			annSJ_ofs << tmpGene_chrName << "\t" << tmpSpliceSitePosPairVec[tmp].first << "\t" 
				<< tmpSpliceSitePosPairVec[tmp].second << "\t" << tmpGene_strand << "\t" << tmpGene_name << endl;
		for(int tmp = 0; tmp < tmpSJnum_novel; tmp++)
			novelSJ_ofs << tmpGene_chrName << "\t" << tmpSpliceSitePosPairVec_novelNotInGeneAnn[tmp].first << "\t" 
				<< tmpSpliceSitePosPairVec_novelNotInGeneAnn[tmp].second << "\t" << tmpGene_strand << "\t" << tmpGene_name << endl;
	}
	log_ofs << "transcriptNumSum: " << gene_transcriptNumSum << endl;
	log_ofs << "annotatedSJnumSum: " << gene_annotatedSJnumSum << endl;
	log_ofs << "novelSJnumSum: " << gene_novelSJnumSum << endl;
	log_ofs.close();
	annSJ_ofs.close();
	novelSJ_ofs.close();
	return 0;
}