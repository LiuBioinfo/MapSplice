// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ExonBoundary2transcriptIdHash_info_H
#define	ExonBoundary2transcriptIdHash_info_H

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
#include <sstream>

using namespace std;

class ExonBoundary2transcriptIdHash_Info
{
private:
	vector< map<int, string> > exonStartPosMapVec_posi;
	vector< map<int, string> > exonStartPosMapVec_nega;
	vector< map<int, string> > exonEndPosMapVec_posi;
	vector< map<int, string> > exonEndPosMapVec_nega;
public:
	ExonBoundary2transcriptIdHash_Info()
	{}

	void initiate_geneAnnEntryFile(Index_Info* indexInfo, string& geneAnnEntryFile)
	{
		int index_chrom_num = indexInfo->return_chromNum();
		for(int tmp = 0; tmp < index_chrom_num; tmp++)
		{
			map<int, string> tmpMap_exonStartPos_posi;
			map<int, string> tmpMap_exonStartPos_nega;
			map<int, string> tmpMap_exonEndPos_posi;
			map<int, string> tmpMap_exonEndPos_nega;
			exonStartPosMapVec_posi.push_back(tmpMap_exonStartPos_posi);
			exonStartPosMapVec_nega.push_back(tmpMap_exonStartPos_nega);
			exonEndPosMapVec_posi.push_back(tmpMap_exonEndPos_posi);
			exonEndPosMapVec_nega.push_back(tmpMap_exonEndPos_nega);					
		}
		ifstream geneAnn_ifs(geneAnnEntryFile.c_str());
		while(!geneAnn_ifs.eof())
		{
			string tmpStr;
			getline(geneAnn_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc_1 = tmpStr.find("\t");
			int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
			int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
			int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
			int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
			int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
			int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
			int tabLoc_8 = tmpStr.find("\t", tabLoc_7 + 1);			
			string tmpChrName = tmpStr.substr(0, tabLoc_1);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
			if(tmpChrNameInt < 0)
				continue;
			string tmpStartPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
			string tmpEndPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
			string tmpStrand = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
			string tmpType = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
			int tmpStartPos = atoi(tmpStartPosStr.c_str());
			int tmpEndPos = atoi(tmpEndPosStr.c_str());
			string tmpTranscriptId = tmpStr.substr(tabLoc_7 + 1, tabLoc_8 - tabLoc_7 - 1);
			if(tmpType == "exon")
			{
				if(tmpStrand == "+")
				{
					map<int,string>::iterator tmpMapIter_startPos, tmpMapIter_endPos;
					tmpMapIter_startPos = exonStartPosMapVec_posi[tmpChrNameInt].find(tmpStartPos);
					if(tmpMapIter_startPos != exonStartPosMapVec_posi[tmpChrNameInt].end())
						(tmpMapIter_startPos->second) += (tmpTranscriptId + ",");
					else
						exonStartPosMapVec_posi[tmpChrNameInt].insert(pair<int,string>(tmpStartPos, (tmpTranscriptId + ",")));
					
					tmpMapIter_endPos = exonEndPosMapVec_posi[tmpChrNameInt].find(tmpEndPos);
					if(tmpMapIter_endPos != exonEndPosMapVec_posi[tmpChrNameInt].end())
						(tmpMapIter_endPos->second) += (tmpTranscriptId + ",");
					else
						exonEndPosMapVec_posi[tmpChrNameInt].insert(pair<int,string>(tmpEndPos, (tmpTranscriptId + ",")));
				}
				else // tmpStrand == "-"
				{
					map<int,string>::iterator tmpMapIter_startPos, tmpMapIter_endPos;
					tmpMapIter_startPos = exonStartPosMapVec_nega[tmpChrNameInt].find(tmpStartPos);
					if(tmpMapIter_startPos != exonStartPosMapVec_nega[tmpChrNameInt].end())
						(tmpMapIter_startPos->second) += (tmpTranscriptId + ",");
					else
						exonStartPosMapVec_nega[tmpChrNameInt].insert(pair<int,string>(tmpStartPos, (tmpTranscriptId + ",")));
					
					tmpMapIter_endPos = exonEndPosMapVec_nega[tmpChrNameInt].find(tmpEndPos);
					if(tmpMapIter_endPos != exonEndPosMapVec_nega[tmpChrNameInt].end())
						(tmpMapIter_endPos->second) += (tmpTranscriptId + ",");
					else					
						exonEndPosMapVec_nega[tmpChrNameInt].insert(pair<int,string>(tmpEndPos, (tmpTranscriptId + ",")));
				}
			}
		}
		geneAnn_ifs.close();
	}

	string searchAndReturnTranscriptId(bool search_exon_startPos_or_endPos_bool, 
		bool strand_posi_or_nega_bool, int tmpChrNameInt, int tmpPos)
	{
		if(search_exon_startPos_or_endPos_bool) // start
		{
			if(strand_posi_or_nega_bool) // posi
			{
				map<int, string>::iterator tmpMapIter = exonStartPosMapVec_posi[tmpChrNameInt].find(tmpPos);
				if(tmpMapIter != exonStartPosMapVec_posi[tmpChrNameInt].end())
					return (tmpMapIter->second);
				else
					return "NULL";
			}
			else // nega
			{
				map<int, string>::iterator tmpMapIter = exonStartPosMapVec_nega[tmpChrNameInt].find(tmpPos);
				if(tmpMapIter != exonStartPosMapVec_nega[tmpChrNameInt].end())
					return (tmpMapIter->second);
				else
					return "NULL";
			}
		} 
		else // end
		{
			if(strand_posi_or_nega_bool) // posi
			{
				map<int, string>::iterator tmpMapIter = exonEndPosMapVec_posi[tmpChrNameInt].find(tmpPos);
				if(tmpMapIter != exonEndPosMapVec_posi[tmpChrNameInt].end())
					return (tmpMapIter->second);
				else
					return "NULL";
			}
			else // nega
			{
				map<int, string>::iterator tmpMapIter = exonEndPosMapVec_nega[tmpChrNameInt].find(tmpPos);
				if(tmpMapIter != exonEndPosMapVec_nega[tmpChrNameInt].end())
					return (tmpMapIter->second);
				else
					return "NULL";
			}
		} 
	}

	void getFusionTranscriptIdPairFromBreakPoint(
		int chrNameInt_gene1, int chrNameInt_gene2,
		int pos_gene1, int pos_gene2, 
		bool strand_for_or_rev_gene1_bool, bool strand_for_or_rev_gene2_bool,
		string& transcriptId_1, string& transcriptId_2)
	{
		// gene1 
		if(strand_for_or_rev_gene1_bool) // +, exonEndPos
			transcriptId_1 = this->searchAndReturnTranscriptId(false, true, chrNameInt_gene1, pos_gene1);
		else // -, exonStartPos
			transcriptId_1 = this->searchAndReturnTranscriptId(true, false, chrNameInt_gene1, pos_gene1);
		// gene2
		if(strand_for_or_rev_gene2_bool) // +, exonStartPos
			transcriptId_2 = this->searchAndReturnTranscriptId(true, true, chrNameInt_gene2, pos_gene2);
		else // -, exonEndPos
			transcriptId_2 = this->searchAndReturnTranscriptId(false, false, chrNameInt_gene2, pos_gene2);
	}
};
#endif