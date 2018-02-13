#ifndef ExonBoundaryHash_info_H
#define	ExonBoundaryHash_info_H

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

class ExonBoundaryHash_Info
{
private:
	vector< set<int> > exonStartPosSetVec_posi;
	vector< set<int> > exonStartPosSetVec_nega;
	vector< set<int> > exonEndPosSetVec_posi;
	vector< set<int> > exonEndPosSetVec_nega;
public:
	ExonBoundaryHash_Info()
	{}

	void initiate_geneAnnEntryFile(Index_Info* indexInfo, string& geneAnnEntryFile)
	{
		int index_chrom_num = indexInfo->return_chromNum();
		for(int tmp = 0; tmp < index_chrom_num; tmp++)
		{
			set<int> tmpSet_exonStartPos_posi;
			set<int> tmpSet_exonStartPos_nega;
			set<int> tmpSet_exonEndPos_posi;
			set<int> tmpSet_exonEndPos_nega;
			exonStartPosSetVec_posi.push_back(tmpSet_exonStartPos_posi);
			exonStartPosSetVec_nega.push_back(tmpSet_exonStartPos_nega);
			exonEndPosSetVec_posi.push_back(tmpSet_exonEndPos_posi);
			exonEndPosSetVec_nega.push_back(tmpSet_exonEndPos_nega);					
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
			if(tmpType == "exon")
			{
				if(tmpStrand == "+")
				{
					exonStartPosSetVec_posi[tmpChrNameInt].insert(tmpStartPos);
					exonEndPosSetVec_posi[tmpChrNameInt].insert(tmpEndPos);
				}
				else // tmpStrand == "-"
				{
					exonStartPosSetVec_nega[tmpChrNameInt].insert(tmpStartPos);
					exonEndPosSetVec_nega[tmpChrNameInt].insert(tmpEndPos);
				}
			}
		}
		geneAnn_ifs.close();
	}

	bool confirm(bool search_startPos_or_endPos_bool, bool strand_posi_or_nega_bool,
		int tmpChrNameInt, int tmpPos)
	{
		if(search_startPos_or_endPos_bool) // start
		{
			if(strand_posi_or_nega_bool) // posi
			{
				if(exonStartPosSetVec_posi[tmpChrNameInt].find(tmpPos) 
					!= exonStartPosSetVec_posi[tmpChrNameInt].end())
					return true;
				else
					return false;
			}
			else // nega
			{
				if(exonStartPosSetVec_nega[tmpChrNameInt].find(tmpPos) 
					!= exonStartPosSetVec_nega[tmpChrNameInt].end())
					return true;
				else
					return false;
			}
		} 
		else // end
		{
			if(strand_posi_or_nega_bool) // posi
			{
				if(exonEndPosSetVec_posi[tmpChrNameInt].find(tmpPos) 
					!= exonEndPosSetVec_posi[tmpChrNameInt].end())
					return true;
				else
					return false;
			}
			else // nega
			{
				if(exonEndPosSetVec_nega[tmpChrNameInt].find(tmpPos) 
					!= exonEndPosSetVec_nega[tmpChrNameInt].end())
					return true;
				else
					return false;
			}
		} 
	}

	bool hangOverExonStartPosOrNot(Index_Info* indexInfo, int tmpChrNameInt,
		int buffer_startPos, int buffer_endPos)
	{
		for(int tmpPos = buffer_startPos; tmpPos <= buffer_endPos; tmpPos ++)
		{
			if((exonStartPosSetVec_posi[tmpChrNameInt].find(tmpPos) != exonStartPosSetVec_posi[tmpChrNameInt].end())
				||(exonStartPosSetVec_nega[tmpChrNameInt].find(tmpPos) != exonStartPosSetVec_nega[tmpChrNameInt].end()))
				return true;
		}
		return false;
	}

	bool hangOverExonEndPosOrNot(Index_Info* indexInfo, int tmpChrNameInt,
		int buffer_startPos, int buffer_endPos)
	{
		for(int tmpPos = buffer_startPos; tmpPos <= buffer_endPos; tmpPos ++)
		{
			if((exonEndPosSetVec_posi[tmpChrNameInt].find(tmpPos) != exonEndPosSetVec_posi[tmpChrNameInt].end())
				||(exonEndPosSetVec_nega[tmpChrNameInt].find(tmpPos) != exonEndPosSetVec_nega[tmpChrNameInt].end()))
				return true;
		}
		return false;
	}	
};
#endif