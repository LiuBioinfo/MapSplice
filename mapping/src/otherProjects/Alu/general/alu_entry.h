// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALU_ENTRY_H
#define ALU_ENTRY_H

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
#include "../../../general/index_info.h"

using namespace std;

class Alu_Entry
{
private:
	int chrNameInt;
	int startPos;
	int endPos;
	string strand;
	string geneId;
	string transcriptId;

public:

	void initiate(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, 
		string& tmpStrand, string& tmpGeneId, string& tmpTranscriptId)
	{
		chrNameInt = tmpChrNameInt;
		startPos = tmpStartPos;
		endPos = tmpEndPos;
		strand = tmpStrand;
		geneId = tmpGeneId;
		transcriptId = tmpTranscriptId;
	}	

	void initiate(string tmpChrName, int tmpStartPos, int tmpEndPos, 
		string& tmpStrand, string& tmpGeneId, string& tmpTranscriptId, Index_Info* indexInfo)
	{
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		this->initiate(tmpChrNameInt, tmpStartPos, tmpEndPos, tmpStrand, tmpGeneId, tmpTranscriptId);
	}
	
	string returnAnnEntry(Index_Info* indexInfo)
	{
		string chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string tmpStr;
		tmpStr = chrNameStr + "\t" + int_to_str(startPos) + "\t" + int_to_str(endPos)
			+ "\t" + strand + "\t" + geneId + "\t" + transcriptId;
		return tmpStr;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnStartPos()
	{
		return startPos;
	}

	int returnEndPos()
	{
		return endPos;
	}

	string returnStrand()
	{
		return strand;
	}

	string returnGeneId()
	{
		return geneId;
	}

	string returnTranscriptId()
	{
		return transcriptId;
	}
};

#endif