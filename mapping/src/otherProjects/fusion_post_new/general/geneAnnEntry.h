#ifndef GENEANNENTRY_H
#define GENEANNENTRY_H

#include <stdio.h>
#include <stdlib.h>
//#include "../../../general/index_info.h"
#include "../../../general/index_info.h"

using namespace std;

class GeneAnnEntry_Info
{
private:
	int chrNameInt;
	int startPos;
	int endPos;
	string strand;
	string source;
	string featureType;
	string geneName;
	string geneId;
	string transcriptName;
	string transcriptId;
public:
	GeneAnnEntry_Info()
	{}

	string returnAnnEntry(Index_Info* indexInfo)
	{
		string chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		string tmpStr;
		tmpStr = chrNameStr + "\t" + int_to_str(startPos) + "\t" + int_to_str(endPos)
			+ "\t" + strand + "\t" + featureType + "\t" + geneName + "\t" + geneId 
			+ "\t" + transcriptName + "\t" + transcriptId + "\t" + source;
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

	string returnSource()
	{
		return source;
	}

	string returnFeatureType()
	{
		return featureType;
	}

	string returnGeneName()
	{
		return geneName;
	}

	string returnGeneId()
	{
		return geneId;
	}

	string returnTranscriptId()
	{
		return transcriptId;
	}

	string returnTranscriptName()
	{
		return transcriptName;
	}

	void initiate(int tmpChrNameInt, int tmpStartPos, int tmpEndPos, string& tmpStrand, string& tmpSource, 
		string& tmpFeatureType, string& tmpGeneName, string& tmpGeneId, string& tmpTranscriptName, string& tmpTranscriptId)
	{
		chrNameInt = tmpChrNameInt;
		startPos = tmpStartPos;
		endPos = tmpEndPos;
		strand = tmpStrand;
		source = tmpSource;
		featureType = tmpFeatureType;
		geneName = tmpGeneName;
		geneId = tmpGeneId;
		transcriptName = tmpTranscriptName;
		transcriptId = tmpTranscriptId;
	}
};


#endif