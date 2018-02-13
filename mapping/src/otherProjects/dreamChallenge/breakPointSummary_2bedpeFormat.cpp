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
#include "../../otherProjects/indexing/utils/build_snpMerIndex.h"
#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../phase1/arrayQueue_phase1.h"
#include "../../phase2/arrayQueue_phase2.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "../../otherProjects/incorporateGenomicVariants/general/syntheticSNPtransSeq_info.h"
#include "../../otherProjects/incorporateGenomicVariants/general/learnedCandiSNPhash_info_vec.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"
#include "../../phase2/spliceJunctionHash_info.h"
#include "../../phase2/unmapEnd_info.h"
#include "../../phase2/unfixedHead.h"
#include "../../phase2/unfixedTail.h"
#include "../../phase2/incompleteLongHead.h"
#include "../../phase2/incompleteLongTail.h"
#include "../../phase2/sam2junc.h"
#include "../../fixHeadTail.h"
#include "../../phase2/fixOneEndUnmapped.h"
#include "../../fixPhase1.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"
#include "general/peSam_info.h"
#include "general/fixFusion_peRead.h"
#include "general/fixFusionResult_info.h"
#include "breakPointDetermination_dreamChallengeOnly.h"
#include "exonBoundary2transcriptIdHash_info.h"

using namespace std;

typedef map<int, pair<int,int> > BreakPoint2FusionInfoVecIndexMap;

time_t nowtime;
struct tm *local;

bool withinTheSameTranscriptOrNot(string& tmpTranscriptId_gene1, string& tmpTranscriptId_gene2)
{
	vector<string> tmpTranscriptIdVec_gene1;
	vector<string> tmpTranscriptIdVec_gene2;
	int commaLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmpNextCommaLoc = tmpTranscriptId_gene1.find(",", commaLoc);
		if(tmpNextCommaLoc == string::npos)
			break;
		tmpTranscriptIdVec_gene1.push_back(tmpTranscriptId_gene1.substr(commaLoc, tmpNextCommaLoc - commaLoc));
		commaLoc = tmpNextCommaLoc + 1;
	}
	commaLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tmpNextCommaLoc = tmpTranscriptId_gene2.find(",", commaLoc);
		if(tmpNextCommaLoc == string::npos)
			break;
		tmpTranscriptIdVec_gene2.push_back(tmpTranscriptId_gene2.substr(commaLoc, tmpNextCommaLoc - commaLoc));
		commaLoc = tmpNextCommaLoc + 1;
	}
	for(int tmp = 0; tmp < tmpTranscriptIdVec_gene1.size(); tmp++)
	{
		string tmp_tmpTranscriptId_gene1 = tmpTranscriptIdVec_gene1[tmp];
		for(int tmp2 = 0; tmp2 < tmpTranscriptIdVec_gene2.size(); tmp2 ++)
		{
			string tmp_tmpTranscriptId_gene2 = tmpTranscriptIdVec_gene2[tmp2];
			if(tmp_tmpTranscriptId_gene1 == tmp_tmpTranscriptId_gene2)
				return true;
		}
	}
	return false;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexDir" << endl;
		cout << "#2 geneAnn" << endl;
		cout << "#3 inputBreakPointFile_upstreamHead" << endl;
		cout << "#4 inputBreakPointFile_downstreamTail" << endl;
		cout << "#5 outputFilePrefix" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string geneAnnEntryFile = argv[2];
	string inputBreakPointFile_upstreamHead = argv[3];
	string inputBreakPointFile_downstreamTail = argv[4];
	string outputFilePrefix = argv[5];

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load indexInfo" << endl;
	cout << "initiate indexInfo ..." << endl;	
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	indexInfo->initiate_withoutLoadingSeq();
	parameter_ifs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading indexInfo" << endl;
	cout << endl << "[" << asctime(local) << "... start to load gene ann file" << endl;
	ExonBoundary2transcriptIdHash_Info exonBoundary2transcriptIdHashInfo;
	exonBoundary2transcriptIdHashInfo.initiate_geneAnnEntryFile(indexInfo, geneAnnEntryFile);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading gene ann file" << endl;
	cout << endl << "[" << asctime(local) << "... start to initiate breakPoint2FusionInfoVecIndexMapVecVec" << endl;
	vector< vector<BreakPoint2FusionInfoVecIndexMap> > breakPoint2FusionInfoVecIndexMapVecVec;
	for(int tmpChr_1 = 0; tmpChr_1 < chromNum; tmpChr_1 ++)
	{
		vector<BreakPoint2FusionInfoVecIndexMap> tmpBreakPoint2FusionInfoVecIndexMapVec;
		for(int tmpChr_2 = 0; tmpChr_2 < chromNum; tmpChr_2 ++)
		{
			BreakPoint2FusionInfoVecIndexMap tmpMap;
			tmpBreakPoint2FusionInfoVecIndexMapVec.push_back(tmpMap);
		}
		breakPoint2FusionInfoVecIndexMapVecVec.push_back(tmpBreakPoint2FusionInfoVecIndexMapVec);
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of initiating breakPoint2FusionInfoVecIndexMapVecVec" << endl;
	cout << endl << "[" << asctime(local) << "... start to merge breakPoint files" << endl;
	string breakPointFile_merged_raw = outputFilePrefix + "breakPoint.merged.raw";
	string cmd_cat_breakPointFile = "cat " + inputBreakPointFile_upstreamHead + " " + inputBreakPointFile_downstreamTail
		+ " > " + breakPointFile_merged_raw;
	system(cmd_cat_breakPointFile.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of merging breakPoint files" << endl;
	cout << endl << "[" << asctime(local) << "... start to build breakPoint2FusionInfoVecIndexMapVecVec" << endl;
	vector<int> chrNameIntVec_gene1;
	vector<int> breakPointPosVec_gene1;
	vector<int> chrNameIntVec_gene2;
	vector<int> breakPointPosVec_gene2;
	vector<string> strandVec_gene1;
	vector<string> strandVec_gene2;
	vector<string> transcriptIdVec_gene1;
	vector<string> transcriptIdVec_gene2;
	ifstream breakPoint_merged_raw_ifs(breakPointFile_merged_raw.c_str());
	int tmpFusionInfoIndex = 0;
	while(!breakPoint_merged_raw_ifs.eof())
	{
		string tmpStr;
		getline(breakPoint_merged_raw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t"); int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1); int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1); int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
		string tmpChrName_gene1 = tmpStr.substr(0, tabLoc_1);
		string tmpPosStr_gene1 = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpChrName_gene2 = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpPosStr_gene2 = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		
		int tmpChrNameInt_gene1 = indexInfo->convertStringToInt(tmpChrName_gene1);
		int tmpPos_gene1 = atoi(tmpPosStr_gene1.c_str());
		int tmpChrNameInt_gene2 = indexInfo->convertStringToInt(tmpChrName_gene2);
		int tmpPos_gene2 = atoi(tmpPosStr_gene2.c_str());		
		string tmpStrand_gene1 = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		string tmpStrand_gene2 = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		bool tmpForOrRev_bool_gene1, tmpForOrRev_bool_gene2;
		if(tmpStrand_gene1 == "+")
			tmpForOrRev_bool_gene1 = true;
		else if(tmpStrand_gene1 == "-")
			tmpForOrRev_bool_gene1 = false;
		else
		{
			cout << "invalid tmpStrand_gene1: " << tmpStrand_gene1 << endl;
			exit(1);
		}
		if(tmpStrand_gene2 == "+")
			tmpForOrRev_bool_gene2 = true;
		else if(tmpStrand_gene2 == "-")
			tmpForOrRev_bool_gene2 = false;
		else
		{
			cout << "invalid tmpStrand_gene2: " << tmpStrand_gene2 << endl;
			exit(1);
		}
		if((tmpChrNameInt_gene1 < 0)||(tmpChrNameInt_gene2 < 0))		
			continue;
		string tmpTranscriptId_gene1, tmpTranscriptId_gene2;
		exonBoundary2transcriptIdHashInfo.getFusionTranscriptIdPairFromBreakPoint(
			tmpChrNameInt_gene1, tmpChrNameInt_gene2, tmpPos_gene1, tmpPos_gene2, 
			tmpForOrRev_bool_gene1, tmpForOrRev_bool_gene2, tmpTranscriptId_gene1, tmpTranscriptId_gene2);
		if(withinTheSameTranscriptOrNot(tmpTranscriptId_gene1, tmpTranscriptId_gene2))
			continue;
		if(tmpChrNameInt_gene1 == tmpChrNameInt_gene2)
		{
			int tmpDistance = tmpPos_gene1 - tmpPos_gene2;
			if((tmpDistance <= 500000)&&(tmpDistance >= -500000))
				continue;
		}
		if(((breakPoint2FusionInfoVecIndexMapVecVec[tmpChrNameInt_gene1])[tmpChrNameInt_gene2]).find(tmpPos_gene1)
			== ((breakPoint2FusionInfoVecIndexMapVecVec[tmpChrNameInt_gene1])[tmpChrNameInt_gene2]).end()) // new fusion info
		{
			// put 2 fusionInfoVec
			chrNameIntVec_gene1.push_back(tmpChrNameInt_gene1);
			breakPointPosVec_gene1.push_back(tmpPos_gene1);
			chrNameIntVec_gene2.push_back(tmpChrNameInt_gene2);
			breakPointPosVec_gene2.push_back(tmpPos_gene2);
			strandVec_gene1.push_back(tmpStrand_gene1);
			strandVec_gene2.push_back(tmpStrand_gene2);
			transcriptIdVec_gene1.push_back(tmpTranscriptId_gene1);
			transcriptIdVec_gene2.push_back(tmpTranscriptId_gene2);
			// insert 2 mapVecVec
			((breakPoint2FusionInfoVecIndexMapVecVec[tmpChrNameInt_gene1])[tmpChrNameInt_gene2]).insert(
				pair<int, pair<int,int> >(tmpPos_gene1, pair<int,int>(tmpPos_gene2, tmpFusionInfoIndex)));
			tmpFusionInfoIndex ++;
		}
	}
	breakPoint_merged_raw_ifs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of building breakPoint2FusionInfoVecIndexMapVecVec" << endl;
	cout << endl << "[" << asctime(local) << "... start to print bedpe file and fusion breakpoint file" << endl;
	string outputFile_breakpoint = outputFilePrefix + "breakpoint.summary";
	string outputFile_bedpe = outputFilePrefix + "fusion.bedpe";
	ofstream breakpoint_ofs(outputFile_breakpoint.c_str());
	ofstream bedpe_ofs(outputFile_bedpe.c_str());
	for(int tmp = 0; tmp < chrNameIntVec_gene1.size(); tmp++)
	{
		int tmpChrNameInt_gene1 = chrNameIntVec_gene1[tmp];
		int tmpPos_gene1 = breakPointPosVec_gene1[tmp];
		int tmpChrNameInt_gene2 = chrNameIntVec_gene2[tmp];
		int tmpPos_gene2 = breakPointPosVec_gene2[tmp];		
		string tmpStrand_gene1 = strandVec_gene1[tmp];
		string tmpTranscriptId_gene1 = transcriptIdVec_gene1[tmp];
		string tmpStrand_gene2 = strandVec_gene2[tmp];
		string tmpTranscriptId_gene2 = transcriptIdVec_gene2[tmp];
		string tmpChrName_gene1 = indexInfo->returnChrNameStr(tmpChrNameInt_gene1);
		string tmpChrName_gene2 = indexInfo->returnChrNameStr(tmpChrNameInt_gene2);
		breakpoint_ofs << tmpChrName_gene1 << "\t" << tmpPos_gene1 << "\t"
			<< tmpChrName_gene2 << "\t" << tmpPos_gene2 << "\t"
			<< tmpTranscriptId_gene1 << "\t" << tmpTranscriptId_gene2 << "\t"
			<< tmpStrand_gene1 << "\t" << tmpStrand_gene2 << endl; 
		if(tmpStrand_gene1 == "+")
			bedpe_ofs << tmpChrName_gene1.substr(3) << "\t" << tmpPos_gene1 - 1 << "\t" << tmpPos_gene1 << "\t";
		else
			bedpe_ofs << tmpChrName_gene1.substr(3) << "\t" << tmpPos_gene1 << "\t" << tmpPos_gene1 + 1 << "\t";

		if(tmpStrand_gene2 == "+")
			bedpe_ofs << tmpChrName_gene2.substr(3) << "\t" << tmpPos_gene2 << "\t" << tmpPos_gene2 + 1 << "\t";
		else
			bedpe_ofs << tmpChrName_gene2.substr(3) << "\t" << tmpPos_gene2 - 1 << "\t" << tmpPos_gene2 << "\t";			
	
		bedpe_ofs <<  tmpTranscriptId_gene1 << "-" << tmpTranscriptId_gene2 << "\t0";
		if(tmpStrand_gene1 == "+")
			bedpe_ofs << "\t1";
		else
			bedpe_ofs << "\t-1";
		if(tmpStrand_gene2 == "+")
			bedpe_ofs << "\t1" << endl;
		else
			bedpe_ofs << "\t-1" << endl;	
	}
	bedpe_ofs.close();
	breakpoint_ofs.close();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of printing bedpe file and fusion breakpoint file" << endl;
	cout << endl << "[" << asctime(local) << "... All jobs done!" << endl;
	delete indexInfo;
	return 0;
}