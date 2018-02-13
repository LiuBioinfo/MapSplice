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
#include "../../../otherProjects/indexing/utils/build_snpMerIndex.h"
#include "../../../general/extractUnmapAlignment2ReadFile.h"
#include "../../../phase1/arrayQueue_phase1.h"
#include "../../../phase2/arrayQueue_phase2.h"
#include "../../../stats_info.h"
#include "../../../constantDefinitions.h"
#include "../../../general/option_info.h"
#include "../../../general/read_block_test.h"
#include "../../../general/bwtmap_info.h"
#include "../../../general/DoubleAnchorScore.h"
#include "../../../general/sbndm.h"
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"
#include "../../../general/enhanced_suffix_array_info.h"
#include "../../../general/annotation_info.h"
#include "../../../phase1/repeatRegion.h"
#include "../../../general/segmentMapping.h"
#include "../../../general/splice_info.h"
#include "../../../general/fixGapRelationParameters.h"
#include "../../../general/read_info.h"
#include "../../../general/seg_info.h"
#include "../../../general/fixDoubleAnchorNWDP_info.h"
#include "../../../general/fixDoubleAnchorMatch_info.h"
#include "../../../general/fixDoubleAnchorInsertion_info.h"
#include "../../../general/fixDoubleAnchorDeletion_info.h"
#include "../../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../../general/fixDoubleAnchorSplice_info.h"
#include "../../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../../general/path_info.h"
#include "../../../general/gap_info.h"
#include "../../../otherProjects/incorporateGenomicVariants/general/syntheticSNPtransSeq_info.h"
#include "../../../otherProjects/incorporateGenomicVariants/general/learnedCandiSNPhash_info_vec.h"
#include "../../../general/align_info.h"
#include "../../../general/peAlign_info.h"
#include "../../../general/groupSeg_info.h"
#include "../../../general/alignInferJunctionHash_info_vec.h"
#include "../../../phase2/spliceJunctionHash_info.h"
#include "../../../phase2/unmapEnd_info.h"
#include "../../../phase2/unfixedHead.h"
#include "../../../phase2/unfixedTail.h"
#include "../../../phase2/incompleteLongHead.h"
#include "../../../phase2/incompleteLongTail.h"
#include "../../../phase2/sam2junc.h"
#include "../../../fixHeadTail.h"
#include "../../../phase2/fixOneEndUnmapped.h"
#include "../../../fixPhase1.h"
#include "../../../general/readSeqPreProcessing.h"
#include "../../../general/headerSection_info.h"
#include "../../../general/otherFunc2.h"
#include "../general/peSam_info.h"
#include "../general/fixFusion_peRead.h"
#include "../general/fixFusionResult_info.h"
using namespace std;

time_t nowtime;
struct tm *local;

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolder GeneAnnEntryFile threads_num inputIncompletePeSamFile outputFolder" << endl;
		exit(1);
	}
	int normalRecordNum = 1000000;
	string threads_num_str = argv[3];
	int threads_num = atoi(threads_num_str.c_str());
	bool fasta_or_fastq_bool = true;
	bool SE_or_PE_bool = false;
	int maximum_allowed_mismatchNum = 3;
	int minimum_allowed_insertionLength = -1; 
	int maximum_allowed_insertionLength = 0;
	bool geneAnnEntryBoundaryOnly_bool = false;
	bool geneAnnIncorporated_bool = true; 
	bool insertionAllowed_bool = true;

	/////////////////////////////////////////////////////
	bool Do_cirRNA = false;
	bool annotation_provided_bool = false;
	bool Do_annotation_only_bool = false;
	bool Do_extendHeadTail_phase1 = true;
	bool checkQualSeqForReadSegSeq = false;
	Annotation_Info* annotationInfo = new Annotation_Info();
	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}
	////////////////////////////////////////////////////
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
   	string settingsLogStr = outputFolderStr + "/settings.log";
   	ofstream log_ofs(settingsLogStr.c_str());
   	string runtimeLogStr = outputFolderStr + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	log_ofs << "CommandLine:" << endl;
   	for(int tmp = 1; tmp < argc; tmp++)
   		log_ofs << "\t" << argv[tmp];
   	log_ofs << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load whole genome index ......" << endl << endl; 
    string indexStr = argv[1];
    string preIndexArrayPreStr = indexStr; 	
    preIndexArrayPreStr.append("/");
    indexStr.append("/");
	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); 
	ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); 
	ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); 
	ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; 
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); 
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; 
	preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; 
	preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, log_ofs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	parameter_file_ifs.close();
	chrom_bit_file_ifs.close();
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
 	
	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 	

	ifstream SA_file_ifs(SA_file.c_str(),ios::binary);
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);

	nowtime = time(NULL); local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... start to load enhanced Suffix Array ......" << endl << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	SA_file_ifs.close();
	lcpCompress_file_ifs.close();
	childTab_file_ifs.close();
	verifyChild_file_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL); local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... all index loaded ......" << endl << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load simplified gene annotation file" << endl;
	string inputGeneAnnEntryFile = argv[2];	
	GeneAnnEntry_Hash_Info tmpGeneAnnHashInfo;
	tmpGeneAnnHashInfo.initiate_geneAnnEntryArea2infoIndexMapVec(chromNum);
	tmpGeneAnnHashInfo.loadGeneAnn(inputGeneAnnEntryFile, indexInfo);
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... end of loading simplified gene annotation file" << endl;


	string inputIncompleteUniquePairedAlignmentPath = argv[4];
	ifstream incompleteUniquePairedAlignment_ifs(inputIncompleteUniquePairedAlignmentPath.c_str());
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	int readTotalNum = 0;

	vector<string> inputSamStr1Vec(normalRecordNum);
	vector<string> inputSamStr2Vec(normalRecordNum);
	vector<string> outputClippedSam_unfixedHeadAtUpstreamEnd(normalRecordNum);
	vector<string> outputClippedSam_unfixedTailAtDownstreamEnd(normalRecordNum);
	vector<string> outputMainBodySam_upstreamEnd(normalRecordNum);
	vector<string> outputMainBodySam_downstreamEnd(normalRecordNum);
	vector<string> outputBreakPointPair_upstreamEnd(normalRecordNum);
	vector<string> outputBreakPointPair_downstreamEnd(normalRecordNum);

	string nonFusionSam_file = outputFolderStr + "/nonFusion.sam";
	string fusionSam_file = outputFolderStr + "/fusion.sam";
	string breakPoint_file_raw = outputFolderStr + "/fusion.breakPoint.raw";
	ofstream nonFusionSam_ofs(nonFusionSam_file.c_str());
	ofstream fusionSam_ofs(fusionSam_file.c_str());
	ofstream breakPoint_ofs_raw(breakPoint_file_raw.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Mapping process starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Mapping process starts ......" << endl << endl; 

	string line1, line2;
	for(tmpTurn = 0; 
		//tmpTurn <= 300     //used to control # of rounds to process
		; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		//cout << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		
		realRecordNum = normalRecordNum;
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if(incompleteUniquePairedAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
    		}
    		getline(incompleteUniquePairedAlignment_ifs, line1); // readName_1
    		if(incompleteUniquePairedAlignment_ifs.eof())
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;
    		}
    		inputSamStr1Vec[recordNumTmp] = line1;
    		getline(incompleteUniquePairedAlignment_ifs, line2);		
    		inputSamStr2Vec[recordNumTmp] = line2;
		}
		
		readTotalNum += realRecordNum;

		runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		runtime_log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;

		// cout << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		// cout << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;
		// cout << "realRecordNum: " << realRecordNum << endl;
		// cout << "threads_num: " << threads_num << endl;

		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();
			//cout << "tmpOpenMP: " << tmpOpenMP << endl;
			PeSam_Info tmpPeSamInfo;
			bool initiatePeSamInfo_bool = tmpPeSamInfo.initiateWith2samStr(inputSamStr1Vec[tmpOpenMP], inputSamStr2Vec[tmpOpenMP], indexInfo);
			if(!initiatePeSamInfo_bool)
				continue;
			//cout << "tmpSam_1: " << inputSamStr1Vec[tmpOpenMP] << endl;
			//cout << "tmpSam_2: " << inputSamStr2Vec[tmpOpenMP] << endl;
			PE_Read_Info tmpPeReadInfo;
			tmpPeReadInfo.initiateReadInfo(tmpPeSamInfo.returnRawReadName_1(), tmpPeSamInfo.returnRawReadName_2(), tmpPeSamInfo.returnRawReadSeq_1(), 
				tmpPeSamInfo.returnRawReadSeq_2(), tmpPeSamInfo.returnRawQualSeq_1(), tmpPeSamInfo.returnRawQualSeq_2(), fasta_or_fastq_bool, SE_or_PE_bool);
			bool unique_or_not_bool = tmpPeSamInfo.unique_or_not_bool();
			if(!unique_or_not_bool)
				continue;
			FixFusionResult_Info tmpFixFusionResultInfo;
			tmpFixFusionResultInfo.initiate_PE();
			bool unfixedHeadAtUpstreamEnd_exists_bool = tmpPeSamInfo.return_unfixedHeadAtUpstreamEnd_exists_bool();
			bool unfixedTailAtDownstreamEnd_exists_bool = tmpPeSamInfo.return_unfixedTailAtDownstreamEnd_exists_bool();
			//cout << "unfixedHeadAtUpstreamEnd_exists_bool: " << unfixedHeadAtUpstreamEnd_exists_bool << endl;
			//cout << "unfixedTailAtDownstreamEnd_exists_bool: " << unfixedTailAtDownstreamEnd_exists_bool << endl; 
			outputClippedSam_unfixedHeadAtUpstreamEnd[tmpOpenMP] = "";			
			outputClippedSam_unfixedTailAtDownstreamEnd[tmpOpenMP] = "";		
			outputMainBodySam_upstreamEnd[tmpOpenMP] = "";
			outputMainBodySam_downstreamEnd[tmpOpenMP] = "";
			outputBreakPointPair_upstreamEnd[tmpOpenMP] = "";
			outputBreakPointPair_downstreamEnd[tmpOpenMP] = "";
			string tmpNullStr = "*";
			// fix fusion at unfixedHeadAtUpstreamEnd
			if(unfixedHeadAtUpstreamEnd_exists_bool)
			{
				PE_Read_Info readInfo_unfixedLeftReadHead;
				FixPhase1Info fixPhase1Info_unfixedLeftReadHead;
				PE_Read_Alignment_Info peAlignInfo_unfixedLeftReadHead;
				string tmpReadName_upstreamEnd = tmpPeSamInfo.returnReadName_upstreamEnd();
				string tmpUnfixedLeftReadHead_readSeq = tmpPeSamInfo.returnUnfixedUpstreamReadHead_readSeq();
				string tmpUnfixedLeftReadHead_qualSeq = tmpPeSamInfo.returnUnfixedUpstreamReadHead_qualSeq(fasta_or_fastq_bool);
				readInfo_unfixedLeftReadHead.initiateReadInfo(tmpReadName_upstreamEnd, tmpNullStr, tmpUnfixedLeftReadHead_readSeq, tmpNullStr, tmpUnfixedLeftReadHead_qualSeq, tmpNullStr, fasta_or_fastq_bool, true);
				fixPhase1Info_unfixedLeftReadHead.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, 
					preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo_unfixedLeftReadHead, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
				fixPhase1Info_unfixedLeftReadHead.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedLeftReadHead, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
				fixPhase1Info_unfixedLeftReadHead.fixPhase1_gapInfo(readInfo_unfixedLeftReadHead, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);			
				peAlignInfo_unfixedLeftReadHead.initiatePeAlignInfo(fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor1, fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm1, 
					fixPhase1Info_unfixedLeftReadHead.pathInfo_Nor2, fixPhase1Info_unfixedLeftReadHead.pathInfo_Rcm2, indexInfo, true);			
				peAlignInfo_unfixedLeftReadHead.chooseBestAlignment_final_SE();
				bool validForFusionPostAnalysis_bool_unfixedLeftReadHead = peAlignInfo_unfixedLeftReadHead.return_validForFusionPostAnalysis_bool();
				//cout << "validForFusionPostAnalysis_bool_unfixedLeftReadHead: " << validForFusionPostAnalysis_bool_unfixedLeftReadHead << endl;
				if(validForFusionPostAnalysis_bool_unfixedLeftReadHead)
				{
					FixFusion_PeRead tmpFixFusionPeReadInfo_unfixedLeftReadHead;
					tmpFixFusionPeReadInfo_unfixedLeftReadHead.initiate_clippedSegGlobalMap(tmpPeSamInfo, peAlignInfo_unfixedLeftReadHead, true, tmpPeReadInfo, indexInfo);
					tmpFixFusionPeReadInfo_unfixedLeftReadHead.fixFusion(tmpPeReadInfo, tmpGeneAnnHashInfo, indexInfo, geneAnnEntryBoundaryOnly_bool, 
						geneAnnIncorporated_bool, insertionAllowed_bool, maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength);
					bool tmpFusionFixed_unfixedLeftReadHead_success_bool = tmpFixFusionPeReadInfo_unfixedLeftReadHead.return_fusionFixed_bool();
					//cout << "tmpFusionFixed_unfixedLeftReadHead_success_bool: " << tmpFusionFixed_unfixedLeftReadHead_success_bool << endl;
					if(tmpFusionFixed_unfixedLeftReadHead_success_bool)
						tmpFixFusionResultInfo.updateFixedFusionResult_unfixedHeadAtUpstreamEnd(readInfo_unfixedLeftReadHead, 
							peAlignInfo_unfixedLeftReadHead, tmpFixFusionPeReadInfo_unfixedLeftReadHead, fasta_or_fastq_bool, indexInfo);
				}
				fixPhase1Info_unfixedLeftReadHead.memoryFree();
				peAlignInfo_unfixedLeftReadHead.memoryFree();
			}
			//cout << "ends of fixing fusion at unfixedHeadAtUpstreamEnd" << endl;
			// fix fusion at unfixedTailAtDownstreamEnd
			if(unfixedTailAtDownstreamEnd_exists_bool)
			{
				PE_Read_Info readInfo_unfixedRightReadTail;
				FixPhase1Info fixPhase1Info_unfixedRightReadTail;
				PE_Read_Alignment_Info peAlignInfo_unfixedRightReadTail;
				string tmpReadName_downstreamEnd = tmpPeSamInfo.returnReadName_downstreamEnd();
				string tmpUnfixedRightReadTail_readSeq = tmpPeSamInfo.returnUnfixedDownstreamReadTail_readSeq();
				string tmpUnfixedRightReadTail_qualSeq = tmpPeSamInfo.returnUnfixedDownstreamReadTail_qualSeq(fasta_or_fastq_bool);
				readInfo_unfixedRightReadTail.initiateReadInfo(tmpReadName_downstreamEnd, tmpNullStr, tmpUnfixedRightReadTail_readSeq, tmpNullStr, tmpUnfixedRightReadTail_qualSeq, tmpNullStr, fasta_or_fastq_bool, true);
				fixPhase1Info_unfixedRightReadTail.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, 
					preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo_unfixedRightReadTail, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, true);
				fixPhase1Info_unfixedRightReadTail.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo_unfixedRightReadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, true);
				fixPhase1Info_unfixedRightReadTail.fixPhase1_gapInfo(readInfo_unfixedRightReadTail, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, true);
				peAlignInfo_unfixedRightReadTail.initiatePeAlignInfo(fixPhase1Info_unfixedRightReadTail.pathInfo_Nor1, fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm1, 
					fixPhase1Info_unfixedRightReadTail.pathInfo_Nor2, fixPhase1Info_unfixedRightReadTail.pathInfo_Rcm2, indexInfo, true);
				peAlignInfo_unfixedRightReadTail.chooseBestAlignment_final_SE();
				bool validForFusionPostAnalysis_bool_unfixedRightReadTail = peAlignInfo_unfixedRightReadTail.return_validForFusionPostAnalysis_bool();
				//cout << "validForFusionPostAnalysis_bool_unfixedRightReadTail: " << validForFusionPostAnalysis_bool_unfixedRightReadTail << endl;
				if(validForFusionPostAnalysis_bool_unfixedRightReadTail)
				{
					//cout << "initiating FixFusion_PeRead ...." << endl;
					FixFusion_PeRead tmpFixFusionPeReadInfo_unfixedRightReadTail;
					tmpFixFusionPeReadInfo_unfixedRightReadTail.initiate_clippedSegGlobalMap(tmpPeSamInfo, peAlignInfo_unfixedRightReadTail, false, tmpPeReadInfo, indexInfo);
					//cout << "start to fixFusion" << endl;
					tmpFixFusionPeReadInfo_unfixedRightReadTail.fixFusion(tmpPeReadInfo, tmpGeneAnnHashInfo, indexInfo, geneAnnEntryBoundaryOnly_bool, 
						geneAnnIncorporated_bool, insertionAllowed_bool, maximum_allowed_mismatchNum, minimum_allowed_insertionLength, maximum_allowed_insertionLength);
					//cout << "end of fixing fusion" << endl;
					bool tmpFusionFixed_unfixedRightReadTail_success_bool = tmpFixFusionPeReadInfo_unfixedRightReadTail.return_fusionFixed_bool();
					//cout << "tmpFusionFixed_unfixedRightReadTail_success_bool: " << tmpFusionFixed_unfixedRightReadTail_success_bool << endl;
					if(tmpFusionFixed_unfixedRightReadTail_success_bool)
						tmpFixFusionResultInfo.updateFixedFusionResult_unfixedTailAtDownstreamEnd(readInfo_unfixedRightReadTail, 
							peAlignInfo_unfixedRightReadTail, tmpFixFusionPeReadInfo_unfixedRightReadTail, fasta_or_fastq_bool, indexInfo);
				}
				fixPhase1Info_unfixedRightReadTail.memoryFree();
				peAlignInfo_unfixedRightReadTail.memoryFree();
			}
			//cout << "ends of fixing fusion at unfixedTailAtDownstreamEnd" << endl;
			//cout << "start to update_unfixedEndSam_mainBodySam_breakPointPairStr..." << endl;
			tmpFixFusionResultInfo.update_unfixedEndSam_mainBodySam_breakPointPairStr(tmpPeSamInfo, tmpPeReadInfo,
				outputClippedSam_unfixedHeadAtUpstreamEnd[tmpOpenMP], outputClippedSam_unfixedTailAtDownstreamEnd[tmpOpenMP],
				outputMainBodySam_upstreamEnd[tmpOpenMP], outputMainBodySam_downstreamEnd[tmpOpenMP], outputBreakPointPair_upstreamEnd[tmpOpenMP], 
				outputBreakPointPair_downstreamEnd[tmpOpenMP], indexInfo, tmpGeneAnnHashInfo);
		}

		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		runtime_log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		//cout << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		//cout << "start to output ... turn: " << tmpTurn+1 << endl;

		for(int tmp = 0; tmp < realRecordNum; tmp++)
		{
			if((outputClippedSam_unfixedHeadAtUpstreamEnd[tmp] == "")&&(outputClippedSam_unfixedTailAtDownstreamEnd[tmp] == ""))
				nonFusionSam_ofs << inputSamStr1Vec[tmp] << endl << inputSamStr2Vec[tmp] << endl;
			else
			{
				if(outputClippedSam_unfixedHeadAtUpstreamEnd[tmp] != "")
					fusionSam_ofs << outputClippedSam_unfixedHeadAtUpstreamEnd[tmp] << endl;
				fusionSam_ofs << outputMainBodySam_upstreamEnd[tmp] << endl << outputMainBodySam_downstreamEnd[tmp] << endl;
				if(outputClippedSam_unfixedTailAtDownstreamEnd[tmp] != "")
				 	fusionSam_ofs << outputClippedSam_unfixedTailAtDownstreamEnd[tmp] << endl;
				if(outputBreakPointPair_upstreamEnd[tmp] != "")
					breakPoint_ofs_raw << outputBreakPointPair_upstreamEnd[tmp] << endl;
				if(outputBreakPointPair_downstreamEnd[tmp] != "")
					breakPoint_ofs_raw << outputBreakPointPair_downstreamEnd[tmp] << endl;
			}
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		runtime_log_ofs << endl << "[" << asctime(local) << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;
		//cout << endl << "[" << asctime(local) << "finish output in phase1, turn: " << tmpTurn + 1 << endl << endl;		
	}
	log_ofs << "readTotalNum: " << readTotalNum << endl;
	breakPoint_ofs_raw.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Mapping process ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Mapping process ends ......" << endl << endl; 

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Generating break point file starts ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Generating break point file starts ......" << endl << endl;


	ifstream breakPointRaw_ifs(breakPoint_file_raw.c_str());
	string breakPoint_file = outputFolderStr + "/fusion.breakPoint";
	string breakPoint_file_filtered = outputFolderStr + "/fusion.breakPoint.filtered";
	ofstream breakPoint_ofs(breakPoint_file.c_str());
	ofstream breakPoint_ofs_filtered(breakPoint_file_filtered.c_str());
	while(!breakPointRaw_ifs.eof())
	{
		string tmpStr;
		getline(breakPointRaw_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr);
		string tmpStrand_1 = tmpFieldVec[5];
		string tmpStrand_2 = tmpFieldVec[6];
		if((tmpStrand_1 == "O")||(tmpStrand_1 == "X")||(tmpStrand_2 == "O")||(tmpStrand_2 == "X"))
			breakPoint_ofs_filtered << tmpStr << endl;
		else
			breakPoint_ofs << tmpStr << endl;
	}
	breakPoint_ofs_filtered.close();
	breakPoint_ofs.close();
	breakPointRaw_ifs.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... Generating break point file ends ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... Generating break point file ends ......" << endl << endl;

	delete annotationInfo;
	nonFusionSam_ofs.close();
	fusionSam_ofs.close();
	incompleteUniquePairedAlignment_ifs.close();
	runtime_log_ofs.close();
	free(preIndexMapLengthArray);
	free(preIndexIntervalStartArray);
	free(preIndexIntervalEndArray);	
	delete indexInfo;
	free(sa);
	free(childTab);
	free(lcpCompress);
	free(verifyChild);
	free(chrom);
   	log_ofs.close();
	return 0;
}