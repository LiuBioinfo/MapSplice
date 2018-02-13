// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXHEADTAIL_PARALLELPROCESSES_H
#define FIXHEADTAIL_PARALLELPROCESSES_H

using namespace std;

void io_stage_fixHeadTail_separateThreadForIO(
	ifstream& inputRecord_ifs, AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue, bool& endOfFile_bool,
	bool& endOfProcessing_bool, int inputReadNumInBatchArray_fixHeadTail,
	int inputTimePerc_fixHeadTail, int outputTimePerc_fixHeadTail, ofstream& log_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs,
	ofstream& input_log_ofs, ofstream& output_log_ofs)
{
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;	

	int tmpBatchIndex_input = 0;
	int tmpBatchIndex_output = 0;
	int tmpThread = omp_get_thread_num();
	log_ofs << "input thread: " << tmpThread << endl;
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	input_log_ofs << endl << tmpTimeStr << "... input of FixHeadTail starts ......" << endl << endl; 
	string line1_1st, line2_1st, line3_1st, line4_1st, line5_1st, line6_1st, line7_1st, line8_1st, line9_1st, line10_1st, line11_1st;
	getline(inputRecord_ifs, line11_1st);
	getline(inputRecord_ifs, line1_1st);
	getline(inputRecord_ifs, line2_1st);
	getline(inputRecord_ifs, line3_1st);
	getline(inputRecord_ifs, line4_1st);
	getline(inputRecord_ifs, line5_1st);
	getline(inputRecord_ifs, line6_1st);
	getline(inputRecord_ifs, line7_1st);
	getline(inputRecord_ifs, line8_1st);
	getline(inputRecord_ifs, line9_1st);
	getline(inputRecord_ifs, line10_1st);
	alignInfoInputQueue->initiateWith1stAlignInfo(line1_1st, line2_1st, line3_1st, line4_1st, line5_1st, 
		line6_1st, line7_1st, line8_1st, line9_1st, line10_1st, input_log_ofs);
	tmpBatchIndex_input ++;
	bool input_stage_end_bool = false;
	bool output_stage_end_bool = false;
	while(1)
	{
		if(!input_stage_end_bool)
		{
			for(int tmp = 0; tmp < inputTimePerc_fixHeadTail; tmp++)
			{
				if(inputRecord_ifs.eof())
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					input_log_ofs << endl << tmpTimeStr << "... end of input of FixHeadTail ......" << endl << endl;  
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;					
				}
				string line11;
				getline(inputRecord_ifs, line11);
				if(inputRecord_ifs.eof())
				{
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;						
				}
				string line1, line2, line3, line4, line5, line6, line7, line8, line9, line10;
				getline(inputRecord_ifs, line1);
				getline(inputRecord_ifs, line2);
				getline(inputRecord_ifs, line3);
				getline(inputRecord_ifs, line4);
				getline(inputRecord_ifs, line5);
				getline(inputRecord_ifs, line6);
				getline(inputRecord_ifs, line7);
				getline(inputRecord_ifs, line8);
				getline(inputRecord_ifs, line9);
				getline(inputRecord_ifs, line10);
				alignInfoInputQueue->getAlignInfoFromInputFile(line1, line2, line3, line4, line5, line6, line7, line8, 
					line9, line10, inputReadNumInBatchArray_fixHeadTail, input_log_ofs, tmpBatchIndex_input);
			}
		}
		if(!output_stage_end_bool)
		{
			for(int tmp = 0; tmp < outputTimePerc_fixHeadTail; tmp++)
			{
				if(fixHeadTailResultQueue->atLeast3Node())
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					tmpBatchIndex_output ++;
					output_log_ofs << endl << tmpTimeStr << "... output regular array 2 file of fixHeadTail starts ......" << endl;  						
					output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
					fixHeadTailResultQueue->outputFrontResultArray(OutputSamFile_fixHeadTail_complete_pair_ofs,
						OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs,
						OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs);
					fixHeadTailResultQueue->popFromResultQueue();
				}
				else if(fixHeadTailResultQueue->only2Node())
				{
					if(endOfProcessing_bool)
					{
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						tmpBatchIndex_output ++;
						output_log_ofs << endl << tmpTimeStr << "... output last array 2 file of fixHeadTail starts ......" << endl;  						
						output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
						fixHeadTailResultQueue->outputFrontResultArray(OutputSamFile_fixHeadTail_complete_pair_ofs,
							OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs,
							OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs);
						fixHeadTailResultQueue->popFromResultQueue();
					}
					else
						continue;
				}
				else if(fixHeadTailResultQueue->resultQueueEmpty())
				{
					if(endOfProcessing_bool)
					{
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						output_log_ofs << endl << tmpTimeStr << "... end of output of fixHeadTail ......" << endl << endl;  
						output_stage_end_bool = true;
						break;	
					}
				}
				else
				{
					continue;
				}
			}
		}
		if(input_stage_end_bool && output_stage_end_bool)
			break;
	}
	log_ofs << "end of input/output !"  << endl;
}

void process_function_fixHeadTail(
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixHeadTail_Array* tmpResultFixHeadTailArray, int tmpBatchArraySize, 
	int threadNumForProcess, bool fasta_or_fastq_bool, Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom, vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress, vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild, Index_Info* indexInfo, SJhash_Info* SJ,
	bool Do_extendHeadTail_fixHeadTail, bool annotation_provided_bool, 
	bool Do_annotation_only_bool, Annotation_Info* annotationInfo,
	bool checkQualSeqForReadSegSeq, bool checkQualSeqForShortAnchorSeqToTargetMap,
	bool spliceJunctionHashExists, bool Do_fixHeadTail_remapping,
	bool Do_fixHeadTail_greedyMapping, bool Do_fixHeadTail_remappingAndTargetMapping,
	bool Do_fixHeadTail_remappingAgain, ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	unsigned int* sa_SNP, BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool)
{
	omp_set_num_threads(threadNumForProcess);
	#pragma omp parallel for schedule(dynamic)
	for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
	{
		////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
		int threadNO = omp_get_thread_num();
		string tmpReadNameAlignNum_1 = alignInfoInputQueue->returnFrontNodeReadNameAlignNum_1(tmpOpenMP);
		string tmpReadSeq_1= alignInfoInputQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
		string tmpReadQualSeq_1 = alignInfoInputQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
		string tmpReadNameAlignNum_2 = alignInfoInputQueue->returnFrontNodeReadNameAlignNum_2(tmpOpenMP);
		string tmpReadSeq_2= alignInfoInputQueue->returnFrontNodeReadSeq_2(tmpOpenMP);
		string tmpReadQualSeq_2 = alignInfoInputQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
		string tmpAlignInfo_Nor1 = alignInfoInputQueue->returnFrontNodeAlignInfo_Nor1(tmpOpenMP);
		string tmpAlignInfo_Rcm1 = alignInfoInputQueue->returnFrontNodeAlignInfo_Rcm1(tmpOpenMP);
		string tmpAlignInfo_Nor2 = alignInfoInputQueue->returnFrontNodeAlignInfo_Nor2(tmpOpenMP);
		string tmpAlignInfo_Rcm2 = alignInfoInputQueue->returnFrontNodeAlignInfo_Rcm2(tmpOpenMP);

		int multiMapSeg_maxLength = 0;
		PE_Read_Info peReadInfo;// = new PE_Read_Info();
		PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
		peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(tmpReadNameAlignNum_1, 
			tmpReadSeq_1, tmpReadQualSeq_1, tmpReadNameAlignNum_2, tmpReadSeq_2, tmpReadQualSeq_2,
			tmpAlignInfo_Nor1, tmpAlignInfo_Rcm1, tmpAlignInfo_Nor2, tmpAlignInfo_Rcm2, peReadInfo, 
			indexInfo, fasta_or_fastq_bool, SE_or_PE_bool, multiMapSeg_maxLength);		
		#ifdef MAP_INFO
		cout << "readName_1: " << peReadInfo.returnReadName_1() << endl;
		cout << "readName_2: " << peReadInfo.returnReadName_2() << endl;
		cout << "start fixHeadTail: " << endl;
		cout << "PeAlignInfo:" << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
		#endif
		FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
		if(Do_fixHeadTail_remapping)
		{
			fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(peReadInfo, peAlignInfo, SJ, 
				secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab,
				secondLevelDetChild, spliceJunctionHashExists, indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);	
		}
		#ifdef MAP_INFO
		cout << "after remapping:" << endl;
		cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
		#endif
		if(Do_fixHeadTail_greedyMapping)
		{
			if(do_segMap2snpMer_bool)
			{	
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly_includeSNPseqMap(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa, secondLevelLcpCompress,
					secondLevelChildTab, secondLevelDetChild, spliceJunctionHashExists, indexInfo, 
					Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, 
					annotationInfo, MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool, sa_SNP, 
					lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq);	
			}
			else
			{	
				fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
					peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa,
					secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
					spliceJunctionHashExists, indexInfo, Do_extendHeadTail_fixHeadTail,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool);	
			}
		}
		#ifdef MAP_INFO
		cout << "after greedyMapping:" << endl;
		cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
		#endif
		if(Do_fixHeadTail_remappingAndTargetMapping)
		{
			fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
				peReadInfo, peAlignInfo, SJ, secondLevelChrom, secondLevelSa,
				secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
				spliceJunctionHashExists, indexInfo, Do_extendHeadTail_fixHeadTail,
				checkQualSeqForShortAnchorSeqToTargetMap, SE_or_PE_bool);	
		}
		#ifdef MAP_INFO
		cout << "after remapping And Target Mapping:" << endl;
		cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;				
		#endif
		if(Do_fixHeadTail_remappingAgain)
		{
			fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
				peReadInfo, peAlignInfo, SJ, 
				secondLevelChrom,
				secondLevelSa,
				secondLevelLcpCompress,
				secondLevelChildTab,
				secondLevelDetChild,
				spliceJunctionHashExists,
				indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);
		}
		//fixHeadTailInfo->fixHeadTail_extend2end(peReadInfo, peAlignInfo, indexInfo);
		fixHeadTailInfo->fixHeadTail_extend2end_finalStepForAligner(
			peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);
			fixHeadTailInfo->fixHeadTail_extend2end_fixIndel(
			peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);							
		#ifdef MAP_INFO
		cout << "after extending: " << endl;
		cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;								
		#endif
				// remove duplicate mismatch
		peAlignInfo->removeDuplicateMismatch(SE_or_PE_bool);
		peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
		bool pairExistsBool = peAlignInfo->finalPairExistsBool();							

		if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
		{
			bool allFinalPairAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();	
			bool unique_bool = peAlignInfo->checkUniqueOrMulti();
			statsInfo->increPairedNum_fixHeadTail(threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
					
			bool align_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
			if(align_lowScore_bool)
			{
				statsInfo->increLowScoreComplete_fixHeadTail(threadNO, allFinalPairAlignmentCompleteBool, unique_bool);
				tmpResultFixHeadTailArray->insert_peAlignSamVec_pair_lowScore(
					peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
			}
			else if(allFinalPairAlignmentCompleteBool)
				tmpResultFixHeadTailArray->insert_peAlignSamVec_complete_pair(
					peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
			else
				tmpResultFixHeadTailArray->insert_peAlignSamVec_incomplete_pair(
					peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);	
		}
		else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
		{
			if(!outputUnpairedSAM_bool)
			{	
				bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
				statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool);
				if(allUnpairAlignmentCompleteBool)
					tmpResultFixHeadTailArray->insert_peAlignSamVec_complete_unpair(
						peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
				else	
					tmpResultFixHeadTailArray->insert_peAlignSamVec_incomplete_unpair(
						peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
			}
			else
			{
				peAlignInfo->chooseBestAlignment_final_PEasSE();
				bool allUnpairAlignmentCompleteBool_final = peAlignInfo->allFinalUnpairedAlignmentCompleted();
				statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool_final);				
				if(outputUnpairedSAM_bothEndsUniqueMappedOnly_bool)
				{
					if(allUnpairAlignmentCompleteBool_final)
						tmpResultFixHeadTailArray->insert_peAlignSamVec_complete_unpair(
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
								peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
					else
						tmpResultFixHeadTailArray->insert_peAlignSamVec_incomplete_unpair(
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
								peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
				}
				else
				{
					if(allUnpairAlignmentCompleteBool_final)
						tmpResultFixHeadTailArray->insert_peAlignSamVec_complete_unpair(
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
								peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
					else
						tmpResultFixHeadTailArray->insert_peAlignSamVec_incomplete_unpair(
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
								peReadInfo, fasta_or_fastq_bool), tmpOpenMP);					
				}
			}
		}
		delete fixHeadTailInfo;
		peAlignInfo->memoryFree();
		delete peAlignInfo;
	}	
}

void process_stage_fixHeadTail_separateThreadForIO(
	AlignInfoInput_Array_Queue* alignInfoInputQueue, 
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue, bool& endOfFile_bool, 
	bool& endOfProcessing_bool, int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom, vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress, vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild, Index_Info* indexInfo,  SJhash_Info* SJ,
	bool Do_extendHeadTail_fixHeadTail, bool annotation_provided_bool,  bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo, bool checkQualSeqForReadSegSeq,
	bool checkQualSeqForShortAnchorSeqToTargetMap, bool spliceJunctionHashExists,	
	bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
	bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
	ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	unsigned int* sa_SNP,  BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool)
{
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;	
	int tmpBatchIndex = 0;
	while(1)
	{
		int tmpBatchArraySize = 0;
		if(alignInfoInputQueue->atLeast3Node())
		{
			mapping_log_ofs << "at least 3 node in fixHeadTail" << endl;
			tmpBatchArraySize = alignInfoInputQueue->returnFrontNodeSize();
		}
		else if(alignInfoInputQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << "... only 2 node -- end of file in fixHeadTail" << endl;
				tmpBatchArraySize = alignInfoInputQueue->returnFrontNodeSize();
			}
			else
			{
				//mapping_log_ofs << "... only 2 node -- not end of file" << endl;
				continue;
			}
		}
		else if(alignInfoInputQueue->inputQueueEmpty())
		{
			//mapping_log_ofs << "inputQueueEmpty" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << "... inputQueueEmpty -- end of file in  fixHeadTail" << endl;
				break;
			}
			else
			{
				//mapping_log_ofs << "... inputQueueEmpty -- not end of file" << endl;
				continue;
			}
		}
		else
		{
			//mapping_log_ofs << "other cases" << endl;
			continue;
		}
			
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		tmpBatchIndex ++;
		mapping_log_ofs << endl << "*******************************************************************************"
			<< endl << tmpTimeStr << "... mapping a batch of reads starts in fixHeadTail ......" << endl << endl;  
		mapping_log_ofs << "In fixHeadTail tmpBatchIndex: " << tmpBatchIndex << endl;			
		mapping_log_ofs << "In fixHeadTail threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "In fixHeadTail tmpBatchArraySize: " << tmpBatchArraySize << endl;			

		Result_FixHeadTail_Array* tmpResultFixHeadTailArray = new Result_FixHeadTail_Array(tmpBatchArraySize);
		process_function_fixHeadTail(
			alignInfoInputQueue, tmpResultFixHeadTailArray, tmpBatchArraySize,
			threadNumForProcess, fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa,
			secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, SJ,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, 
			annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
			spliceJunctionHashExists, Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
			Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain,
			mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
			verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << endl << "*******************************************************************************" 
			<< endl << tmpTimeStr << "... mapping a batch of input alignInfo ends in fixHeadTail ......" << endl << endl; 		
		fixHeadTailResultQueue->pushBack2ResultArrayQueue(tmpResultFixHeadTailArray);
		alignInfoInputQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}

void process_stage_fixHeadTail_separateThreadForIO(
	AlignInfoInput_Array_Queue* alignInfoInputQueue, 
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue, bool& endOfFile_bool, 
	bool& endOfProcessing_bool, int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom, vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress, vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild, Index_Info* indexInfo,  SJhash_Info* SJ,
	bool Do_extendHeadTail_fixHeadTail, bool annotation_provided_bool,  bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo, bool checkQualSeqForReadSegSeq,
	bool checkQualSeqForShortAnchorSeqToTargetMap, bool spliceJunctionHashExists,	
	bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
	bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
	ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
	// provided SNPs
	unsigned int* sa_SNP,  BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool,
	// learned SNPs
	unsigned int* sa_SNP_learned,  BYTE* lcpCompress_SNP_learned, unsigned int* childTab_SNP_learned, char* chrom_SNP_learned, BYTE* verifyChild_SNP_learned, 
	Index_Info* indexInfo_SNP_learned, int SNPlocInSyntheticSNPseq_learned, bool do_segMap2snpMer_bool_learned, bool snpLearned_success_bool)
{
	if(do_segMap2snpMer_bool && (!do_segMap2snpMer_bool_learned))
	{	
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, SJ,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq,
			checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists,	Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
			Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP,  lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);
	}
	else if((!do_segMap2snpMer_bool) && do_segMap2snpMer_bool_learned)
	{
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, SJ,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq,
			checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists,	Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
			Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, snpLearned_success_bool);
	}
	else if((!do_segMap2snpMer_bool) && (!do_segMap2snpMer_bool_learned))
	{
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, SJ,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq,
			checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists,	Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
			Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP,  lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, // does not matter as SNPmerMapping won't be applied 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, // does not matter as SNPmerMapping won't be applied 
			false);
	}
	else
	{
		cout << "segMap2SNPmer_phase2_bool == true && segMap2SNPmer_phase2_bool_learned == true" << endl;
		cout << "Not allowed for MPS for now" << endl;
		exit(1);
	}		
}

void process_stage_fixHeadTail_separateThreadForIO(
	AlignInfoInput_Array_Queue* alignInfoInputQueue, 
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue, bool& endOfFile_bool, 
	bool& endOfProcessing_bool, int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom, vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress, vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild, Index_Info* indexInfo,
	////////////////////////////////////////////////////
	SJhash_Info* SJ_shared, SJhash_Info* SJ_tmpData,
	////////////////////////////////////////////////////
	bool Do_extendHeadTail_fixHeadTail, bool annotation_provided_bool,  bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
	bool checkQualSeqForReadSegSeq, bool checkQualSeqForShortAnchorSeqToTargetMap, 
	////////////////////////////////////////////////////
	bool spliceJunctionHashExists_shared, bool spliceJunctionHashExists_tmpData,
	////////////////////////////////////////////////////
	bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping, bool Do_fixHeadTail_remappingAndTargetMapping, 
	bool Do_fixHeadTail_remappingAgain, ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
	// provided SNPs
	unsigned int* sa_SNP,  BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool,
	// learned SNPs
	unsigned int* sa_SNP_learned,  BYTE* lcpCompress_SNP_learned, unsigned int* childTab_SNP_learned, char* chrom_SNP_learned, BYTE* verifyChild_SNP_learned, 
	Index_Info* indexInfo_SNP_learned, int SNPlocInSyntheticSNPseq_learned, bool do_segMap2snpMer_bool_learned, bool snpLearned_success_bool)
{
	if(spliceJunctionHashExists_shared && (!spliceJunctionHashExists_tmpData))
	{	
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo, 
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo,
			////////////////////////////////////////////////////
			SJ_shared,
			////////////////////////////////////////////////////
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, 
			////////////////////////////////////////////////////
			true, //spliceJunctionHashExists_shared == true,
			////////////////////////////////////////////////////
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else if((!spliceJunctionHashExists_shared) && spliceJunctionHashExists_tmpData)
	{
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo, 
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo,
			////////////////////////////////////////////////////
			SJ_tmpData,
			////////////////////////////////////////////////////
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, 
			////////////////////////////////////////////////////
			true, //spliceJunctionHashExists_tmpData == true,
			////////////////////////////////////////////////////
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else if((!spliceJunctionHashExists_shared) && (!spliceJunctionHashExists_tmpData))
	{
		process_stage_fixHeadTail_separateThreadForIO(
			alignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo, 
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo,
			////////////////////////////////////////////////////
			SJ_shared, // does not matter, as it won't be used
			////////////////////////////////////////////////////
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, 
			////////////////////////////////////////////////////
			false, // spliceJunctionHashExists_shared == false, spliceJunctionHashExists_tmpData == false
			////////////////////////////////////////////////////
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else
	{
		cout << "spliceJunctionHashExists_shared == true && spliceJunctionHashExists_tmpData == true" << endl;
		cout << "Not allowed for MPS for now" << endl;
		exit(1);
	}		
}

void io_process_fixHeadTail_allThreadsSharedByBothStage(
	ifstream& inputRecord_ifs, AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue,
	int inputReadNumInBatchArray_fixHeadTail, ofstream& log_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs,
	ofstream& input_log_ofs, ofstream& output_log_ofs,
	int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, SJhash_Info* SJ,
	bool Do_extendHeadTail_fixHeadTail,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	bool checkQualSeqForReadSegSeq,
	bool checkQualSeqForShortAnchorSeqToTargetMap,
	bool spliceJunctionHashExists,
	bool Do_fixHeadTail_remapping,
	bool Do_fixHeadTail_greedyMapping,
	bool Do_fixHeadTail_remappingAndTargetMapping,
	bool Do_fixHeadTail_remappingAgain,
	ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	unsigned int* sa_SNP, BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, char* chrom_SNP,
	BYTE* verifyChild_SNP, Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool)
{
	int tmpBatchIndex = 0;
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;	
	bool endOfFile_bool = false;
	while(1)
	{
		if(endOfFile_bool)
			break;
		// input stage
		tmpBatchIndex ++;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		input_log_ofs << "***********************************************************" << endl;
		input_log_ofs << tmpTimeStr << "... start to input reads ..." << endl;
		input_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		for(int tmp = 0; tmp < inputReadNumInBatchArray_fixHeadTail; tmp++)
		{
			if(inputRecord_ifs.eof())
			{
				nowtime = time(NULL);
				local = localtime(&nowtime);
				tmpTimeStr = asctime(local);
				tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
				input_log_ofs << endl << tmpTimeStr << "... end of input of FixHeadTail ......" << endl << endl;  
				endOfFile_bool = true;
				break;					
			}
			string line11;
			getline(inputRecord_ifs, line11);
			if(inputRecord_ifs.eof())
			{
				endOfFile_bool = true;
				break;						
			}
			string line1, line2, line3, line4, line5, line6, line7, line8, line9, line10;
			getline(inputRecord_ifs, line1);
			getline(inputRecord_ifs, line2);
			getline(inputRecord_ifs, line3);
			getline(inputRecord_ifs, line4);
			getline(inputRecord_ifs, line5);
			getline(inputRecord_ifs, line6);
			getline(inputRecord_ifs, line7);
			getline(inputRecord_ifs, line8);
			getline(inputRecord_ifs, line9);
			getline(inputRecord_ifs, line10);
			if(tmp == 0)
				alignInfoInputQueue->initiateWith1stAlignInfo(line1, line2, line3, line4, line5,  line6, line7, line8, line9, line10, input_log_ofs);
			else
				alignInfoInputQueue->pushBack2currentFrontAlignInfoArray(line1, line2, line3, line4, line5,  line6, line7, line8, line9, line10);
		}
		if(alignInfoInputQueue->inputQueueEmpty())
		{
			input_log_ofs << "inputQueueEmpty ! end of fixHeadTail !" << endl;
			break;			
		}
		int tmpBatchArraySize = alignInfoInputQueue->returnFrontNodeSize();
		input_log_ofs << "tmpBatchArraySize: " << tmpBatchArraySize << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		input_log_ofs << tmpTimeStr << "... end of inputting reads ..." << endl << endl;

		// process stage
		mapping_log_ofs << "***********************************************************" << endl;
		mapping_log_ofs << tmpTimeStr << "... start to map reads ......" << endl;
		mapping_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		mapping_log_ofs << "ThreadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "BatchArraySize: " << tmpBatchArraySize << endl;	
		Result_FixHeadTail_Array* tmpResultFixHeadTailArray = new Result_FixHeadTail_Array(tmpBatchArraySize);
		process_function_fixHeadTail(
			alignInfoInputQueue, tmpResultFixHeadTailArray, tmpBatchArraySize, threadNumForProcess, 
			fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, 
			secondLevelChildTab, secondLevelDetChild, indexInfo, SJ, Do_extendHeadTail_fixHeadTail, 
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, 
			checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists, Do_fixHeadTail_remapping, 
			Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, 
			mapping_log_ofs, SE_or_PE_bool, 
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << tmpTimeStr << "... end of mapping reads ......" << endl << endl;	
		fixHeadTailResultQueue->pushBack2ResultArrayQueue(tmpResultFixHeadTailArray);
		alignInfoInputQueue->popFromReadQueue();
		// output stage
		output_log_ofs << "***********************************************************" << endl;
		output_log_ofs << tmpTimeStr << "... start to output alignment ......" << endl;
		output_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		fixHeadTailResultQueue->outputFrontResultArray(OutputSamFile_fixHeadTail_complete_pair_ofs, 
			OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs, 
			OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		output_log_ofs << tmpTimeStr << "... end of outputting alignment ......" << endl << endl;
		fixHeadTailResultQueue->popFromResultQueue();
	}
}

void io_process_fixHeadTail_allThreadsSharedByBothStage(
	ifstream& inputRecord_ifs, AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue,
	int inputReadNumInBatchArray_fixHeadTail, ofstream& log_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs,
	ofstream& input_log_ofs, ofstream& output_log_ofs,
	int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, SJhash_Info* SJ,
	bool Do_extendHeadTail_fixHeadTail,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	bool checkQualSeqForReadSegSeq,
	bool checkQualSeqForShortAnchorSeqToTargetMap,
	bool spliceJunctionHashExists,
	bool Do_fixHeadTail_remapping,
	bool Do_fixHeadTail_greedyMapping,
	bool Do_fixHeadTail_remappingAndTargetMapping,
	bool Do_fixHeadTail_remappingAgain,
	ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
	// provided SNPs
	unsigned int* sa_SNP, BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool,
	// learned SNPs
	unsigned int* sa_SNP_learned, BYTE* lcpCompress_SNP_learned, unsigned int* childTab_SNP_learned, char* chrom_SNP_learned, BYTE* verifyChild_SNP_learned, 
	Index_Info* indexInfo_SNP_learned, int SNPlocInSyntheticSNPseq_learned, bool do_segMap2snpMer_bool_learned, bool snpLearned_success_bool)
{
	if(do_segMap2snpMer_bool && (!do_segMap2snpMer_bool_learned))
	{	
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs,
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs,
			OutputSamFile_fixHeadTail_complete_unpair_ofs, OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
			OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, SJ, Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
			checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists, Do_fixHeadTail_remapping,
			 Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			 outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);
	}
	else if((!do_segMap2snpMer_bool) && do_segMap2snpMer_bool_learned)
	{
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs,
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs,
			OutputSamFile_fixHeadTail_complete_unpair_ofs, OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
			OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, SJ, Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
			checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists, Do_fixHeadTail_remapping,
			 Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			 outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, snpLearned_success_bool);
	}
	else if((!do_segMap2snpMer_bool) && (!do_segMap2snpMer_bool_learned))
	{	
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs,
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs,
			OutputSamFile_fixHeadTail_complete_unpair_ofs, OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
			OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, fasta_or_fastq_bool,
			statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, SJ, Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
			checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists, Do_fixHeadTail_remapping,
			 Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			 outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, // does not matter as SNPmerMapping won't be applied
			indexInfo_SNP, SNPlocInSyntheticSNPseq, // does not matter as SNPmerMapping won't be applied
			false);
	}	
	else
	{
		cout << "segMap2SNPmer_phase2_bool == true && segMap2SNPmer_phase2_bool_learned == true" << endl;
		cout << "Not allowed for MPS for now" << endl;
		exit(1);
	}		
}

void io_process_fixHeadTail_allThreadsSharedByBothStage(
	ifstream& inputRecord_ifs, AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue,
	int inputReadNumInBatchArray_fixHeadTail, ofstream& log_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
	ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
	ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs,
	ofstream& input_log_ofs, ofstream& output_log_ofs,
	int threadNumForProcess, bool fasta_or_fastq_bool,
	Stats_Info* statsInfo, vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	////////////////////////////////////////////////////////////
	SJhash_Info* SJ_shared, SJhash_Info* SJ_tmpData,
	////////////////////////////////////////////////////////////
	bool Do_extendHeadTail_fixHeadTail,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	bool checkQualSeqForReadSegSeq,
	bool checkQualSeqForShortAnchorSeqToTargetMap,
	////////////////////////////////////////////////////////////
	bool spliceJunctionHashExists_shared, bool spliceJunctionHashExists_tmpData,
	////////////////////////////////////////////////////////////
	bool Do_fixHeadTail_remapping,
	bool Do_fixHeadTail_greedyMapping,
	bool Do_fixHeadTail_remappingAndTargetMapping,
	bool Do_fixHeadTail_remappingAgain,
	ofstream& mapping_log_ofs, bool SE_or_PE_bool,
	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
	// provided SNPs
	unsigned int* sa_SNP, BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, char* chrom_SNP, BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool,
	// learned SNPs
	unsigned int* sa_SNP_learned, BYTE* lcpCompress_SNP_learned, unsigned int* childTab_SNP_learned, char* chrom_SNP_learned, BYTE* verifyChild_SNP_learned, 
	Index_Info* indexInfo_SNP_learned, int SNPlocInSyntheticSNPseq_learned, bool do_segMap2snpMer_bool_learned, bool snpLearned_success_bool)
{
	if(spliceJunctionHashExists_shared && (!spliceJunctionHashExists_tmpData))
	{
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs, 
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs, 
			OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, 
			fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
			SJ_shared,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
			true, //spliceJunctionHashExists_shared == true,
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned,
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else if((!spliceJunctionHashExists_shared) && spliceJunctionHashExists_tmpData)
	{
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs, 
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs, 
			OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, 
			fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
			SJ_tmpData,
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
			true, // spliceJunctionHashExists_tmpData == true,
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else if((!spliceJunctionHashExists_shared) && (!spliceJunctionHashExists_tmpData))
	{
		io_process_fixHeadTail_allThreadsSharedByBothStage(
			inputRecord_ifs, alignInfoInputQueue, fixHeadTailResultQueue, inputReadNumInBatchArray_fixHeadTail, log_ofs, 
			OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs, 
			OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs, output_log_ofs, threadNumForProcess, 
			fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
			SJ_shared, // does not matter, as it won't be used, 
			Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
			false, // spliceJunctionHashExists_shared == false, spliceJunctionHashExists_tmpData == false
			Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNPs
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			// learned SNPs
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned,
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, do_segMap2snpMer_bool_learned, snpLearned_success_bool);
	}
	else
	{
		cout << "spliceJunctionHashExists_shared == true && spliceJunctionHashExists_tmpData == true" << endl;
		cout << "Not allowed for MPS for now" << endl;
		exit(1);
	}		
}


#endif