// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXONEENDUNMAPPED_PARALLELPROCESSES_H
#define FIXONEENDUNMAPPED_PARALLELPROCESSES_H

using namespace std;

void io_stage_fixOneEndUnmapped_separateThreadForIO(
	ifstream& inputRecord_ifs,
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,	
	int inputReadNumInBatchArray_fixOneEndUnmapped,
	int inputTimePerc_fixOneEndUnmapped,
	int outputTimePerc_fixOneEndUnmapped,
	ofstream& log_ofs,
	ofstream& OutputSamFile_oneEndMapped_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& OutputSamFile_oneEndMapped_unpair_ofs,
	ofstream& OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,	
	ofstream& input_log_ofs,
	ofstream& output_log_ofs
	)
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
	input_log_ofs << endl << tmpTimeStr << "... input of FixOneEndUnmapped starts ......" << endl << endl; 

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

	alignInfoInputQueue->initiateWith1stAlignInfo(
		line1_1st, line2_1st, line3_1st, line4_1st, line5_1st, 
		line6_1st, line7_1st, line8_1st, line9_1st, line10_1st, input_log_ofs);
	tmpBatchIndex_input ++;
	bool input_stage_end_bool = false;
	bool output_stage_end_bool = false;
	while(1)
	{
		if(!input_stage_end_bool)
		{
			for(int tmp = 0; tmp < inputTimePerc_fixOneEndUnmapped; tmp++)
			{
				if(inputRecord_ifs.eof())
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					input_log_ofs << endl << tmpTimeStr << "... end of input of FixOneEndUnmapped ......" << endl << endl;  
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
				string line1, line2, line3, line4, line5, line6,
					line7, line8, line9, line10;
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
				alignInfoInputQueue->getAlignInfoFromInputFile(line1, line2,
					line3, line4, line5, line6, line7, line8, line9, line10,
					inputReadNumInBatchArray_fixOneEndUnmapped, input_log_ofs, tmpBatchIndex_input);
			}
		}
		if(!output_stage_end_bool)
		{
			for(int tmp = 0; tmp < outputTimePerc_fixOneEndUnmapped; tmp++)
			{
				if(fixOneEndUnmappedResultQueue->atLeast3Node())
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					tmpBatchIndex_output ++;
					output_log_ofs << endl << tmpTimeStr << "... output regular array 2 file of fixOneEndUnmapped starts ......" << endl;  						
					output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
					fixOneEndUnmappedResultQueue->outputFrontResultArray(	
						OutputSamFile_oneEndMapped_ofs,
						tmpAlignIncompletePair_ofs,
						OutputSamFile_oneEndMapped_unpair_ofs,
						OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs						
						);
					fixOneEndUnmappedResultQueue->popFromResultQueue();
				}
				else if(fixOneEndUnmappedResultQueue->only2Node())
				{
					if(endOfProcessing_bool)
					{
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						tmpBatchIndex_output ++;
						output_log_ofs << endl << tmpTimeStr << "... output last array 2 file of fixOneEndUnmapped starts ......" << endl;  						
						output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
						fixOneEndUnmappedResultQueue->outputFrontResultArray(
							OutputSamFile_oneEndMapped_ofs, tmpAlignIncompletePair_ofs,
							OutputSamFile_oneEndMapped_unpair_ofs,
							OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs);
						fixOneEndUnmappedResultQueue->popFromResultQueue();
					}
					else
						continue;
				}
				else if(fixOneEndUnmappedResultQueue->resultQueueEmpty())
				{
					if(endOfProcessing_bool)
					{
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						output_log_ofs << endl << tmpTimeStr << "... end of output of fixOneEndUnmapped ......" << endl << endl;  
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

void process_function_fixOneEndUnmapped(
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array* tmpResultFixOneEndUnmappedArray,
	int tmpBatchArraySize,
	int threadNumForProcess,
	
	bool fasta_or_fastq_bool,
	Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	bool Do_extendHeadTail_fixOneEndUnmapped,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	int spliceJunctionDistanceMax,
	bool checkQualSeqForReadSegSeq,
	
	ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,

	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
	
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool)
{
		omp_set_num_threads(threadNumForProcess);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
		{	
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				//tmpRecordNum_oneEndUnmapped ++;
				int threadNO = omp_get_thread_num();
				string tmpReadNameAlignNum_1 = alignInfoInputQueue->returnFrontNodeReadNameAlignNum_1(tmpOpenMP);
				string tmpReadSeq_1 = alignInfoInputQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
				string tmpReadQualSeq_1 = alignInfoInputQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
				string tmpReadNameAlignNum_2 = alignInfoInputQueue->returnFrontNodeReadNameAlignNum_2(tmpOpenMP);
				string tmpReadSeq_2 = alignInfoInputQueue->returnFrontNodeReadSeq_2(tmpOpenMP);
				string tmpReadQualSeq_2 = alignInfoInputQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
				string tmpAlignInfo_Nor1 = alignInfoInputQueue->returnFrontNodeAlignInfo_Nor1(tmpOpenMP);
				string tmpAlignInfo_Rcm1 = alignInfoInputQueue->returnFrontNodeAlignInfo_Rcm1(tmpOpenMP);
				string tmpAlignInfo_Nor2 = alignInfoInputQueue->returnFrontNodeAlignInfo_Nor2(tmpOpenMP);
				string tmpAlignInfo_Rcm2 = alignInfoInputQueue->returnFrontNodeAlignInfo_Rcm2(tmpOpenMP);

				int multiMapSeg_maxLength = 0;

				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixOneEndUnmapped_getline(
					tmpReadNameAlignNum_1, tmpReadSeq_1, tmpReadQualSeq_1,
					tmpReadNameAlignNum_2, tmpReadSeq_2, tmpReadQualSeq_2,
					tmpAlignInfo_Nor1, tmpAlignInfo_Rcm1, //line9StrVec[tmpOpenMP],
					tmpAlignInfo_Nor2, tmpAlignInfo_Rcm2, peReadInfo, 
					indexInfo, fasta_or_fastq_bool, SE_or_PE_bool, multiMapSeg_maxLength);		

				FixOneEndUnmappedInfo* fixOneEndUnmappedInfo = new FixOneEndUnmappedInfo();
				fixOneEndUnmappedInfo->fixOneEndUnmapped(peReadInfo, peAlignInfo,
					secondLevelChrom,
					secondLevelSa,
					secondLevelLcpCompress,
					secondLevelChildTab,
					secondLevelDetChild,
					indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
					checkQualSeqForReadSegSeq);

				//#ifdef PERSONALIZED_CHR_SEQ
				if(do_segMap2snpMer_bool)
				{
					fixOneEndUnmappedInfo->fixOneEndUnmapped_includeSNPseqMap(peReadInfo, peAlignInfo,
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2,
						checkQualSeqForReadSegSeq,
						sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP,
						verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq);	
				}
				//#endif

				peAlignInfo->alignmentFilter_fixOneEndUnmapped_SJpenalty(peReadInfo.returnReadSeqLength_1(),
					peReadInfo.returnReadSeqLength_2());

				bool pairExistsBool = peAlignInfo->finalPairExistsBool();
				bool allAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();
				bool allUnpairedAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
				bool unique_bool = peAlignInfo->checkUniqueOrMulti();

				statsInfo->increNum_fixUnpaired(threadNO, pairExistsBool, allAlignmentCompleteBool, 
					allUnpairedAlignmentCompleteBool, unique_bool);

				if(pairExistsBool && allAlignmentCompleteBool) // some pair exists, all completed, print out paired SAM info
				{
					bool completeAlign_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
					if(completeAlign_lowScore_bool)
					{
						statsInfo->increLowScoreComplete_fixUnpair(threadNO, unique_bool);				
						tmpResultFixOneEndUnmappedArray->insert_PeAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair(
								peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
					}
					else
						tmpResultFixOneEndUnmappedArray->insert_PeAlignSamVec_fixUnpair(
							peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
				}
				else if(pairExistsBool && (!allAlignmentCompleteBool)) // pair exists, incomplete
					tmpResultFixOneEndUnmappedArray->insert_PeAlignInfoVec_fixUnpair(peAlignInfo->getTmpAlignInfoForFinalPair(
						peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(), peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(), 
						peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(), fasta_or_fastq_bool), tmpOpenMP);
				else			
				{
					if(!outputUnpairedSAM_bool)
						tmpResultFixOneEndUnmappedArray->insert_PeAlignSamVec_unpair_fixUnpair(
							peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
					else
					{
						bool partialAlignmentExistsOrNot_bool = peAlignInfo->partialAlignmentExists();
						if(partialAlignmentExistsOrNot_bool)
							tmpResultFixOneEndUnmappedArray->insert_PeAlignInfoVec_fixUnpair(peAlignInfo->getTmpAlignInfo(
								peReadInfo.returnReadName_1(), peReadInfo.returnReadName_2(), peReadInfo.returnReadSeq_1(), peReadInfo.returnReadSeq_2(),
								peReadInfo.returnReadQual_1(), peReadInfo.returnReadQual_2(), fasta_or_fastq_bool), tmpOpenMP);

						else
						{
							peAlignInfo->chooseBestAlignment_final_PEasSE();
							if(outputUnpairedSAM_bothEndsUniqueMappedOnly_bool)
								tmpResultFixOneEndUnmappedArray->insert_PeAlignSamVec_unpair_fixUnpair(
									peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_outputUniqueUnpairedBothEndsMappedOnly(
										peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
							else
								tmpResultFixOneEndUnmappedArray->insert_PeAlignSamVec_unpair_fixUnpair(
									peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot_allCases(
										peReadInfo, fasta_or_fastq_bool), tmpOpenMP);
						}
					}
				}
				delete fixOneEndUnmappedInfo;
				peAlignInfo->memoryFree();
				delete peAlignInfo;
		}
}

void process_stage_fixOneEndUnmapped_separateThreadForIO(
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,
	int threadNumForProcess,	
	bool fasta_or_fastq_bool,
	Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	bool Do_extendHeadTail_fixOneEndUnmapped,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	int spliceJunctionDistanceMax,
	bool checkQualSeqForReadSegSeq,
	ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,

	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool)
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
			mapping_log_ofs << "at least 3 node in fixOneEndUnmapped" << endl;
			tmpBatchArraySize = alignInfoInputQueue->returnFrontNodeSize();
		}
		else if(alignInfoInputQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << "... only 2 node -- end of file in fixOneEndUnmapped" << endl;
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
				mapping_log_ofs << "... inputQueueEmpty -- end of file in  fixOneEndUnmapped" << endl;
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
			<< endl << tmpTimeStr << "... mapping a batch of reads starts in fixOneEndUnmapped......" << endl;  
		mapping_log_ofs << "In fixOneEndUnmapped tmpBatchIndex: " << tmpBatchIndex << endl;
		mapping_log_ofs << "In fixOneEndUnmapped threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "In fixOneEndUnmapped tmpBatchArraySize: " << tmpBatchArraySize << endl;			
		Result_FixOneEndUnmapped_Array* tmpResultFixOneEndUnmappedArray = new Result_FixOneEndUnmapped_Array(tmpBatchArraySize);
		
		process_function_fixOneEndUnmapped(
			alignInfoInputQueue,
			tmpResultFixOneEndUnmappedArray,
			tmpBatchArraySize,
			threadNumForProcess,
			fasta_or_fastq_bool,
			statsInfo,
			secondLevelChrom,
			secondLevelSa,
			secondLevelLcpCompress,
			secondLevelChildTab,
			secondLevelDetChild,
			indexInfo, 
			Do_extendHeadTail_fixOneEndUnmapped,
			annotation_provided_bool, 
			Do_annotation_only_bool, 
			annotationInfo,
			spliceJunctionDistanceMax,
			checkQualSeqForReadSegSeq,
			mapping_log_ofs,
			SE_or_PE_bool,

			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

			sa_SNP, 
			lcpCompress_SNP,
			childTab_SNP, 
			chrom_SNP,
			verifyChild_SNP, 
			indexInfo_SNP,
			SNPlocInSyntheticSNPseq,
			do_segMap2snpMer_bool);

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << endl << tmpTimeStr << "... mapping a batch of input alignInfo ends in fixOneEndUnmapped......" 
			<< endl << "*******************************************************************************" << endl << endl; 		
		fixOneEndUnmappedResultQueue->pushBack2ResultArrayQueue(tmpResultFixOneEndUnmappedArray);
		alignInfoInputQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}

void process_stage_fixOneEndUnmapped_separateThreadForIO(
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,
	int threadNumForProcess,	
	bool fasta_or_fastq_bool,
	Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	bool Do_extendHeadTail_fixOneEndUnmapped,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	int spliceJunctionDistanceMax,
	bool checkQualSeqForReadSegSeq,
	ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,

	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	// provided SNP
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool,
	// learned SNP
	unsigned int* sa_SNP_learned,
	BYTE* lcpCompress_SNP_learned,
	unsigned int* childTab_SNP_learned,
	char* chrom_SNP_learned,
	BYTE* verifyChild_SNP_learned,
	Index_Info* indexInfo_SNP_learned,
	int SNPlocInSyntheticSNPseq_learned,
	bool do_segMap2snpMer_bool_learned,
	bool snpLearned_success_bool)
{
	if(do_segMap2snpMer_bool && (!do_segMap2snpMer_bool_learned))
	{
		process_stage_fixOneEndUnmapped_separateThreadForIO(alignInfoInputQueue, fixOneEndUnmappedResultQueue,
			endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo,
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, 
			annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNP
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);
	}
	else if((!do_segMap2snpMer_bool) && do_segMap2snpMer_bool_learned)
	{
		process_stage_fixOneEndUnmapped_separateThreadForIO(alignInfoInputQueue, fixOneEndUnmappedResultQueue,
			endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo,
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, 
			annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// learned SNP
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned,
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, snpLearned_success_bool);	
	}
	else if((!do_segMap2snpMer_bool) && (!do_segMap2snpMer_bool_learned))
	{
		process_stage_fixOneEndUnmapped_separateThreadForIO(alignInfoInputQueue, fixOneEndUnmappedResultQueue,
			endOfFile_bool, endOfProcessing_bool, threadNumForProcess, fasta_or_fastq_bool, statsInfo,
			secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild,
			indexInfo, Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, 
			annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq, mapping_log_ofs, SE_or_PE_bool,
			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNP
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

void io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(
	ifstream& inputRecord_ifs,
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue,
	int inputReadNumInBatchArray_fixOneEndUnmapped,
	ofstream& log_ofs,
	ofstream& OutputSamFile_oneEndMapped_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& OutputSamFile_oneEndMapped_unpair_ofs,
	ofstream& OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,	
	ofstream& input_log_ofs,
	ofstream& output_log_ofs,
	int threadNumForProcess,	
	bool fasta_or_fastq_bool,
	Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	bool Do_extendHeadTail_fixOneEndUnmapped,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	int spliceJunctionDistanceMax,
	bool checkQualSeqForReadSegSeq,
	ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,

	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool)
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
		for(int tmp = 0; tmp < inputReadNumInBatchArray_fixOneEndUnmapped; tmp++)
		{
			if(inputRecord_ifs.eof())
			{
				nowtime = time(NULL);
				local = localtime(&nowtime);
				tmpTimeStr = asctime(local);
				tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
				input_log_ofs << endl << tmpTimeStr << "... end of input of FixOneEndUnmapped ......" << endl << endl;
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
				alignInfoInputQueue->initiateWith1stAlignInfo(line1, line2, line3, 
					line4, line5,  line6, line7, line8, line9, line10, input_log_ofs);
			else
				alignInfoInputQueue->pushBack2currentFrontAlignInfoArray(line1, line2, line3, 
					line4, line5,  line6, line7, line8, line9, line10);
		}
		if(alignInfoInputQueue->inputQueueEmpty())
		{
			input_log_ofs << "inputQueueEmpty ! end of fixHeadTail !" << endl;
			//cout << "readQueueEmpty ! end of phase 1 !" << endl;
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
		// nowtime = time(NULL);
		// local = localtime(&nowtime);
		// tmpTimeStr = asctime(local);
		// tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << "***********************************************************" << endl;
		mapping_log_ofs << tmpTimeStr << "... start to map reads ......" << endl;
		mapping_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		mapping_log_ofs << "ThreadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "BatchArraySize: " << tmpBatchArraySize << endl;	
		Result_FixOneEndUnmapped_Array* tmpResultFixOneEndUnmappedArray = new Result_FixOneEndUnmapped_Array(tmpBatchArraySize);
		process_function_fixOneEndUnmapped(
			alignInfoInputQueue,
			tmpResultFixOneEndUnmappedArray,
			tmpBatchArraySize,
			threadNumForProcess,
			fasta_or_fastq_bool,
			statsInfo,
			secondLevelChrom,
			secondLevelSa,
			secondLevelLcpCompress,
			secondLevelChildTab,
			secondLevelDetChild,
			indexInfo, 
			Do_extendHeadTail_fixOneEndUnmapped,
			annotation_provided_bool, 
			Do_annotation_only_bool, 
			annotationInfo,
			spliceJunctionDistanceMax,
			checkQualSeqForReadSegSeq,
			mapping_log_ofs,
			SE_or_PE_bool,

			outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

			sa_SNP, 
			lcpCompress_SNP,
			childTab_SNP, 
			chrom_SNP,
			verifyChild_SNP, 
			indexInfo_SNP,
			SNPlocInSyntheticSNPseq,
			do_segMap2snpMer_bool);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << tmpTimeStr << "... end of mapping reads ......" << endl << endl;	
		fixOneEndUnmappedResultQueue->pushBack2ResultArrayQueue(tmpResultFixOneEndUnmappedArray);
		alignInfoInputQueue->popFromReadQueue();
		// output stage
		// nowtime = time(NULL);
		// local = localtime(&nowtime);
		// tmpTimeStr = asctime(local);
		// tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		output_log_ofs << "***********************************************************" << endl;
		output_log_ofs << tmpTimeStr << "... start to output alignment ......" << endl;
		output_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		fixOneEndUnmappedResultQueue->outputFrontResultArray(	
			OutputSamFile_oneEndMapped_ofs,
			tmpAlignIncompletePair_ofs,
			OutputSamFile_oneEndMapped_unpair_ofs,
			OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		output_log_ofs << tmpTimeStr << "... end of outputting alignment ......" << endl << endl;
		fixOneEndUnmappedResultQueue->popFromResultQueue();
	}
}

void io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(
	ifstream& inputRecord_ifs,
	AlignInfoInput_Array_Queue* alignInfoInputQueue,
	Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue,
	int inputReadNumInBatchArray_fixOneEndUnmapped,
	ofstream& log_ofs,
	ofstream& OutputSamFile_oneEndMapped_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& OutputSamFile_oneEndMapped_unpair_ofs,
	ofstream& OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,	
	ofstream& input_log_ofs,
	ofstream& output_log_ofs,
	int threadNumForProcess,	
	bool fasta_or_fastq_bool,
	Stats_Info* statsInfo,
	vector<char*>& secondLevelChrom,
	vector<unsigned int*>& secondLevelSa,
	vector<BYTE*>& secondLevelLcpCompress,
	vector<unsigned int*>& secondLevelChildTab,
	vector<BYTE*>& secondLevelDetChild,
	Index_Info* indexInfo, 
	bool Do_extendHeadTail_fixOneEndUnmapped,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,
	int spliceJunctionDistanceMax,
	bool checkQualSeqForReadSegSeq,
	ofstream& mapping_log_ofs,
	bool SE_or_PE_bool,

	bool outputUnpairedSAM_bool, bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,

	// provided SNP
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool,
	// learned SNP
	unsigned int* sa_SNP_learned, 
	BYTE* lcpCompress_SNP_learned,
	unsigned int* childTab_SNP_learned, 
	char* chrom_SNP_learned,
	BYTE* verifyChild_SNP_learned, 
	Index_Info* indexInfo_SNP_learned,
	int SNPlocInSyntheticSNPseq_learned,
	bool do_segMap2snpMer_bool_learned,
	bool snpLearned_success_bool)
{
	if(do_segMap2snpMer_bool && (!do_segMap2snpMer_bool_learned))
	{
		io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(inputRecord_ifs, alignInfoInputQueue, fixOneEndUnmappedResultQueue, 
			inputReadNumInBatchArray_fixOneEndUnmapped, log_ofs, OutputSamFile_oneEndMapped_ofs, tmpAlignIncompletePair_ofs,
			OutputSamFile_oneEndMapped_unpair_ofs, OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs, input_log_ofs,
			output_log_ofs, threadNumForProcess, fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa,
			secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq,
			mapping_log_ofs, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNP
			sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
			indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool);
	}
	else if((!do_segMap2snpMer_bool) && do_segMap2snpMer_bool_learned)
	{
		io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(inputRecord_ifs, alignInfoInputQueue, fixOneEndUnmappedResultQueue, 
			inputReadNumInBatchArray_fixOneEndUnmapped, log_ofs, OutputSamFile_oneEndMapped_ofs, tmpAlignIncompletePair_ofs,
			OutputSamFile_oneEndMapped_unpair_ofs, OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs, input_log_ofs,
			output_log_ofs, threadNumForProcess, fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa,
			secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq,
			mapping_log_ofs, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// learned SNP
			sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned,
			indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, snpLearned_success_bool);
	}
	else if((!do_segMap2snpMer_bool) && (!do_segMap2snpMer_bool_learned))
	{
		io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(inputRecord_ifs, alignInfoInputQueue, fixOneEndUnmappedResultQueue, 
			inputReadNumInBatchArray_fixOneEndUnmapped, log_ofs, OutputSamFile_oneEndMapped_ofs, tmpAlignIncompletePair_ofs,
			OutputSamFile_oneEndMapped_unpair_ofs, OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs, input_log_ofs,
			output_log_ofs, threadNumForProcess, fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa,
			secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, Do_extendHeadTail_fixOneEndUnmapped,
			annotation_provided_bool, Do_annotation_only_bool, annotationInfo, spliceJunctionDistanceMax, checkQualSeqForReadSegSeq,
			mapping_log_ofs, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
			// provided SNP
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

#endif