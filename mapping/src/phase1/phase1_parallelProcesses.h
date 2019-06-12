// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef PHASE1_PARALLELPROCESSES_H
#define PHASE1_PARALLELPROCESSES_H

using namespace std;

//Truncates read name to just the part before first whitespace (for Cufflinks compatibility)
//Performed by default, command line option --keepFullReadName overrides this
void truncateReadName(string &readName) {
    int pos;
    if ((pos = readName.find(" ")) != string::npos) {
        readName = readName.substr(0, pos);
    }
}

void process_function_phase1(
	Read_Array_Queue* readArrayQueue, Result_Array* tmpResultArray,
	int tmpBatchArraySize, int threadNumForProcess, unsigned int* sa, 
	BYTE* lcpCompress, unsigned int* childTab, char* chrom,
	BYTE* verifyChild, Index_Info* indexInfo, int* preIndexMapLengthArray, 
	unsigned int* preIndexIntervalStartArray, unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec, bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1, bool annotation_provided_bool, 
	bool Do_annotation_only_bool, Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only, bool Do_Phase1_Only,	
	Stats_Info* statsInfo, bool fasta_or_fastq_bool, ofstream& mapping_log_ofs, 
	bool checkQualSeqForReadSegSeq, bool SE_or_PE_bool, unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, char* chrom_SNP, 
	BYTE* verifyChild_SNP, Index_Info* indexInfo_SNP, int SNPlocInSyntheticSNPseq, 
	bool do_segMap2snpMer_bool,
	LearnedCandiSNPhash_Info_Vec& mismatchHashInfoVec,
	bool collectMismatchfromAlignment_bool)
{
		omp_set_num_threads(threadNumForProcess);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
		{	
			int threadNO = omp_get_thread_num();
			
			string tmpReadName_1 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_1(tmpOpenMP);
			string tmpReadSeq_1 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
			string tmpReadName_2 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_2(tmpOpenMP);
			string tmpReadSeq_2 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_2(tmpOpenMP);	
			string tmpReadQualSeq_1 = "";
			string tmpReadQualSeq_2 = "";
			//cout << "tmpReadName_1: " << tmpReadName_1 << endl;
			//cout << "tmpReadName_2: " << tmpReadName_2 << endl;
			if(!fasta_or_fastq_bool)
			{	
				tmpReadQualSeq_1 = readArrayQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
				tmpReadQualSeq_2 = readArrayQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
			}

 			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(tmpReadName_1, tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2, fasta_or_fastq_bool, SE_or_PE_bool);	
    		FixPhase1Info fixPhase1Info;
    		//cout << "start to do fixPhase1_segInfo ..." << endl;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);						

			//#ifdef PERSONALIZED_CHR_SEQ
			if(do_segMap2snpMer_bool)
			{
				fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
					childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
				fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
					childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			}
			//#endif
			//cout << "start to do fixPhase1_pathInfo ..." << endl;
			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);
			//cout << "start to do fixPhase1_gapInfo ..." << endl;
			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);
			//cout << "start to do initiatePeAlignInfo ..." << endl;
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo(fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);		
			//cout << "start to do candi alignment selection ..." << endl;
			if(Do_Phase1_Only)
			{
				//if(SE_or_PE_bool)
				//	peAlignInfo.chooseBestAlignment_final_SE();
				//else
					peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{	
				//if(SE_or_PE_bool)
				//	peAlignInfo.alignmentFilter_fixPhase1_SE(readInfo.returnReadLength_SE());
				//else
					peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
						readInfo.returnReadSeqLength_2());
			}
			int overlapLength_top2candiAlignment = peAlignInfo.getMaxOverlapLengthInAlignmentPair(
				readInfo.returnReadSeqLength_1(), readInfo.returnReadSeqLength_2());
			//cout << endl << endl << "peAlignInfo: \n" << peAlignInfo.returnPeAlignInfoStr() << endl << endl;
			//cout << "start to do mismatch collection ..." << endl;
			int tmpReadLength_end1 = readInfo.returnReadSeqLength_1();
			int tmpReadLength_end2 = readInfo.returnReadSeqLength_2();
			if(collectMismatchfromAlignment_bool)
				peAlignInfo.collectMismatch(mismatchHashInfoVec, threadNO, indexInfo, tmpReadLength_end1, tmpReadLength_end2);
			//cout << "start to do output_phase1" << endl;
			peAlignInfo.output_phase1(outputDirectlyBool_Phase1Only, Do_Phase1_Only, tmpResultArray, repeatRegionInfoVec, 
				readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo, fasta_or_fastq_bool, overlapLength_top2candiAlignment);
			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();
		}
}

void io_stage_phase1_separateThreadForIO(ifstream& input_ifs_1,
	ifstream& input_ifs_2,
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,
	int inputReadNumInBatchArray_phase1,
	int inputTimePerc_phase1,
	int outputTimePerc_phase1,
	ofstream& log_ofs,
	InputReadPreProcess* readPreProcessInfo,
	ofstream& tmpAlignCompleteRead_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& tmpAlignOneEndUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
	ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
	ofstream& tmpAlignIncompletePair_SAM_ofs,
	ofstream& input_log_ofs,
	ofstream& output_log_ofs,
	bool fasta_or_fastq_bool,
	bool SE_or_PE_bool,
	int& tmp_total_read_num,
	bool keepFullReadName_bool
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
	input_log_ofs << endl << tmpTimeStr << "... input of Phase1 starts ......" << endl << endl;  

	string readNameStr_1, readNameStr_2, readSeqStr_1, readSeqStr_2;	
	getline(input_ifs_1, readNameStr_1);
	getline(input_ifs_1, readSeqStr_1);
	getline(input_ifs_2, readNameStr_2);
	getline(input_ifs_2, readSeqStr_2);
	
	if (!keepFullReadName_bool) {
	    truncateReadName(readNameStr_1);
	    truncateReadName(readNameStr_2);
	}
	
	if(fasta_or_fastq_bool)
	{	
		readArrayQueue->initiateWith1stRead(readNameStr_1, readNameStr_2, readSeqStr_1, readSeqStr_2, input_log_ofs);
		tmp_total_read_num ++;
		tmpBatchIndex_input ++;
	}
	else
	{
		string readCommentStr_1, readCommentStr_2, readQualSeq_1, readQualSeq_2;
		getline(input_ifs_1, readCommentStr_1);
		getline(input_ifs_1, readQualSeq_1);
		getline(input_ifs_2, readCommentStr_2);
		getline(input_ifs_2, readQualSeq_2);
		readArrayQueue->initiateWith1stRead_fq(readNameStr_1, readNameStr_2, readSeqStr_1, readSeqStr_2, readQualSeq_1, readQualSeq_2, input_log_ofs);		
		tmp_total_read_num ++;
		tmpBatchIndex_input ++;
	}

	bool input_stage_end_bool = false;
	bool output_stage_end_bool = false;
	while(1)
	{
		if(!input_stage_end_bool)
		{
			// read fa/fq file and insert the readArray into readArrayList
			for(int tmp = 0; tmp < inputTimePerc_phase1; tmp++)
			{
				if((input_ifs_1.eof())||(input_ifs_2.eof()))
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					input_log_ofs << endl << tmpTimeStr << "... end of input ......" << endl << endl;  
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadName_1, tmpReadName_2;
				getline(input_ifs_1, tmpReadName_1);
				getline(input_ifs_2, tmpReadName_2); 
				
				if (!keepFullReadName_bool) {
	                truncateReadName(tmpReadName_1);
	                truncateReadName(tmpReadName_2);
	            }
				
				if((input_ifs_1.eof())||(input_ifs_2.eof())||(tmpReadName_1 == "")||(tmpReadName_2 == ""))
				{
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadSeq_1, tmpReadSeq_2;
				getline(input_ifs_1, tmpReadSeq_1);
				getline(input_ifs_2, tmpReadSeq_2);

				if(fasta_or_fastq_bool)
				{
					readArrayQueue->getSeqFromInputFile(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
					tmp_total_read_num ++;
				}
				else
				{
					string tmpCommentStr_1, tmpCommentStr_2, tmpReadQualSeq_1, tmpReadQualSeq_2;
					getline(input_ifs_1, tmpCommentStr_1);
					getline(input_ifs_1, tmpReadQualSeq_1);
					getline(input_ifs_2, tmpCommentStr_2);
					getline(input_ifs_2, tmpReadQualSeq_2);
					readArrayQueue->getSeqFromInputFile_fq(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,						
						tmpReadQualSeq_1,
						tmpReadQualSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
					tmp_total_read_num ++;
				}
			}
		}
		if(!output_stage_end_bool)
		{	
			for(int tmp = 0; tmp < outputTimePerc_phase1; tmp++)
			{		
				if(resultArrayQueue->atLeast3Node())
				{
					nowtime = time(NULL);
					local = localtime(&nowtime);
					tmpTimeStr = asctime(local);
					tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
					tmpBatchIndex_output ++;
					output_log_ofs << endl << tmpTimeStr << "... output regular array to file starts ......" << endl;  					
					output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
					resultArrayQueue->outputFrontResultArray(
						tmpAlignCompleteRead_ofs,
						tmpAlignIncompletePair_ofs,
						tmpAlignOneEndUnmapped_ofs,
						tmpAlignBothEndsUnmapped_ofs,
						tmpAlignBothEndsUnmapped_lowScore_ofs,
						tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
						tmpAlignIncompletePair_SAM_ofs);
					resultArrayQueue->popFromResultQueue();
				}
				else if(resultArrayQueue->only2Node())
				{
					if(endOfProcessing_bool)
					{	
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						tmpBatchIndex_output ++;
						output_log_ofs << endl << tmpTimeStr << "... output last array to file starts ......" << endl;
						output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
						resultArrayQueue->outputFrontResultArray(
							tmpAlignCompleteRead_ofs,
							tmpAlignIncompletePair_ofs,
							tmpAlignOneEndUnmapped_ofs,
							tmpAlignBothEndsUnmapped_ofs,
							tmpAlignBothEndsUnmapped_lowScore_ofs,
							tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
							tmpAlignIncompletePair_SAM_ofs);
						resultArrayQueue->popFromResultQueue();
					}
					else
						continue;
				}
				else if(resultArrayQueue->resultQueueEmpty())
				{
					if(endOfProcessing_bool)
					{
						nowtime = time(NULL);
						local = localtime(&nowtime);
						tmpTimeStr = asctime(local);
						tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
						output_log_ofs << endl << tmpTimeStr << "... end of output ......" << endl << endl;  
						output_stage_end_bool = true;
						break;	
					}
				}
				else
				{
					//cout << "other cases in output stage" << endl;
					continue;
				}
			}
		}
		if(input_stage_end_bool && output_stage_end_bool)
			break;
	}
	log_ofs << "end of input/output " << endl;
}

void process_stage_phase1_separateThreadForIO(
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,	
	int threadNumForProcess,
	unsigned int* sa, 
	BYTE* lcpCompress,
	unsigned int* childTab, 
	char* chrom,
	BYTE* verifyChild, 
	Index_Info* indexInfo,
	int* preIndexMapLengthArray, 
	unsigned int* preIndexIntervalStartArray,
	unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec,
	bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only,
	bool Do_Phase1_Only,	
	Stats_Info* statsInfo, 
	bool fasta_or_fastq_bool,
	ofstream& mapping_log_ofs,//, vector<ofstream*>& mapping_log_ofs_vec
	bool checkQualSeqForReadSegSeq,
	bool SE_or_PE_bool,
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq,
	bool do_segMap2snpMer_bool,
	LearnedCandiSNPhash_Info_Vec& mismatchHashInfoVec,
	bool collectMismatchfromAlignment_bool)
{
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;	
	int tmpBatchIndex = 0;
	while(1)
	{
		int tmpBatchArraySize = 0;
		if(readArrayQueue->atLeast3Node())
		{
			mapping_log_ofs << endl << "at least 3 node" << endl;
			tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
		}
		else if(readArrayQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... only 2 node -- end of file" << endl;
				tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
			}
			else
			{
				//mapping_log_ofs << "... only 2 node -- not end of file" << endl;
				continue;
			}
		}
		else if(readArrayQueue->readQueueEmpty())
		{
			//mapping_log_ofs << "readQueueEmpty" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... readQueueEmpty -- end of file" << endl;
				break;
			}
			else
			{
				//mapping_log_ofs << "... readQueueEmpty -- not end of file" << endl;
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
		mapping_log_ofs << "*******************************************************************************" 
			<< endl << tmpTimeStr << "... mapping a batch of reads starts ......" << endl;  
		tmpBatchIndex ++;
		mapping_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl;
		mapping_log_ofs << "threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "tmpBatchArraySize: " << tmpBatchArraySize << endl;			

		Result_Array* tmpResultArray = new Result_Array(tmpBatchArraySize);
		process_function_phase1(
			readArrayQueue, tmpResultArray, tmpBatchArraySize, threadNumForProcess, 
			sa, lcpCompress, childTab, chrom, verifyChild, indexInfo,
			preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
			repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
			Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only,
			Do_Phase1_Only,	statsInfo, fasta_or_fastq_bool, mapping_log_ofs,
			checkQualSeqForReadSegSeq, SE_or_PE_bool, sa_SNP, lcpCompress_SNP, childTab_SNP, 
			chrom_SNP, verifyChild_SNP, indexInfo_SNP, SNPlocInSyntheticSNPseq, do_segMap2snpMer_bool,
			mismatchHashInfoVec, collectMismatchfromAlignment_bool);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << endl << "*******************************************************************************"
			<< endl << tmpTimeStr << "... mapping a batch of reads ends ......" << endl << endl << endl; 		
		resultArrayQueue->pushBack2ResultArrayQueue(tmpResultArray);
		readArrayQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}

void io_process_phase1_allThreadsSharedByBothStage(
	ifstream& input_ifs_1, ifstream& input_ifs_2, Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue, int inputReadNumInBatchArray_phase1,
	ofstream& log_ofs, InputReadPreProcess* readPreProcessInfo,
	ofstream& tmpAlignCompleteRead_ofs, ofstream& tmpAlignIncompletePair_ofs,
	ofstream& tmpAlignOneEndUnmapped_ofs, ofstream& tmpAlignBothEndsUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
	ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
	ofstream& tmpAlignIncompletePair_SAM_ofs, ofstream& input_log_ofs,
	ofstream& output_log_ofs, bool fasta_or_fastq_bool,
	bool SE_or_PE_bool, int threadNumForProcess, unsigned int* sa, 
	BYTE* lcpCompress, unsigned int* childTab, char* chrom,
	BYTE* verifyChild, Index_Info* indexInfo,
	int* preIndexMapLengthArray, unsigned int* preIndexIntervalStartArray,
	unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec, bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1, bool annotation_provided_bool, 
	bool Do_annotation_only_bool, Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only, bool Do_Phase1_Only,	
	Stats_Info* statsInfo, ofstream& mapping_log_ofs,
	bool checkQualSeqForReadSegSeq, unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP, unsigned int* childTab_SNP, 
	char* chrom_SNP, BYTE* verifyChild_SNP, Index_Info* indexInfo_SNP, 
	int SNPlocInSyntheticSNPseq, bool do_segMap2snpMer_bool,
	LearnedCandiSNPhash_Info_Vec& mismatchHashInfoVec,
	bool collectMismatchfromAlignment_bool, int& tmp_total_read_num,
	bool keepFullReadName_bool)
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
		for(int tmp = 0; tmp < inputReadNumInBatchArray_phase1; tmp++)
		{
			if((input_ifs_1.eof())||(input_ifs_2.eof()))
			{
				nowtime = time(NULL);
				local = localtime(&nowtime);
				input_log_ofs << endl << "[" << asctime(local) << "... end of input ......" << endl << endl;  
				endOfFile_bool = true;
				break;
			}
			string tmpReadName_1, tmpReadName_2;
			getline(input_ifs_1, tmpReadName_1);
			getline(input_ifs_2, tmpReadName_2);
			
			if (!keepFullReadName_bool) {
	            truncateReadName(tmpReadName_1);
	            truncateReadName(tmpReadName_2);
	        }
			
			if((input_ifs_1.eof())||(input_ifs_2.eof()))
			{
				nowtime = time(NULL);
				local = localtime(&nowtime);
				input_log_ofs << endl << "[" << asctime(local) << "... end of input ......" << endl << endl;  
				endOfFile_bool = true;
				break;
			}			
			string tmpReadName_1_afterProcess = tmpReadName_1.substr(1);
			string tmpReadName_2_afterProcess = tmpReadName_2.substr(1);
			string tmpReadSeq_1, tmpReadSeq_2;
			getline(input_ifs_1, tmpReadSeq_1);
			getline(input_ifs_2, tmpReadSeq_2);
			string tmpReadSeq_1_afterProcess = readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1);
			string tmpReadSeq_2_afterProcess = readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2);
			if(fasta_or_fastq_bool)
			{
				tmp_total_read_num ++;
				if(tmp == 0)
					readArrayQueue->initiateWith1stRead(tmpReadName_1_afterProcess, tmpReadName_2_afterProcess, 
						tmpReadSeq_1_afterProcess, tmpReadSeq_2_afterProcess, input_log_ofs);					
				else
					readArrayQueue->pushBack2currentFrontReadArray(tmpReadName_1_afterProcess, tmpReadName_2_afterProcess,
						tmpReadSeq_1_afterProcess, tmpReadSeq_2_afterProcess);
			}
			else
			{
				tmp_total_read_num ++;
				string tmpCommentStr_1, tmpCommentStr_2, tmpReadQualSeq_1, tmpReadQualSeq_2;
				getline(input_ifs_1, tmpCommentStr_1);
				getline(input_ifs_1, tmpReadQualSeq_1);
				getline(input_ifs_2, tmpCommentStr_2);
				getline(input_ifs_2, tmpReadQualSeq_2);
				if(tmp == 0)
					readArrayQueue->initiateWith1stRead_fq(tmpReadName_1_afterProcess, tmpReadName_2_afterProcess, 
						tmpReadSeq_1_afterProcess, tmpReadSeq_2_afterProcess,
						tmpReadQualSeq_1, tmpReadQualSeq_2, input_log_ofs);					
				else
					readArrayQueue->pushBack2currentFrontReadArray_fq(tmpReadName_1_afterProcess, tmpReadName_2_afterProcess,
						 tmpReadSeq_1_afterProcess, tmpReadSeq_2_afterProcess, tmpReadQualSeq_1, tmpReadQualSeq_2);
			}
		}
		if(readArrayQueue->readQueueEmpty())
		{
			input_log_ofs << "readQueueEmpty ! end of phase 1 !" << endl;
			break;
		}
		int tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
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
		Result_Array* tmpResultArray = new Result_Array(tmpBatchArraySize);
		process_function_phase1(
			readArrayQueue,
			tmpResultArray,
			tmpBatchArraySize,
			threadNumForProcess,
			sa, 
			lcpCompress,
			childTab, 
			chrom,
			verifyChild, 
			indexInfo,
			preIndexMapLengthArray, 
			preIndexIntervalStartArray,
			preIndexIntervalEndArray,
			repeatRegionInfoVec,
			Do_cirRNA, 
			Do_extendHeadTail_phase1,
			annotation_provided_bool, 
			Do_annotation_only_bool, 
			annotationInfo,			
			outputDirectlyBool_Phase1Only,
			Do_Phase1_Only,	
			statsInfo, 
			fasta_or_fastq_bool,
			mapping_log_ofs,
			checkQualSeqForReadSegSeq,
			SE_or_PE_bool,
			sa_SNP, 
			lcpCompress_SNP,
			childTab_SNP, 
			chrom_SNP,
			verifyChild_SNP, 
			indexInfo_SNP,
			SNPlocInSyntheticSNPseq,
			do_segMap2snpMer_bool,
			mismatchHashInfoVec,
			collectMismatchfromAlignment_bool);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		mapping_log_ofs << tmpTimeStr << "... end of mapping reads ......" << endl << endl;	
		resultArrayQueue->pushBack2ResultArrayQueue(tmpResultArray);
		readArrayQueue->popFromReadQueue();
		// output stage 
		output_log_ofs << "***********************************************************" << endl;
		output_log_ofs << tmpTimeStr << "... start to output alignment ......" << endl;
		output_log_ofs << "BatchIndex :" << tmpBatchIndex << endl;
		resultArrayQueue->outputFrontResultArray(
			tmpAlignCompleteRead_ofs, tmpAlignIncompletePair_ofs, tmpAlignOneEndUnmapped_ofs,
			tmpAlignBothEndsUnmapped_ofs, tmpAlignBothEndsUnmapped_lowScore_ofs,
			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs, tmpAlignIncompletePair_SAM_ofs);
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		output_log_ofs << tmpTimeStr << "... end of outputting alignment ......" << endl << endl;
		resultArrayQueue->popFromResultQueue();
	}
}

#endif
