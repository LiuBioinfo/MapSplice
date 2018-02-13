// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
/*    
 *    mps_batch.cpp
 *	  MapSplice3
 *
 *    Copyright (C) 2016 University of Kentucky and
 *                       Xinan Liu
 *
 *    Authors: Xinan Liu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>
//#include "switch.h"
#include "general/extractUnmapAlignment2ReadFile.h"
//#include "phase1/arrayQueue.h"
#include "phase1/arrayQueue_phase1.h"
#include "phase2/arrayQueue_phase2.h"
#include "stats_info.h"
#include "constantDefinitions.h"
#include "general/option_info.h"
#include "general/read_block_test.h"
#include "general/bwtmap_info.h"
#include "general/DoubleAnchorScore.h"
#include "general/sbndm.h"
#include "general/otherFunc.h"
#include "general/index_info.h"
#include "general/enhanced_suffix_array_info.h"
#include "general/annotation_info.h"
#include "phase1/repeatRegion.h"
#include "general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "general/splice_info.h"
#include "general/fixGapRelationParameters.h"
#include "general/read_info.h"
#include "general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "general/fixDoubleAnchorNWDP_info.h"
#include "general/fixDoubleAnchorMatch_info.h"
#include "general/fixDoubleAnchorInsertion_info.h"
#include "general/fixDoubleAnchorDeletion_info.h"
#include "general/fixDoubleAnchorSplice_complicate_info.h"
#include "general/fixDoubleAnchorSplice_info.h"
#include "general/fixDoubleAnchorCirRNA_info.h"
#include "general/path_info.h"
#include "general/gap_info.h"
#include "otherProjects/incorporateGenomicVariants/general/syntheticSNPtransSeq_info.h"
#include "otherProjects/incorporateGenomicVariants/general/learnedCandiSNPhash_info_vec.h"
#include "general/align_info.h"
#include "general/peAlign_info.h"
#include "general/groupSeg_info.h"
#include "general/alignInferJunctionHash_info_vec.h"
#include "phase2/spliceJunctionHash_info.h"
#include "phase2/unmapEnd_info.h"
#include "phase2/unfixedHead.h"
#include "phase2/unfixedTail.h"
#include "phase2/incompleteLongHead.h"
#include "phase2/incompleteLongTail.h"
#include "phase2/sam2junc.h"
#include "fixHeadTail.h"
#include "phase2/fixOneEndUnmapped.h"
#include "fixPhase1.h"
#include "general/readSeqPreProcessing.h"
#include "general/headerSection_info.h"
#include "general/otherFunc2.h"
#include "phase1/phase1_parallelProcesses.h"
#include "phase2/fixOneEndUnmapped_parallelProcesses.h"
#include "phase2/fixHeadTail_parallelProcesses.h"
#include "otherProjects/indexing/utils/build_snpMerIndex.h"
#include "general/batch_manager/batch_manager_info.h"

using namespace std;  

int main(int argc, char**argv)
{
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;

	bool segMap2SNPmer_phase2_bool_learned = false;//true; // NOTE: to debug,
	bool snpLearned_success_bool = segMap2SNPmer_phase2_bool_learned;
	int SNPmerLength = 201;
	int SNPlocInSyntheticSNPseq = (SNPmerLength - 1)/2;
	int SNPlocInSyntheticSNPseq_learned = (SNPmerLength - 1)/2;

	int normalRecordNum_1stMapping = 500000;
	int normalRecordNum_fixOneEndUnmapped = 500000;
	int normalRecordNum_fixHeadTail = 500000;

    bool checkQualSeqForShortAnchorSeqToTargetMap = false;//true;
    bool checkQualSeqForReadSegSeq = false;//true;
	bool outputUnpairedSAM_bool = false;
	bool outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = false;
	//outputUnpairedSAM_bool = true;
	#ifdef OUTPUT_UNPAIRED_REGULAR
	outputUnpairedSAM_bool = true;
	outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = false;
	#endif	

	#ifdef OUTPUTUNIQUEUNPAIREDBOTHENDSMAPPEDSAM
	outputUnpairedSAM_bool = true;
	outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = true;
	#endif 

	#ifdef DETECT_CIRCULAR_RNA
	if(outputUnpairedSAM_bool)
	{
		cout << "for now, MPS3_circRNA only supports SE reads." << endl;
		cout << "MPS3_circRNA does not support 'outputUnpairedSAM' " << endl;
		exit(1);
	}
	#endif 

    /////////////////   get option from command line ////////////////////
    Option_Info* optionInfo = new Option_Info();
    optionInfo->getOpt_long(argc, argv);

    string outputFolder_path_batch = optionInfo->outputFolder_path;
    string mkdir_outputFolder_path_batch_command = "mkdir -p " + outputFolder_path_batch;
    system(mkdir_outputFolder_path_batch_command.c_str());
    string log_batch_file = outputFolder_path_batch + "/batch.log";
    string settings_batch_file = outputFolder_path_batch + "/batch.settings";
    string process_batch_file = outputFolder_path_batch + "/batch.process";
    string runtime_batch_file = outputFolder_path_batch + "/batch.runtime";
    ofstream log_batch_ofs(log_batch_file.c_str());
    ofstream settings_batch_ofs(settings_batch_file.c_str());
    ofstream process_batch_ofs(process_batch_file.c_str());
    ofstream runtime_batch_ofs(runtime_batch_file.c_str());

    bool backSplice_search_bool = optionInfo->return_backSplice_search_bool();
    bool fusion_search_bool = optionInfo->return_fusion_search_bool();
    bool reportJunc_bool = optionInfo->return_report_junc_bool();
    if(fusion_search_bool)
    {
		outputUnpairedSAM_bool = true;
		outputUnpairedSAM_bothEndsUniqueMappedOnly_bool = true;   	
    }
    #ifdef DETECT_CIRCULAR_RNA
    if(!backSplice_search_bool)
    {
    	cout << "backSplice_search_bool: false" << endl;
    	cout << "but DETECT_CIRCULAR_RNA is defined" << endl;
    	cout << "exiting ......" << endl;
    	exit(1);
    }
    #endif
    cout << "backSplice_search_bool: " << backSplice_search_bool << endl;
    log_batch_ofs << "backSplice_search_bool: " << backSplice_search_bool << endl;
    settings_batch_ofs << "backSplice_search_bool: " << backSplice_search_bool << endl;
    cout << "fusion_search_bool: " << fusion_search_bool << endl;
    log_batch_ofs << "fusion_search_bool: " << fusion_search_bool << endl;
    settings_batch_ofs << "fusion_search_bool: " << fusion_search_bool << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... MPS starts ......" << endl;  
	log_batch_ofs << endl << tmpTimeStr << "... MPS starts ......" << endl;
	runtime_batch_ofs << endl << tmpTimeStr << "... MPS starts ......" << endl; 

	bool sharingSplicingContext_bool = optionInfo->return_sharing_splicing_context_bool();
	bool separatingSplicingContext_bool = (!sharingSplicingContext_bool);
	bool sharingSNPcontext_bool = optionInfo->return_sharing_SNP_context_bool();
	bool separatingSNPcontext_bool = (!sharingSNPcontext_bool);
    bool sharedThreadForIO_bool = optionInfo->return_sharedThreadForIO_bool();
    bool separateThreadForIO_bool = (!sharedThreadForIO_bool);

    settings_batch_ofs << "sharingSplicingContext_bool: " << sharingSplicingContext_bool << endl;
    settings_batch_ofs << "sharingSNPcontext_bool: " << sharingSNPcontext_bool << endl;
    settings_batch_ofs << "sharedThreadForIO_bool: " << sharedThreadForIO_bool << endl;
    log_batch_ofs << "sharingSplicingContext_bool: " << sharingSplicingContext_bool << endl;
    log_batch_ofs << "sharingSNPcontext_bool: " << sharingSNPcontext_bool << endl;
    log_batch_ofs << "sharedThreadForIO_bool: " << sharedThreadForIO_bool << endl;  
	cout << "sharingSplicingContext_bool: " << sharingSplicingContext_bool << endl;
    cout << "sharingSNPcontext_bool: " << sharingSNPcontext_bool << endl;
    cout << "sharedThreadForIO_bool: " << sharedThreadForIO_bool << endl;  
    bool extractUnmapAlignment2ReadFile_bool = optionInfo->return_extractUnmapAlignment2ReadFile_bool();
    string inputSamFilePathStr = optionInfo->returnInputSamFilePath();
    if(extractUnmapAlignment2ReadFile_bool)
    	unmappedAlignemnt2ReadFile(inputSamFilePathStr);

    ////////////////////////////// check SE or PE reads  ///////////////////////////////
    bool SE_or_PE_bool = optionInfo->returnSEorPE_bool();
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
   	bool annotation_provided_bool = optionInfo->annotation_provided_bool;
   	bool Do_annotation_only_bool = false;//annotation_provided_bool;
	bool Do_Phase1_Only = optionInfo->Do_phase1_only_bool;
	bool outputAlignInfoAndSamForAllPairedAlignmentBool = false;
	//outputAlignInfoAndSamForAllPairedAlignmentBool = true;
	bool removeAllIntermediateFilesBool = false;
	//removeAllIntermediateFilesBool = true;
	bool DoSam2JuncBool = false;
	DoSam2JuncBool = true;
	bool load2ndLevelIndexBool = false;
	load2ndLevelIndexBool = true;
	bool load2ndLevelIndexBool_compressedSize = false;
	load2ndLevelIndexBool_compressedSize = true;
	bool DoMappingOnUnmapEndReadsBool = false;
	DoMappingOnUnmapEndReadsBool = true;
	bool DoMappingOnUnfixedHeadTailAlignmentBool = false;
	DoMappingOnUnfixedHeadTailAlignmentBool = true;
	bool outputDirectlyBool_Phase1Only = true;
	outputDirectlyBool_Phase1Only = false;
	bool Do_cirRNA = true;
	Do_cirRNA = false;

	//bool Do_extendHeadTail = true;
	//Do_extendHeadTail = false;

	bool Do_extendHeadTail_phase1 = true;
	//Do_extendHeadTail_phase1 = false;
	bool Do_extendHeadTail_fixOneEndUnmapped = true;
	//Do_extendHeadTail_fixOneEndUnmapped = false;
	bool Do_extendHeadTail_fixHeadTail = true;
	//Do_extendHeadTail_fixHeadTail = false;	
	bool Do_fixHeadTail_remapping = true;
	//Do_fixHeadTail_remapping = false;
	bool Do_fixHeadTail_greedyMapping = true;
	//Do_fixHeadTail_greedyMapping = false;
	bool Do_fixHeadTail_remappingAndTargetMapping = true;
	//Do_fixHeadTail_remappingAndTargetMapping = false;
	bool Do_fixHeadTail_remappingAgain = true;
	Do_fixHeadTail_remappingAgain = false;
	bool Do_fixHeadTail_extend2end_finalStep = true;
	//Do_fixHeadTail_extend2end_finalStep = false;
	bool Do_fixHeadTail_extend2end_fixIndel = true;
	//Do_fixHeadTail_extend2end_fixIndel = false;

	// Fix ME:
	// #ifdef DETECT_CIRCULAR_RNA
	// Do_fixHeadTail_remapping = false;
	// Do_fixHeadTail_greedyMapping = false;
	// Do_fixHeadTail_remappingAndTargetMapping = false;
	// #endif
	if(Do_Phase1_Only)
	{
		settings_batch_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		load2ndLevelIndexBool_compressedSize = false;
		DoMappingOnUnmapEndReadsBool = false;
		DoMappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		settings_batch_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;//false;
		load2ndLevelIndexBool = true;//false;
		load2ndLevelIndexBool_compressedSize = true;//false;
		DoMappingOnUnmapEndReadsBool = true;//false;
		DoMappingOnUnfixedHeadTailAlignmentBool = true;//false;
	}
	
	optionInfo->outputSwitchInfo(Do_Phase1_Only, outputAlignInfoAndSamForAllPairedAlignmentBool,
		removeAllIntermediateFilesBool, Do_cirRNA, outputDirectlyBool_Phase1Only, 
		normalRecordNum_1stMapping, normalRecordNum_fixOneEndUnmapped,
		normalRecordNum_fixHeadTail, Do_extendHeadTail_phase1, 
		Do_extendHeadTail_fixOneEndUnmapped, Do_extendHeadTail_fixHeadTail, 
		Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
		Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, settings_batch_ofs);
	//////////////////////////////////////////////        LOAD INDEX         ////////////////////////////////////////////////////////////////
    string InputReadFile = optionInfo->read_file_path_1;//read sample, exacted from fastq file every time
    if(SE_or_PE_bool)
    	InputReadFile = optionInfo->read_file_path_SE;
    string InputReadFile_PE = optionInfo->read_file_path_2;// another end read for pair-end reads
	int threads_num = optionInfo->threads_num;//atoi(threadsNumStr.c_str());
	if((threads_num == 1)&& separateThreadForIO_bool)
	{
		cout << "threads_num == 1, but separateThreadForIO_bool == true" << endl;
		exit(1);
	}
	#ifdef MAP_INFO
	threads_num = 1;
	#endif
	omp_set_num_threads(threads_num);

	bool InputAsFastq = (!(optionInfo->fasta_or_fastq_bool));
	bool fasta_or_fastq_bool = optionInfo->fasta_or_fastq_bool;
	cout << "fasta_or_fastq_bool: " << fasta_or_fastq_bool << endl;
	if(checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_batch_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start to load whole genome index ......" << endl; 
	log_batch_ofs << endl << tmpTimeStr << "... start to load whole genome index ......" << endl; 
	runtime_batch_ofs << endl << tmpTimeStr << "... start to load whole genome index ......" << endl; 

    //cout << "start to load preIndex ..." << endl;
    log_batch_ofs << "start to load preIndex ..." << endl;
    string indexStr = optionInfo->global_index_file_path_prefix; //argv[6];
    string preIndexArrayPreStr = indexStr;
    string secondLevelIndexStr = optionInfo->local_index_file_path_prefix; //argv[7];
    preIndexArrayPreStr.append("/");
    indexStr.append("/");
    secondLevelIndexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); 
	ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); 
	ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); 
	ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); 
	preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); 
	preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	log_batch_ofs << "finish loading preIndex ..." << endl;
	
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); 
	ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_batch_ofs);
	settings_batch_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_batch_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	settings_batch_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	settings_batch_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	log_batch_ofs << "finish loading chromosomes" << endl;
 	
	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; 
    unsigned int *childTab;
    BYTE *lcpCompress;
    BYTE *verifyChild;
    sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE));
	verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE));
	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";

	log_batch_ofs << "start to load SA, lcpCompress, childTab, detChild ....." << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_batch_ofs << "All index files loaded" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... all index loaded ......" << endl; 
	log_batch_ofs << endl << tmpTimeStr << "... all index loaded ......" << endl;
	runtime_batch_ofs << endl << tmpTimeStr << "... all index loaded ......" << endl;	

	//////////////////////////////////////////////////
	HeaderSection_Info* headerSectionInfo = new HeaderSection_Info(indexInfo);	
	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	//cout << endl << tmpTimeStr << "... start to load annotation file (SJs)......" << endl; 
	log_batch_ofs << endl << tmpTimeStr << "... start to load annotation file (SJs) ......" << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path; // junction files
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////	
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_shared = new AlignInferJunctionHash_Info();
	bool SJalignInferHash_provided_bool = optionInfo->spliceJunctionAlignInferHash_provided_bool;
	//cout << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	log_batch_ofs << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo_shared->initiateAlignInferJunctionInfo(chromNum);
	if(SJalignInferHash_provided_bool)
	{	
		//cout << "start to insert SJ into SJmap" << endl;
		string tmpInputJuncFile = optionInfo->spliceJunctionAlignInferHash_file_path;
		alignInferJunctionHashInfo_shared->insertJuncFromJuncFile_onlyChrNamePos(tmpInputJuncFile, indexInfo);
		//cout << "start to output SJ map" << endl;
	}
	//////////////////////////////////////////////////       finish LOADing INDEX      ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	bool segMap2SNPmer_phase_defined_bool = optionInfo->segMap2SNPmer_phase_defined_bool;
	bool segMap2SNPmer_bothPhase_or_phase2Only_bool = optionInfo->segMap2SNPmer_bothPhase_or_phase2Only_bool;
	bool segMap2SNPmer_phase1_bool, segMap2SNPmer_phase2_bool;
	if(segMap2SNPmer_phase_defined_bool)
	{
		if(segMap2SNPmer_bothPhase_or_phase2Only_bool)
		{
			segMap2SNPmer_phase1_bool = true;
			segMap2SNPmer_phase2_bool = true;			
		}
		else
		{	
			segMap2SNPmer_phase1_bool = false;
			segMap2SNPmer_phase2_bool = true;
		}
		segMap2SNPmer_phase2_bool_learned = false;
		snpLearned_success_bool = false;
	}
	else
	{
		segMap2SNPmer_phase1_bool = false;
		segMap2SNPmer_phase2_bool = false;
	}
	cout << "segMap2SNPmer_phase_defined_bool: " << segMap2SNPmer_phase_defined_bool << endl;
	cout << "segMap2SNPmer_phase1_bool: " << segMap2SNPmer_phase1_bool << endl;
	cout << "segMap2SNPmer_phase2_bool: " << segMap2SNPmer_phase2_bool << endl;
	cout << "segMap2SNPmer_phase2_bool_learned: " << segMap2SNPmer_phase2_bool_learned << endl;
	cout << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
	log_batch_ofs << "segMap2SNPmer_phase_defined_bool: " << segMap2SNPmer_phase_defined_bool << endl;
	log_batch_ofs << "segMap2SNPmer_phase1_bool: " << segMap2SNPmer_phase1_bool << endl;
	log_batch_ofs << "segMap2SNPmer_phase2_bool: " << segMap2SNPmer_phase2_bool << endl;
	log_batch_ofs << "segMap2SNPmer_phase2_bool_learned: " << segMap2SNPmer_phase2_bool_learned << endl;
	log_batch_ofs << "snpLearned_success_bool: " << snpLearned_success_bool << endl; 

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////   initiate batch manager info  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	log_batch_ofs << "start to initiate batchManagerInfo" << endl;
	//cout << "start to initiate batchManagerInfo" << endl;
	Batch_Manager_Info tmpBatchManagerInfo;
	bool readFileIn_listFile_or_commandLine_bool = optionInfo->return_readFileIn_listFile_or_commandLine_bool();
	if(readFileIn_listFile_or_commandLine_bool)
		tmpBatchManagerInfo.initiate_readFileVec_withReadFileListFile(InputReadFile, InputReadFile_PE);
	else
		tmpBatchManagerInfo.initaite_readFileVec_withReadFileInCommandLine(InputReadFile, InputReadFile_PE);
	tmpBatchManagerInfo.initiate_dataSetNum_accordingToReadFileVec();
    string outputDirStr = optionInfo->outputFolder_path; //argv[3];
    tmpBatchManagerInfo.initiate_create_resultFolderVec_accordingToReadFileVec(outputDirStr);
    tmpBatchManagerInfo.initiate_SNPlearnedSuccessBoolVec_withDataSetNum();
    tmpBatchManagerInfo.initiateStatsInfoVec_withDataSetNum(SE_or_PE_bool, threads_num);
    tmpBatchManagerInfo.initaite_totalReadNumVec_withDataSetNum();
    int batch_dataSet_num = tmpBatchManagerInfo.return_dataSet_num();
  	log_batch_ofs << "end of initiating batchManagerInfo " << endl;
  	log_batch_ofs << "batch_dataSet_num: " << batch_dataSet_num << endl;
  	//cout << "end of initiating batchManagerInfo " << endl;
  	cout << "batch_dataSet_num: " << batch_dataSet_num << endl;  	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start of phase1 for all dataSet ..." << endl;
	runtime_batch_ofs << tmpTimeStr << "... start of phase1 for all dataSet ..." << endl;
    for(int tmpData = 0; tmpData < batch_dataSet_num; tmpData ++)
    {
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << tmpTimeStr << "... start of phase1 for dataSet " << tmpData + 1 << endl;
		runtime_batch_ofs << tmpTimeStr << "... start of phase1 for dataSet " << tmpData + 1 << endl;
		//tmpData_log_phase1_ofs << endl << tmpTimeStr << "... 1st mapping process starts " << tmpData + 1 << endl;

    	string tmpData_outputDirStr = tmpBatchManagerInfo.returnResultFolder_withIndexInVec(tmpData);
    	string tmpData_outputDirStr_logs = tmpData_outputDirStr + "/logs";
	   	string tmpData_mkdirOutputCommand_log = "mkdir -p " + tmpData_outputDirStr_logs;
   		system(tmpData_mkdirOutputCommand_log.c_str());

    	string tmpData_log_phase1_fie = tmpData_outputDirStr + "/logs/log_phase1.txt";
    	string tmpData_settings_phase1_file = tmpData_outputDirStr + "/logs/settings_phase1.txt";
    	string tmpData_stats_phase1_file = tmpData_outputDirStr + "/logs/stats_phase1.txt";
    	ofstream tmpData_log_phase1_ofs(tmpData_log_phase1_fie.c_str());
    	ofstream tmpData_settings_phase1_ofs(tmpData_settings_phase1_file.c_str());
    	ofstream tmpData_stats_phase1_ofs(tmpData_stats_phase1_file.c_str());
    	int tmpData_readTotalNum = 0;
		//////////////////////////////////////////////////   #BEGIN#    Personal Genome           ////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////   #BEGIN#    Personal Genome           ////////////////////////////////////////////////////////////////
		Index_Info* tmpData_indexInfo_SNP = new Index_Info();
		char *tmpData_chrom_SNP;
	    unsigned int *tmpData_sa_SNP; 
	    BYTE *tmpData_lcpCompress_SNP;
	    unsigned int *tmpData_childTab_SNP;
	    BYTE *tmpData_verifyChild_SNP;
	    string tmpData_SNPfilePath;
		if(segMap2SNPmer_phase1_bool)
		{	
			tmpData_SNPfilePath = tmpBatchManagerInfo.returnSNPfile_withIndexInVec(tmpData);
			string tmpData_SNP_seq_index_path = tmpBatchManagerInfo.returnSNPmerIndexFolderPath_withIndexInVec(tmpData);
			//cout << "SNPfilePath: " << endl << tmpData_SNPfilePath << endl;
			tmpData_settings_phase1_ofs << "SNPfilePath: " << endl << tmpData_SNPfilePath << endl;
			// Xinan: for now, do not replace ref bases with alternate bases. In the future, 
			// MPS3 should incorporate SNPs when doing sequence matching and canonical splice site detection
			indexInfo->insertSNP2chromStr(tmpData_SNPfilePath, tmpData_log_phase1_ofs);
			//cout << "start to load indexes" << endl;
			string tmpData_indexStr_SNP = tmpData_SNP_seq_index_path;
			//cout << "tmpData_indexStr_SNP: " << endl << tmpData_indexStr_SNP << endl;
			tmpData_settings_phase1_ofs << "tmpData_indexStr_SNP: " << endl << tmpData_indexStr_SNP << endl;
			tmpData_indexStr_SNP.append("/");
			string tmpData_SA_file_SNP = tmpData_indexStr_SNP + "_SA"; 
			ifstream tmpData_SA_file_ifs_SNP(tmpData_SA_file_SNP.c_str(),ios::binary); 
			string tmpData_lcpCompress_file_SNP = tmpData_indexStr_SNP + "_lcpCompress"; 
			ifstream tmpData_lcpCompress_file_ifs_SNP(tmpData_lcpCompress_file_SNP.c_str(),ios::binary);
			string tmpData_childTab_file_SNP = tmpData_indexStr_SNP + "_childTab"; 
			ifstream tmpData_childTab_file_ifs_SNP(tmpData_childTab_file_SNP.c_str(),ios::binary);
			string tmpData_verifyChild_file_SNP = tmpData_indexStr_SNP + "_detChild"; 
			ifstream tmpData_verifyChild_file_ifs_SNP(tmpData_verifyChild_file_SNP.c_str(),ios::binary);	
			string tmpData_chrom_bit_file_SNP = tmpData_indexStr_SNP + "_chrom"; 
			ifstream tmpData_chrom_bit_file_ifs_SNP(tmpData_chrom_bit_file_SNP.c_str(),ios::binary);
			string tmpData_parameter_file_SNP = tmpData_indexStr_SNP + "_parameter"; 
			ifstream tmpData_parameter_file_ifs_SNP(tmpData_parameter_file_SNP.c_str(),ios::binary);

			tmpData_indexInfo_SNP->initiate(tmpData_parameter_file_ifs_SNP, tmpData_log_phase1_ofs);
			tmpData_chrom_SNP = (char*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
			tmpData_chrom_bit_file_ifs_SNP.read((char*)tmpData_chrom_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
			tmpData_indexInfo_SNP->readGenome(tmpData_chrom_SNP);
			tmpData_indexInfo_SNP->initiate();	
			tmpData_indexInfo_SNP->initiateChrNameIndexArray(1000);
		    tmpData_sa_SNP = (unsigned int*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); 
		    tmpData_SA_file_ifs_SNP.read((char*)tmpData_sa_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_lcpCompress_SNP = (BYTE*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_lcpCompress_file_ifs_SNP.read((char*)tmpData_lcpCompress_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));	 
			tmpData_childTab_SNP = (unsigned int*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); 
			tmpData_childTab_file_ifs_SNP.read((char*)tmpData_childTab_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_verifyChild_SNP = (BYTE*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_verifyChild_file_ifs_SNP.read((char*)tmpData_verifyChild_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
			tmpData_log_phase1_ofs << "SyntheticSNPtransSeq index files loaded" << endl;
		
			tmpData_SA_file_ifs_SNP.close();
			tmpData_lcpCompress_file_ifs_SNP.close();
			tmpData_childTab_file_ifs_SNP.close();
			tmpData_verifyChild_file_ifs_SNP.close();
			tmpData_chrom_bit_file_ifs_SNP.close();
			tmpData_parameter_file_ifs_SNP.close();
		}
		//////////////////////////////////////////////////   #END#    Personal Genome           ////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////   #END#    Personal Genome           ////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
		///////////////////////   initiate log files    ///////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
	    
	    string tmpData_outputDirStr_logs_phase1 = tmpData_outputDirStr_logs + "/phase1_log";
	    string tmpData_outputDirStr_logs_phase2 = tmpData_outputDirStr_logs + "/phase2_log";
	    string tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped = tmpData_outputDirStr_logs_phase2 + "/fixOneEndUnmapped";
	    string tmpData_outputDirStr_logs_phase2_fixHeadTail = tmpData_outputDirStr_logs_phase2 + "/fixHeadTail";
   		string tmpData_mkdirOutputCommand_log_phase1 = "mkdir -p " + tmpData_outputDirStr_logs_phase1;
  	 	system(tmpData_mkdirOutputCommand_log_phase1.c_str());   	
   		string tmpData_mkdirOutputCommand_log_phase2 = "mkdir -p " + tmpData_outputDirStr_logs_phase2;
   		system(tmpData_mkdirOutputCommand_log_phase2.c_str());
   		string tmpData_mkdirOutputCommand_log_phase2_fixOneEndUnmapped = "mkdir -p " + tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped;
   		system(tmpData_mkdirOutputCommand_log_phase2_fixOneEndUnmapped.c_str());
   		string tmpData_mkdirOutputCommand_log_phase2_fixHeadTail = "mkdir -p " + tmpData_outputDirStr_logs_phase2_fixHeadTail;
   		system(tmpData_mkdirOutputCommand_log_phase2_fixHeadTail.c_str());

	   	string tmpData_inputLogStr_phase1 = tmpData_outputDirStr_logs_phase1 + "/input.log";
	   	ofstream tmpData_input_log_ofs_phase1(tmpData_inputLogStr_phase1.c_str());
	   	string tmpData_outputLogStr_phase1 = tmpData_outputDirStr_logs_phase1 + "/output.log";
	   	ofstream tmpData_output_log_ofs_phase1(tmpData_outputLogStr_phase1.c_str());
	   	string tmpData_mappingLogStr_phase1 = tmpData_outputDirStr_logs_phase1 + "/mapping.log";
	   	ofstream tmpData_mapping_log_ofs_phase1(tmpData_mappingLogStr_phase1.c_str());   	
	   	optionInfo->outputOptStr(tmpData_settings_phase1_ofs);
		//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
		/* align main*/
	   	string tmpData_tmpHeadSectionInfo = tmpData_outputDirStr + "/headSectionInfo";
	   	ofstream tmpData_tmpHeadSectionInfo_ofs(tmpData_tmpHeadSectionInfo.c_str());    	
	   	tmpData_tmpHeadSectionInfo_ofs << headerSectionInfo->returnHeaderSectionInfoStr() << endl;

		string tmpData_mkdirOutputCommand_phase1 = "mkdir -p " + tmpData_outputDirStr + "/phase1_output";
		system(tmpData_mkdirOutputCommand_phase1.c_str());
		string tmpData_mkdirOutputCommand_repeatRegionFile = tmpData_mkdirOutputCommand_phase1 + "/repeat_region";
	   	system(tmpData_mkdirOutputCommand_repeatRegionFile.c_str());
		string tmpData_repeatRegionFile = tmpData_outputDirStr + "/phase1_output/repeat_region/repeatRegion";
		ofstream tmpData_repeatRegionFile_ofs(tmpData_repeatRegionFile.c_str());
		string tmpData_mkdirOutputCommand_tmpAlignCompleteRead = tmpData_mkdirOutputCommand_phase1 + "/completePair";
		system(tmpData_mkdirOutputCommand_tmpAlignCompleteRead.c_str());
		string tmpData_tmpAlignCompleteRead_SE = tmpData_outputDirStr + "/phase1_output/completePair/complete_SE.sam";
		ofstream tmpData_tmpAlignCompleteRead_SE_ofs(tmpData_tmpAlignCompleteRead_SE.c_str());
		string tmpData_tmpAlignCompleteRead = tmpData_outputDirStr + "/phase1_output/completePair/completePair.sam";
		ofstream tmpData_tmpAlignCompleteRead_ofs(tmpData_tmpAlignCompleteRead.c_str());
		string tmpData_mkdirOutputCommand_tmpAlignOneEndUnmapped = tmpData_mkdirOutputCommand_phase1 + "/oneEndUnmapped";
		system(tmpData_mkdirOutputCommand_tmpAlignOneEndUnmapped.c_str());
		string tmpData_tmpAlignOneEndUnmapped = tmpData_outputDirStr + "/phase1_output/oneEndUnmapped/oneEndUnmapped";
		if(Do_Phase1_Only)
			tmpData_tmpAlignOneEndUnmapped += ".sam";	
		else
			tmpData_tmpAlignOneEndUnmapped += ".alignInfo";
		ofstream tmpData_tmpAlignOneEndUnmapped_ofs(tmpData_tmpAlignOneEndUnmapped.c_str());

		string tmpData_tmpAlignUnmapped_mappedToRepeatRegionFile_SE = tmpData_outputDirStr + "/phase1_output/repeat_region/unmapped_mappedToRepeatRegion_SE.sam";
		ofstream tmpData_tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs(tmpData_tmpAlignUnmapped_mappedToRepeatRegionFile_SE.c_str());
		string tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = tmpData_outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam";
		ofstream tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs(tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
		string tmpData_mkdirOutputCommand_tmpAlignBothEndsUnmapped = tmpData_mkdirOutputCommand_phase1 + "/bothEndsUnmapped";
		system(tmpData_mkdirOutputCommand_tmpAlignBothEndsUnmapped.c_str());
		string tmpData_tmpAlignUnmapped_SE = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_SE.sam";
		ofstream tmpData_tmpAlignUnmapped_SE_ofs(tmpData_tmpAlignUnmapped_SE.c_str());
		string tmpData_tmpAlignBothEndsUnmapped = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
		ofstream tmpData_tmpAlignBothEndsUnmapped_ofs(tmpData_tmpAlignBothEndsUnmapped.c_str());
		string tmpData_tmpAlignUnmapped_lowScore_SE = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_lowScore_SE.sam";
		ofstream tmpData_tmpAlignUnmapped_lowScore_SE_ofs(tmpData_tmpAlignUnmapped_lowScore_SE.c_str());
		string tmpData_tmpAlignBothEndsUnmapped_lowScore = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
		ofstream tmpData_tmpAlignBothEndsUnmapped_lowScore_ofs(tmpData_tmpAlignBothEndsUnmapped_lowScore.c_str());
		string tmpData_mkdirOutputCommand_tmpAlignIncompletePair = tmpData_mkdirOutputCommand_phase1 + "/incomplete";
		system(tmpData_mkdirOutputCommand_tmpAlignIncompletePair.c_str());
		string tmpData_tmpAlignIncomplete_SE = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete_SE.alignInfo"; 
		ofstream tmpData_tmpAlignIncomplete_SE_ofs(tmpData_tmpAlignIncomplete_SE);
		string tmpData_tmpAlignIncompletePair = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo"; 
		ofstream tmpData_tmpAlignIncompletePair_ofs(tmpData_tmpAlignIncompletePair.c_str());
		string tmpData_tmpAlignIncomplete_SE_SAM = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete_SE.sam";
		ofstream tmpData_tmpAlignIncomplete_SE_SAM_ofs(tmpData_tmpAlignIncomplete_SE_SAM.c_str());
		string tmpData_tmpAlignIncompletePair_SAM = tmpData_outputDirStr + "/phase1_output/incomplete/incompletePair.sam"; 
		ofstream tmpData_tmpAlignIncompletePair_SAM_ofs(tmpData_tmpAlignIncompletePair_SAM.c_str());	
		string tmpData_tmpIntermediateJunctionFile = tmpData_outputDirStr + "/phase2_output/inter.junc";

	    string tmpData_read_file_1 = tmpBatchManagerInfo.returnReadFile_1_withIndexInVec(tmpData);
	    string tmpData_read_file_2 = tmpBatchManagerInfo.returnReadFile_2_withIndexInVec(tmpData);
		ifstream tmpData_inputRead_ifs(tmpData_read_file_1.c_str());
		ifstream tmpData_inputRead_PE_ifs(tmpData_read_file_2.c_str());
    
	    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE, line2_afterProcess, line2_PE_afterProcess;
		vector< RepeatRegion_Info* > tmpData_repeatRegionInfoVec;
		for(int tmp = 0; tmp < threads_num; tmp++)
		{
			RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
			tmpData_repeatRegionInfoVec.push_back(repeatRegionInfo);
		}

		////////////////////    initiating   learning SNPs from phase1 (by default)        ////////////////////////////
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "start to initiate learnedCandiSNPhashInfo for dataSet " << tmpData + 1 << endl; 
		tmpData_log_phase1_ofs << endl << tmpTimeStr << "start to initiate learnedCandiSNPhashInfo for dataSet " << tmpData + 1 << endl;
		LearnedCandiSNPhash_Info_Vec tmpData_mismatchHashInfoVec;
		int mismatchHashInfoVecSize = threads_num;
		tmpData_mismatchHashInfoVec.initiateLearnedCandiSNPhashInfoVec(mismatchHashInfoVecSize, chromNum);

		InputReadPreProcess* tmpData_readPreProcessInfo = new InputReadPreProcess();
		Read_Array_Queue* tmpData_readArrayQueue = new Read_Array_Queue();
		Result_Array_Queue* tmpData_resultArrayQueue = new Result_Array_Queue();
		int tmpData_tmpInputReadNumInBatchArray_phase1 = normalRecordNum_1stMapping;//ReadNumInReadArray_Phase1;
		int tmpData_tmpInputTimeWeight_phase1 = InputTimeWeight_Phase1;
		int tmpData_tmpOutputTimeWeigth_phase1 = OutputTimeWeight_Phase1; 
		bool tmpData_endOfFile_bool = false;
		bool tmpData_endOfProcessing_bool = false;

		if(separateThreadForIO_bool)
		{	
			omp_set_num_threads(2);
			omp_set_nested(1);
		#pragma omp parallel
			{
		#pragma omp sections
				{
		#pragma omp section
					io_stage_phase1_separateThreadForIO(tmpData_inputRead_ifs, tmpData_inputRead_PE_ifs, tmpData_readArrayQueue, tmpData_resultArrayQueue, 
						tmpData_endOfFile_bool, tmpData_endOfProcessing_bool, 
						tmpData_tmpInputReadNumInBatchArray_phase1, tmpData_tmpInputTimeWeight_phase1, tmpData_tmpOutputTimeWeigth_phase1, tmpData_log_phase1_ofs, 
						tmpData_readPreProcessInfo, tmpData_tmpAlignCompleteRead_ofs, tmpData_tmpAlignIncompletePair_ofs, tmpData_tmpAlignOneEndUnmapped_ofs, 
						tmpData_tmpAlignBothEndsUnmapped_ofs, tmpData_tmpAlignBothEndsUnmapped_lowScore_ofs, tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
						tmpData_tmpAlignIncompletePair_SAM_ofs, tmpData_input_log_ofs_phase1, tmpData_output_log_ofs_phase1, fasta_or_fastq_bool, SE_or_PE_bool,
						tmpData_readTotalNum);
		#pragma omp section
					process_stage_phase1_separateThreadForIO(tmpData_readArrayQueue, tmpData_resultArrayQueue, tmpData_endOfFile_bool, tmpData_endOfProcessing_bool, 
						threads_num-1, sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, preIndexIntervalStartArray,
						preIndexIntervalEndArray, tmpData_repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
						Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only, Do_Phase1_Only,	tmpBatchManagerInfo.statsInfoVec[tmpData], 
						fasta_or_fastq_bool, 
						tmpData_mapping_log_ofs_phase1, checkQualSeqForReadSegSeq, SE_or_PE_bool,
						tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, 
						tmpData_indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase1_bool, tmpData_mismatchHashInfoVec, segMap2SNPmer_phase2_bool_learned);
				}
			}
		}
		else
		{
			io_process_phase1_allThreadsSharedByBothStage(tmpData_inputRead_ifs, tmpData_inputRead_PE_ifs, tmpData_readArrayQueue, tmpData_resultArrayQueue,
				tmpData_tmpInputReadNumInBatchArray_phase1, tmpData_log_phase1_ofs, tmpData_readPreProcessInfo, tmpData_tmpAlignCompleteRead_ofs, 
				tmpData_tmpAlignIncompletePair_ofs, tmpData_tmpAlignOneEndUnmapped_ofs, tmpData_tmpAlignBothEndsUnmapped_ofs, tmpData_tmpAlignBothEndsUnmapped_lowScore_ofs, 
				tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs, tmpData_tmpAlignIncompletePair_SAM_ofs, tmpData_input_log_ofs_phase1, tmpData_output_log_ofs_phase1, 
				fasta_or_fastq_bool, SE_or_PE_bool, 
				threads_num, sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, 
				preIndexIntervalStartArray, preIndexIntervalEndArray, tmpData_repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only, Do_Phase1_Only,	tmpBatchManagerInfo.statsInfoVec[tmpData], 
				tmpData_mapping_log_ofs_phase1, 
				checkQualSeqForReadSegSeq, tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, tmpData_indexInfo_SNP, 
				SNPlocInSyntheticSNPseq, segMap2SNPmer_phase1_bool, tmpData_mismatchHashInfoVec, segMap2SNPmer_phase2_bool_learned, tmpData_readTotalNum);
		}
		//log_batch_ofs << "tmpData_readTotalNum: " << tmpData_readTotalNum << " for dataSet " << tmpData + 1 << endl;
		//cout << "tmpData_readTotalNum: " << tmpData_readTotalNum << " for dataSet " << tmpData + 1 << endl;
		tmpBatchManagerInfo.assignTotalReadNum(tmpData_readTotalNum, tmpData);
		//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
		tmpData_repeatRegionFile_ofs << "Repeat Region Info: size = " << tmpData_repeatRegionInfoVec.size() << endl;
		for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
			tmpData_repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, tmpData_repeatRegionFile_ofs);
		for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
		{
			delete tmpData_repeatRegionInfoVec[tmpThread];
			tmpData_repeatRegionInfoVec[tmpThread] = NULL;
		}
		/////////////////////  convert SNPbases back to reference bases    ///////////////////////////////////////
		if(segMap2SNPmer_phase1_bool)
			indexInfo->convertSNPbaseBackToReferenceBase(tmpData_SNPfilePath, tmpData_log_phase1_ofs);
		///////////////////////////////////////////      merge learned candi SNPs      ////////////////////////////////////////////
		if(segMap2SNPmer_phase2_bool_learned)
		{
			LearnedCandiSNPhash_Info tmpData_mismatchHashInfo_merged;
			tmpData_mismatchHashInfo_merged.initiate(chromNum);
			// merge to one
			tmpData_mismatchHashInfoVec.mergeLearnedCandiSNPhashInfoVec2one(tmpData_mismatchHashInfo_merged, indexInfo);
			tmpData_mismatchHashInfoVec.freeMemory();
			string tmpData_dir_SNPmer_learned = tmpData_outputDirStr + "/SNPmer_learned";
		   	string tmpData_mkdir_SNPmer_learned = "mkdir -p " + tmpData_dir_SNPmer_learned;
		   	system(tmpData_mkdir_SNPmer_learned.c_str());
		   	string tmpData_filteredLearnedSNP_file = tmpData_dir_SNPmer_learned + "/filteredLearnedSNP.txt";
		   	int tmpData_called_snp_num = 0;
			int tmpData_filterLearnedSNP_supNumMin = FILTER_LEARNSNP_SUPNUMMIN;
			double tmpData_filterLearnedSNP_ratioMin = FILTER_LEARNSNP_RATIOMIN;
			tmpData_mismatchHashInfo_merged.output_filteredLearnedSNP(tmpData_filteredLearnedSNP_file, tmpData_called_snp_num, 
				tmpData_filterLearnedSNP_supNumMin, tmpData_filterLearnedSNP_ratioMin, indexInfo);
			tmpData_settings_phase1_ofs << "filterLearnedSNP_supNumMin: " << tmpData_filterLearnedSNP_supNumMin << " for dataSet " << tmpData + 1 << endl;
			tmpData_settings_phase1_ofs << "filterLearnedSNP_ratioMin: " << tmpData_filterLearnedSNP_ratioMin << " for dataSet " << tmpData + 1 << endl;
			tmpData_settings_phase1_ofs << "called_snp_num: " << tmpData_called_snp_num << " for dataSet " << tmpData + 1 << endl;
			if(tmpData_called_snp_num == 0)
				tmpBatchManagerInfo.assign_SNPlearnedSuccessBool(false, tmpData);
	    	else
	    		tmpBatchManagerInfo.assign_SNPlearnedSuccessBool(true, tmpData);
    	}
    	else
    		tmpBatchManagerInfo.assign_SNPlearnedSuccessBool(false, tmpData);

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... end of phase1 for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase1_ofs << endl << tmpTimeStr << "... end of phase1" << endl;
		log_batch_ofs << endl << tmpTimeStr << "... end of phase1 for dataSet " << tmpData + 1 << endl;
		runtime_batch_ofs << endl << tmpTimeStr << "... end of phase1 for dataSet " << tmpData + 1 << endl;
		//tmpData_log_phase1_ofs << endl << "**********************************" << endl << "**********************************";
		
		tmpData_inputRead_ifs.close();
		tmpData_inputRead_PE_ifs.close();
		tmpData_log_phase1_ofs.close();
		tmpData_settings_phase1_ofs.close();

	   	tmpData_input_log_ofs_phase1.close();
	   	tmpData_output_log_ofs_phase1.close();
	   	tmpData_mapping_log_ofs_phase1.close();	
 
	   	tmpData_tmpHeadSectionInfo_ofs.close(); 	
		tmpData_repeatRegionFile_ofs.close();
		tmpData_tmpAlignCompleteRead_SE_ofs.close();
		tmpData_tmpAlignCompleteRead_ofs.close();

		tmpData_tmpAlignOneEndUnmapped_ofs.close();
		tmpData_tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs.close();
		tmpData_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs.close();
		tmpData_tmpAlignUnmapped_SE_ofs.close();
		tmpData_tmpAlignBothEndsUnmapped_ofs.close();
		tmpData_tmpAlignUnmapped_lowScore_SE_ofs.close();
		tmpData_tmpAlignBothEndsUnmapped_lowScore_ofs.close();
		tmpData_tmpAlignIncomplete_SE_ofs.close();
		tmpData_tmpAlignIncompletePair_ofs.close();
		tmpData_tmpAlignIncomplete_SE_SAM_ofs.close();
		tmpData_tmpAlignIncompletePair_SAM_ofs.close();

		//cout << "start to initiate and free arrays" << endl;

		delete tmpData_readPreProcessInfo;
		tmpData_readPreProcessInfo = NULL;
		delete tmpData_readArrayQueue;
		tmpData_readArrayQueue = NULL;
		delete tmpData_resultArrayQueue;
		tmpData_resultArrayQueue = NULL;
		if(segMap2SNPmer_phase1_bool)
		{
			free(tmpData_sa_SNP);
			free(tmpData_chrom_SNP);
	    	free(tmpData_lcpCompress_SNP);
	    	free(tmpData_childTab_SNP);
	    	free(tmpData_verifyChild_SNP);
    	}
    	tmpData_sa_SNP = NULL;
    	tmpData_chrom_SNP = NULL;
    	tmpData_lcpCompress_SNP = NULL;
    	tmpData_childTab_SNP = NULL;
    	tmpData_verifyChild_SNP = NULL;
    	delete tmpData_indexInfo_SNP;
    	tmpData_indexInfo_SNP = NULL;
    }

	free(preIndexMapLengthArray); preIndexMapLengthArray = NULL;
	free(preIndexIntervalStartArray); preIndexIntervalStartArray = NULL;
	free(preIndexIntervalEndArray); preIndexIntervalEndArray = NULL;
	free(sa); sa = NULL;
	free(lcpCompress); lcpCompress = NULL;
	free(childTab); childTab = NULL;
	free(verifyChild); verifyChild = NULL;
	free(chrom); chrom = NULL;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";	
	cout << endl << tmpTimeStr << "... end of phase 1 for all datasets ..." << endl; 
	log_batch_ofs << endl << tmpTimeStr << "... end of phase 1 for all datasets ..." << endl; 
	runtime_batch_ofs << endl << tmpTimeStr << "... end of phase 1 for all datasets ..." << endl;

	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";	
	cout << endl << endl << tmpTimeStr << "... start of loading local index ..." << endl;
	log_batch_ofs << endl << endl << tmpTimeStr << "... start of loading local index ..." << endl;
	runtime_batch_ofs << endl << endl << tmpTimeStr << "... start of loading local index ..." << endl; 
	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;
	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;
	if(load2ndLevelIndexBool)
	{
		//log_batch_ofs << "start to load second-level index ..." << endl;		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->returnChromNum(); tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;				
				string inputIndexFileStr = secondLevelIndexStr + "/" + indexInfo->returnChrNameStr(tmpChrNO) + "/" + tmpFileNumStr + "/";
				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);
				string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
				ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
				string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
				ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
				string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
				ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					

				int sizeOfIndex = indexInfo->returnSecondLevelIndexNormalSize() + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);

				bool No_ATGC_Bool = true;
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					char ch = tmpSecondLevelChrom[tmpMallocSpace];
					if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
					{
						No_ATGC_Bool = false;
						break;
					}
				}				
				if(No_ATGC_Bool)
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				secondLevelChrom.push_back(tmpSecondLevelChrom);				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);
				BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);					
				unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab.push_back(tmpSecondLevelChildTab);
				BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild.push_back(tmpSecondLevelDetChild);

				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();
				secondLevelLcpCompress_file_ifs.close();
				secondLevelChildTab_file_ifs.close();
				secondLevelDetChild_file_ifs.close();							
				secondLevelIndexNO ++;
			}
			//log_batch_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
		}
		//log_batch_ofs << "finish loading ALL 2nd-level index !" << endl;
		//log_batch_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
		//loadIndex_end = clock(); 
	}
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... end of loading local index ..." << endl << endl;
	log_batch_ofs << endl << tmpTimeStr << "... end of loading local index ..." << endl << endl;
	runtime_batch_ofs << endl << tmpTimeStr << "... end of loading local index ..." << endl << endl;
	//cout << endl << "**********************************" << endl << "**********************************" << endl;
	//runtime_batch_ofs << endl << "**********************************" << endl << "**********************************" << endl;
   	////////////////////////    	start to do sam2junc for sharing among all the data sets      ///////////////////////////
	SJhash_Info* shared_SJ = new SJhash_Info();
	bool spliceJunctionHashExists_shared = false;
	if(sharingSplicingContext_bool)
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";	
		cout << endl << endl << tmpTimeStr << "... start to do sam2junc for sharing among all the data sets ..." << endl; 
		log_batch_ofs << endl << endl << tmpTimeStr << "... start to do sam2junc for sharing among all the data sets ..." << endl;
		runtime_batch_ofs << endl << endl << tmpTimeStr << "... start to do sam2junc for sharing among all the data sets ..." << endl;
		string shared_junc_file = outputDirStr + "/shared_junc.txt";
		vector<string> toShareJuncSamPathVec;
		for(int tmpData = 0; tmpData < batch_dataSet_num; tmpData ++)
		{
			string tmpData_outputDirStr = tmpBatchManagerInfo.returnResultFolder_withIndexInVec(tmpData);
			if(SE_or_PE_bool)
			{
				string tmpData_tmpAlignCompleteRead_SE = tmpData_outputDirStr + "/phase1_output/completePair/complete_SE.sam";
				string tmpData_tmpAlignIncomplete_SE_SAM = tmpData_outputDirStr +  "/phase1_output/incomplete/incomplete_SE.sam";
				toShareJuncSamPathVec.push_back(tmpData_tmpAlignCompleteRead_SE);
				toShareJuncSamPathVec.push_back(tmpData_tmpAlignIncomplete_SE_SAM);
			}
			else
			{
				string tmpData_tmpAlignCompleteRead = tmpData_outputDirStr + "/phase1_output/completePair/completePair.sam";
				string tmpData_tmpAlignIncompletePair_SAM = tmpData_outputDirStr + "/phase1_output/incomplete/incompletePair.sam";
				//string tmpData_OutputSamFile_oneEndMapped = "/phase2_output/oneEndUnmapped.pairedComplete.sam";
				toShareJuncSamPathVec.push_back(tmpData_tmpAlignCompleteRead);
				toShareJuncSamPathVec.push_back(tmpData_tmpAlignIncompletePair_SAM);
				//toShareJuncSamPathVec.push_back(tmpData_OutputSamFile_oneEndMapped);
			}
		}
		//log_batch_ofs << "start to load SJs in alignments" << endl;
		//cout << "start to load SJs in alignments" << endl;
		AlignInferJunctionHash_Info_Vec* shared_alignInferJunctionHashInfoVec = new AlignInferJunctionHash_Info_Vec();
		int alignInferJuncHashInfoVecSize = threads_num;
		shared_alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(alignInferJuncHashInfoVecSize, chromNum);		
		shared_alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(
			toShareJuncSamPathVec, indexInfo, alignInferJuncHashInfoVecSize, log_batch_ofs);
		shared_alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(alignInferJunctionHashInfo_shared, indexInfo);
		alignInferJunctionHashInfo_shared->outputAlignInferInfoHashInfo_chrNamePos_supportNum(indexInfo, shared_junc_file);
		shared_alignInferJunctionHashInfoVec->freeMemory();
		delete shared_alignInferJunctionHashInfoVec;
		shared_alignInferJunctionHashInfoVec = NULL;
	  	//log_batch_ofs << "start to convert shared sam2alignInferJuncHash 2 SJhashInfo" << endl;
		//cout << "start to convert shared sam2alignInferJuncHash 2 SJhashInfo" << endl;
			
		shared_SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());			
		alignInferJunctionHashInfo_shared->convert2SJhashInfo(shared_SJ, indexInfo);
		int junctionNum_in_alignInferJuncHashInfo_shared = alignInferJunctionHashInfo_shared->returnAlignInferInfoVecSize();		
		if(junctionNum_in_alignInferJuncHashInfo_shared == 0)
			spliceJunctionHashExists_shared = false;
		else
			spliceJunctionHashExists_shared = true;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... end of sam2junc for sharing among all the data sets ......" << endl << endl;
		log_batch_ofs << endl << tmpTimeStr << "... end of sam2junc for sharing among all the data sets ......" << endl << endl;
		runtime_batch_ofs << endl << tmpTimeStr << "... end of sam2junc for sharing among all the data sets ......" << endl << endl;
		//log_batch_ofs << endl << "**********************************" << endl << "**********************************" << endl;
	}
	///////////////////////////////////////////////   do Phase2 mapping   ///////////////////////////////////////////////
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << endl << tmpTimeStr << "... start of phase2 for all datasets ..." << endl;
	log_batch_ofs << endl << endl << tmpTimeStr << "... start of phase2 for all datasets ..." << endl;
	runtime_batch_ofs << endl << endl << tmpTimeStr << "... start of phase2 for all datasets ..." << endl;
	//log_batch_ofs << endl << "**********************************" << endl << "**********************************" << endl;
	for(int tmpData = 0; tmpData < batch_dataSet_num; tmpData ++)
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... start of phase2 mapping for dataSet " << tmpData + 1 << endl;
		log_batch_ofs << endl << tmpTimeStr << "... start of phase2 mapping for dataSet " << tmpData + 1 << endl;
		runtime_batch_ofs << endl << tmpTimeStr << "... start of phase2 mapping for dataSet " << tmpData + 1 << endl;

		/////////////////////////////   build SNPmer index for learned SNPs from phase1    //////////////////////////////
		//nowtime = time(NULL);
		//local = localtime(&nowtime);
		//tmpTimeStr = asctime(local);
		//tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... start to build SNPmer index for SNPs learned from phase1 for dataSet " << tmpData + 1 << endl;
		//log_batch_ofs << endl << tmpTimeStr << "... start to build SNPmer index for SNPs learned from phase1 for dataSet " << tmpData + 1 << endl;
		//log_batch_ofs << endl << "**********************************" << endl << "**********************************" << endl;

		string tmpData_outputDirStr = tmpBatchManagerInfo.returnResultFolder_withIndexInVec(tmpData);
    	string tmpData_log_phase2_file = tmpData_outputDirStr + "/log_phase2.txt";
    	string tmpData_settings_phase2_file = tmpData_outputDirStr + "/settings_phase2.txt";
    	string tmpData_stats_phase2_file = tmpData_outputDirStr + "/stats_phase2.txt";
    	ofstream tmpData_log_phase2_ofs(tmpData_log_phase2_file.c_str());
    	ofstream tmpData_settings_phase2_ofs(tmpData_settings_phase2_file.c_str());
    	ofstream tmpData_stats_phase2_ofs(tmpData_stats_phase2_file.c_str());    

    	bool tmpData_snpLearned_success_bool = tmpBatchManagerInfo.returnSNPlearnedSuccessBool(tmpData);
    	string tmpData_dir_SNPmer_learned, tmpData_dir_SNPmer_index_learned;
    	if(segMap2SNPmer_phase2_bool_learned && tmpData_snpLearned_success_bool)
    	{	
			tmpData_dir_SNPmer_learned = tmpData_outputDirStr + "/SNPmer_learned";
		   	// pull out SNPmers ... 
		   	//cout << "start to pull out SNPmers" << endl;
		   	string tmpData_filteredLearnedSNP_file = tmpData_dir_SNPmer_learned + "/filteredLearnedSNP.txt";
			string tmpData_filteredLearnedSNP_SNPmer_fa_file = tmpData_dir_SNPmer_learned + "/filteredLearnedSNP_SNPmer.fa";
			if(tmpData_snpLearned_success_bool)
			 	indexInfo->insertSNP2chromStr_outputSNPmer(tmpData_filteredLearnedSNP_file, tmpData_filteredLearnedSNP_SNPmer_fa_file, 
			 		tmpData_log_phase2_ofs, SNPmerLength);
			// indexing SNPmers
			//cout << "start to do indexing SNPmers " << endl;
		   	tmpData_dir_SNPmer_index_learned = tmpData_dir_SNPmer_learned + "/index";
		   	//cout << "tmpData_filteredLearnedSNP_SNPmer_fa_file: " << endl << tmpData_filteredLearnedSNP_SNPmer_fa_file << endl;
		   	//cout << "tmpData_dir_SNPmer_index_learned: " << endl << tmpData_dir_SNPmer_index_learned << endl;
		 	string cmd_buildSNPmerIndex_tmpData = "./buildSNPmerIndex_embeded_mps_batch "
		 		+ tmpData_filteredLearnedSNP_SNPmer_fa_file + " " + tmpData_dir_SNPmer_index_learned;
		 	system(cmd_buildSNPmerIndex_tmpData.c_str());
		 	string tmpData_dir_SNPmer_index_learned_success_file = tmpData_dir_SNPmer_index_learned + "/results.txt";
		 	string tmpData_buildSNPmerIndex_result_str;
		 	ifstream tmpData_dir_SNPmer_index_learned_success_ifs(tmpData_dir_SNPmer_index_learned_success_file.c_str());
		 	getline(tmpData_dir_SNPmer_index_learned_success_ifs, tmpData_buildSNPmerIndex_result_str);
		 	tmpData_dir_SNPmer_index_learned_success_ifs.close();
		 	if(tmpData_buildSNPmerIndex_result_str == "success")
		 		tmpData_snpLearned_success_bool = true;
		 	else
		 		tmpData_snpLearned_success_bool = false;

			// Build_snpMerIndex* tmpData_buildSnpMerIndexInfo = new Build_snpMerIndex();
			// if(tmpData_snpLearned_success_bool)
			// 	tmpData_snpLearned_success_bool = tmpData_buildSnpMerIndexInfo->build_snpMerIndex(
			// 		tmpData_filteredLearnedSNP_SNPmer_fa_file, tmpData_dir_SNPmer_index_learned); 
			// delete tmpData_buildSnpMerIndexInfo;
			// tmpData_buildSnpMerIndexInfo = NULL;
			//cout << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
			tmpData_log_phase2_ofs << "tmpData_buildSNPmerIndex_result_str: " << tmpData_buildSNPmerIndex_result_str << endl;
			tmpData_log_phase2_ofs << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
		}
		// loading SNPmers
		Index_Info* tmpData_indexInfo_SNP_learned = new Index_Info();
		char* tmpData_chrom_SNP_learned;
	    unsigned int* tmpData_sa_SNP_learned; 
	    BYTE* tmpData_lcpCompress_SNP_learned;
	    unsigned int* tmpData_childTab_SNP_learned;
	    BYTE* tmpData_verifyChild_SNP_learned;
	    if(tmpData_snpLearned_success_bool)
	    {
		    string tmpData_indexStr_SNP_learned = tmpData_dir_SNPmer_index_learned + "/";
			string tmpData_SA_file_SNP_learned = tmpData_indexStr_SNP_learned + "_SA"; 
			ifstream tmpData_SA_file_ifs_SNP_learned(tmpData_SA_file_SNP_learned.c_str(),ios::binary); 
			string tmpData_lcpCompress_file_SNP_learned = tmpData_indexStr_SNP_learned + "_lcpCompress"; 
			ifstream tmpData_lcpCompress_file_ifs_SNP_learned(tmpData_lcpCompress_file_SNP_learned.c_str(),ios::binary);
			string tmpData_childTab_file_SNP_learned = tmpData_indexStr_SNP_learned + "_childTab"; 
			ifstream tmpData_childTab_file_ifs_SNP_learned(tmpData_childTab_file_SNP_learned.c_str(),ios::binary);
			string tmpData_verifyChild_file_SNP_learned = tmpData_indexStr_SNP_learned + "_detChild"; 
			ifstream tmpData_verifyChild_file_ifs_SNP_learned(tmpData_verifyChild_file_SNP_learned.c_str(),ios::binary);	
			string tmpData_chrom_bit_file_SNP_learned = tmpData_indexStr_SNP_learned + "_chrom"; 
			ifstream tmpData_chrom_bit_file_ifs_SNP_learned(tmpData_chrom_bit_file_SNP_learned.c_str(),ios::binary);
			string tmpData_parameter_file_SNP_learned = tmpData_indexStr_SNP_learned + "_parameter"; 
			ifstream tmpData_parameter_file_ifs_SNP_learned(tmpData_parameter_file_SNP_learned.c_str(),ios::binary);

			tmpData_indexInfo_SNP_learned->initiate(tmpData_parameter_file_ifs_SNP_learned, tmpData_log_phase2_ofs);
			tmpData_chrom_SNP_learned = (char*)malloc((tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(char)); 
			tmpData_chrom_bit_file_ifs_SNP_learned.read((char*)tmpData_chrom_SNP_learned, (tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(char)); 
			tmpData_indexInfo_SNP_learned->readGenome(tmpData_chrom_SNP_learned);
			tmpData_indexInfo_SNP_learned->initiate();	
			tmpData_indexInfo_SNP_learned->initiateChrNameIndexArray(1000);
			tmpData_sa_SNP_learned = (unsigned int*)malloc((tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int)); 
			tmpData_SA_file_ifs_SNP_learned.read((char*)tmpData_sa_SNP_learned, (tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_lcpCompress_SNP_learned = (BYTE*)malloc((tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_lcpCompress_file_ifs_SNP_learned.read((char*)tmpData_lcpCompress_SNP_learned, (tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE));	 
			tmpData_childTab_SNP_learned = (unsigned int*)malloc((tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int)); 
			tmpData_childTab_file_ifs_SNP_learned.read((char*)tmpData_childTab_SNP_learned, (tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_verifyChild_SNP_learned = (BYTE*)malloc((tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_verifyChild_file_ifs_SNP_learned.read((char*)tmpData_verifyChild_SNP_learned, (tmpData_indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE));

			tmpData_SA_file_ifs_SNP_learned.close();
			tmpData_lcpCompress_file_ifs_SNP_learned.close();
			tmpData_childTab_file_ifs_SNP_learned.close();
			tmpData_verifyChild_file_ifs_SNP_learned.close();
			tmpData_chrom_bit_file_ifs_SNP_learned.close();
			tmpData_parameter_file_ifs_SNP_learned.close();	
			//cout << "Self learned SNPmer index files loaded" << endl;
			//tmpData_log_phase2_ofs << "Self learned SNPmer index files loaded" << endl;	    	
	    }
		///////////////////////////////////////////   build SNPmer index for phase2 with provided SNPs    ////////////////////////////////////////////
		Index_Info* tmpData_indexInfo_SNP = new Index_Info();
		char *tmpData_chrom_SNP;
	    unsigned int *tmpData_sa_SNP; 
	    BYTE *tmpData_lcpCompress_SNP;
	    unsigned int *tmpData_childTab_SNP;
	    BYTE *tmpData_verifyChild_SNP;
	    string tmpData_SNPfilePath;		
		if(segMap2SNPmer_phase2_bool && (!segMap2SNPmer_phase1_bool))
		{	
			string tmpData_SNPfilePath = tmpBatchManagerInfo.returnSNPfile_withIndexInVec(tmpData);
			string tmpData_SNP_seq_index_path = tmpBatchManagerInfo.returnSNPmerIndexFolderPath_withIndexInVec(tmpData);
			//cout << "SNPfilePath: " << endl << tmpData_SNPfilePath << endl;
			tmpData_log_phase2_ofs << "SNPfilePath: " << endl << tmpData_SNPfilePath << endl;
			// Xinan: for now, do not replace ref bases with alternate bases. In the future, 
			// MPS3 should incorporate SNPs when doing sequence matching and canonical splice site detection
			indexInfo->insertSNP2chromStr(tmpData_SNPfilePath, tmpData_log_phase2_ofs);
			//cout << "start to load indexes" << endl;
			string tmpData_indexStr_SNP = tmpData_SNP_seq_index_path;
			//cout << "tmpData_indexStr_SNP: " << endl << tmpData_indexStr_SNP << endl;
			tmpData_log_phase2_ofs << "tmpData_indexStr_SNP: " << endl << tmpData_indexStr_SNP << endl;
			tmpData_indexStr_SNP.append("/");
			string tmpData_SA_file_SNP = tmpData_indexStr_SNP + "_SA"; 
			ifstream tmpData_SA_file_ifs_SNP(tmpData_SA_file_SNP.c_str(),ios::binary); 
			string tmpData_lcpCompress_file_SNP = tmpData_indexStr_SNP + "_lcpCompress"; 
			ifstream tmpData_lcpCompress_file_ifs_SNP(tmpData_lcpCompress_file_SNP.c_str(),ios::binary);
			string tmpData_childTab_file_SNP = tmpData_indexStr_SNP + "_childTab"; 
			ifstream tmpData_childTab_file_ifs_SNP(tmpData_childTab_file_SNP.c_str(),ios::binary);
			string tmpData_verifyChild_file_SNP = tmpData_indexStr_SNP + "_detChild"; 
			ifstream tmpData_verifyChild_file_ifs_SNP(tmpData_verifyChild_file_SNP.c_str(),ios::binary);	
			string tmpData_chrom_bit_file_SNP = tmpData_indexStr_SNP + "_chrom"; 
			ifstream tmpData_chrom_bit_file_ifs_SNP(tmpData_chrom_bit_file_SNP.c_str(),ios::binary);
			string tmpData_parameter_file_SNP = tmpData_indexStr_SNP + "_parameter"; 
			ifstream tmpData_parameter_file_ifs_SNP(tmpData_parameter_file_SNP.c_str(),ios::binary);

			tmpData_indexInfo_SNP->initiate(tmpData_parameter_file_ifs_SNP, tmpData_log_phase2_ofs);
			tmpData_chrom_SNP = (char*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
			tmpData_chrom_bit_file_ifs_SNP.read((char*)tmpData_chrom_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
			tmpData_indexInfo_SNP->readGenome(tmpData_chrom_SNP);
			tmpData_indexInfo_SNP->initiate();	
			tmpData_indexInfo_SNP->initiateChrNameIndexArray(1000);
		    tmpData_sa_SNP = (unsigned int*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); 
		    tmpData_SA_file_ifs_SNP.read((char*)tmpData_sa_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_lcpCompress_SNP = (BYTE*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_lcpCompress_file_ifs_SNP.read((char*)tmpData_lcpCompress_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));	 
			tmpData_childTab_SNP = (unsigned int*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); 
			tmpData_childTab_file_ifs_SNP.read((char*)tmpData_childTab_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
			tmpData_verifyChild_SNP = (BYTE*)malloc((tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); 
			tmpData_verifyChild_file_ifs_SNP.read((char*)tmpData_verifyChild_SNP, (tmpData_indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
			//cout << "SyntheticSNPtransSeq index files loaded" << endl;
		}	
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////   Do mapping on one end unmapped Reads    ///////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase2_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts for dataSet " << tmpData + 1 << endl;
		//tmpData_settings_phase2_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts for dataSet " << tmpData + 1 << endl;
		
		string tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped = tmpData_outputDirStr + "/logs/phase2_log/fixOneEndUnmapped";		
	   	string tmpData_inputLogStr_phase2_fixOneEndUnmapped = tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped + "/input.log";
	   	ofstream tmpData_input_log_ofs_phase2_fixOneEndUnmapped(tmpData_inputLogStr_phase2_fixOneEndUnmapped.c_str());
	   	string tmpData_outputLogStr_phase2_fixOneEndUnmapped = tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped + "/output.log";
	   	ofstream tmpData_output_log_ofs_phase2_fixOneEndUnmapped(tmpData_outputLogStr_phase2_fixOneEndUnmapped.c_str());
	   	string tmpData_mappingLogStr_phase2_fixOneEndUnmapped = tmpData_outputDirStr_logs_phase2_fixOneEndUnmapped + "/mapping.log";
	   	ofstream tmpData_mapping_log_ofs_phase2_fixOneEndUnmapped(tmpData_mappingLogStr_phase2_fixOneEndUnmapped.c_str());  

		string tmpData_mkdirOutputCommand_phase2 = "mkdir -p " + tmpData_outputDirStr + "/phase2_output";
		system(tmpData_mkdirOutputCommand_phase2.c_str());
		string tmpData_OutputSamFile_oneEndMapped = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam";
		ofstream tmpData_OutputSamFile_oneEndMapped_ofs(tmpData_OutputSamFile_oneEndMapped.c_str());
		string tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
		ofstream tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs(tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());
		string tmpData_OutputSamFile_oneEndMapped_unpair = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.unpaired.sam";
		ofstream tmpData_OutputSamFile_oneEndMapped_unpair_ofs(tmpData_OutputSamFile_oneEndMapped_unpair.c_str());
		string tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.incompletePair.alignInfo";
		ofstream tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo_ofs(tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo);
		// string tmpData_OutputSamFile_oneEndMapped_alignInfo = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam_alignInfo";
		// ofstream tmpData_OutputSamFile_oneEndMapped_alignInfo_ofs(tmpData_OutputSamFile_oneEndMapped_alignInfo.c_str());	

		AlignInfoInput_Array_Queue* tmpData_alignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixOneEndUnmapped_Array_Queue* tmpData_fixOneEndUnmappedResultQueue = new Result_FixOneEndUnmapped_Array_Queue();
			
		string tmpData_oneEndMappedFileStr = tmpData_outputDirStr + "/phase1_output/oneEndUnmapped/oneEndUnmapped.alignInfo";
		ifstream tmpData_inputRecord_ifs(tmpData_oneEndMappedFileStr.c_str());		
		bool tmpData_endOfFile_bool = false;
		bool tmpData_endOfProcessing_bool = false;
		int tmpData_tmpInputReadNumInBatchArray_fixOneEndUnmapped = normalRecordNum_fixOneEndUnmapped;//ReadNumInReadArray_FixOneEndUnmapped;
		int tmpData_tmpInputTimeWeight_fixOneEndUnmapped = InputTimeWeight_FixOneEndUnmapped;
		int tmpData_tmpOutputTimeWeight_fixOneEndUnmapped = OutputTimeWeight_FixOneEndUnmapped;

		if(separateThreadForIO_bool)
		{	
			omp_set_num_threads(2);
			omp_set_nested(1);
			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					io_stage_fixOneEndUnmapped_separateThreadForIO(
						tmpData_inputRecord_ifs, tmpData_alignInfoInputQueue, tmpData_fixOneEndUnmappedResultQueue, tmpData_endOfFile_bool, 
						tmpData_endOfProcessing_bool, tmpData_tmpInputReadNumInBatchArray_fixOneEndUnmapped, tmpData_tmpInputTimeWeight_fixOneEndUnmapped,
						tmpData_tmpOutputTimeWeight_fixOneEndUnmapped, tmpData_log_phase2_ofs, tmpData_OutputSamFile_oneEndMapped_ofs,
						tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo_ofs, tmpData_OutputSamFile_oneEndMapped_unpair_ofs, 
						tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,
						tmpData_input_log_ofs_phase2_fixOneEndUnmapped, tmpData_output_log_ofs_phase2_fixOneEndUnmapped);
					#pragma omp section
					process_stage_fixOneEndUnmapped_separateThreadForIO(
						tmpData_alignInfoInputQueue, tmpData_fixOneEndUnmappedResultQueue, tmpData_endOfFile_bool, tmpData_endOfProcessing_bool, 
						threads_num-1,
						fasta_or_fastq_bool, tmpBatchManagerInfo.statsInfoVec[tmpData], secondLevelChrom, secondLevelSa, 
						secondLevelLcpCompress, secondLevelChildTab,
						secondLevelDetChild, indexInfo, Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, 
						Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq,
						tmpData_mapping_log_ofs_phase2_fixOneEndUnmapped, SE_or_PE_bool,
						outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
						// provided SNP
						tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, 
						tmpData_indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool, 
						// learned SNP
						tmpData_sa_SNP_learned, tmpData_lcpCompress_SNP_learned, tmpData_childTab_SNP_learned, tmpData_chrom_SNP_learned, tmpData_verifyChild_SNP_learned, 
						tmpData_indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, tmpData_snpLearned_success_bool);
				}
			}
		}
		else
		{
			io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(
				tmpData_inputRecord_ifs, tmpData_alignInfoInputQueue, tmpData_fixOneEndUnmappedResultQueue,	tmpData_tmpInputReadNumInBatchArray_fixOneEndUnmapped,
				tmpData_log_phase2_ofs, tmpData_OutputSamFile_oneEndMapped_ofs, tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo_ofs, 
				tmpData_OutputSamFile_oneEndMapped_unpair_ofs,
				tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs, tmpData_input_log_ofs_phase2_fixOneEndUnmapped,
				tmpData_output_log_ofs_phase2_fixOneEndUnmapped, threads_num, fasta_or_fastq_bool, tmpBatchManagerInfo.statsInfoVec[tmpData], secondLevelChrom,
				secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
				Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
				MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, tmpData_mapping_log_ofs_phase2_fixOneEndUnmapped, SE_or_PE_bool,
				outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
				// provided SNP
				tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, 
				tmpData_indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool,
				// learned SNP
				tmpData_sa_SNP_learned, tmpData_lcpCompress_SNP_learned, tmpData_childTab_SNP_learned, tmpData_chrom_SNP_learned, tmpData_verifyChild_SNP_learned, 
				tmpData_indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, tmpData_snpLearned_success_bool);
		}
		tmpData_alignInfoInputQueue->free(); delete tmpData_alignInfoInputQueue; tmpData_alignInfoInputQueue = NULL;
		tmpData_fixOneEndUnmappedResultQueue->free(); delete tmpData_fixOneEndUnmappedResultQueue; tmpData_fixOneEndUnmappedResultQueue = NULL;
		tmpData_inputRecord_ifs.close();

		tmpData_input_log_ofs_phase2_fixOneEndUnmapped.close();
		tmpData_output_log_ofs_phase2_fixOneEndUnmapped.close();
		tmpData_mapping_log_ofs_phase2_fixOneEndUnmapped.close();

		tmpData_OutputSamFile_oneEndMapped_ofs.close();
		tmpData_OutputSamFile_oneEndMapped_unpair_ofs.close();
		//tmpData_OutputSamFile_oneEndMapped_alignInfo_ofs.close();
		tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo_ofs.close();
		tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs.close();

		// copy files in tmpData_outputSamfile_oneEndMapped_incompletePair_alignInfo to "/phase1_output/incomplete/incomplete.alignInfo"
		string cp_targetFile = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo";
		string cp_sourceFile = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.incompletePair.alignInfo";
		string tmp_cp_cmd = "cat " + cp_sourceFile + " >> " + cp_targetFile;
		system(tmp_cp_cmd.c_str());

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... fixing oneEndUnmapped reads ends for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase2_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads ends for dataSet " << tmpData + 1 << endl;
		//tmpData_log_phase2_ofs << endl << "**********************************************************************************" << endl;
		////////////////////////////////              #BEGIN# sam 2 junc         /////////////////////////////////////////////
		////////////////////////////////              #BEGIN# sam 2 junc         /////////////////////////////////////////////
		////////////////////////////////              #BEGIN# sam 2 junc         /////////////////////////////////////////////
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo_tmpData = new AlignInferJunctionHash_Info();
		SJhash_Info* SJ_tmpData = new SJhash_Info();
		bool spliceJunctionHashExists_tmpData = false;
		if(separatingSplicingContext_bool)
		{	
			nowtime = time(NULL);
			local = localtime(&nowtime);
			tmpTimeStr = asctime(local);
			tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
			tmpData_log_phase2_ofs << endl << tmpTimeStr << "... start to build separate spliceJunction Hash" << endl;
			// generating alignInferJunctionHashInfo_tmpData
			alignInferJunctionHashInfo_tmpData->initiateAlignInferJunctionInfo(chromNum);		
			string juncfile_alignInferHash_tmpData = tmpData_outputDirStr + "/phase2_output/inter.junc.alignInferHash";
			vector<string> tmp2generateJuncAlignmentFileVec_tmpData;
			if(SE_or_PE_bool)
			{	
				string tmpData_tmpAlignCompleteRead_SE = tmpData_outputDirStr + "/phase1_output/completePair/complete_SE.sam";
				string tmpData_tmpAlignIncomplete_SE_SAM = tmpData_outputDirStr +  "/phase1_output/incomplete/incomplete_SE.sam";
				tmp2generateJuncAlignmentFileVec_tmpData.push_back(tmpData_tmpAlignCompleteRead_SE);
				tmp2generateJuncAlignmentFileVec_tmpData.push_back(tmpData_tmpAlignIncomplete_SE_SAM);		
			}
			else
			{
				string tmpData_tmpAlignCompleteRead = tmpData_outputDirStr + "/phase1_output/completePair/completePair.sam";
				string tmpData_tmpAlignIncompletePair_SAM = tmpData_outputDirStr + "/phase1_output/incomplete/incompletePair.sam";
				string tmpData_OutputSamFile_oneEndMapped = tmpData_outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam";
				tmp2generateJuncAlignmentFileVec_tmpData.push_back(tmpData_tmpAlignCompleteRead);
				tmp2generateJuncAlignmentFileVec_tmpData.push_back(tmpData_tmpAlignIncompletePair_SAM);
				tmp2generateJuncAlignmentFileVec_tmpData.push_back(tmpData_OutputSamFile_oneEndMapped);
			}
			AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec_tmpData = new AlignInferJunctionHash_Info_Vec();
			alignInferJunctionHashInfoVec_tmpData->initiateAlignInferJunctionHashInfoVec(threads_num, chromNum);
			alignInferJunctionHashInfoVec_tmpData->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(tmp2generateJuncAlignmentFileVec_tmpData, 
				indexInfo, threads_num, tmpData_log_phase2_ofs);
			alignInferJunctionHashInfoVec_tmpData->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(alignInferJunctionHashInfo_tmpData, indexInfo);
			alignInferJunctionHashInfo_tmpData->outputAlignInferInfoHashInfo_chrNamePos_supportNum(indexInfo, juncfile_alignInferHash_tmpData);
			alignInferJunctionHashInfoVec_tmpData->freeMemory();
			delete alignInferJunctionHashInfoVec_tmpData;	

			// generating SJ_tmpData
			SJ_tmpData->initiateAreaAndStringHash(chromNum);
			int junctionNum_tmpData = 0;
			int junctionNum_in_alignInferJuncHashInfo_tmpData = 0;
			int junctionNum_in_annotation_tmpData = 0;			
			string entryString_tmpData;
			int tabLocation1_tmpData, tabLocation2_tmpData, tabLocation3_tmpData, tabLocation4_tmpData, tabLocation5_tmpData;
			//char entry[500];
			int chrInt_tmpData;
			int spliceStartPos_tmpData;
			int spliceEndPos_tmpData;
			string chrIntString_tmpData;
			string spliceStartPosString_tmpData;
			string spliceEndPosString_tmpData;
			//////////////////////  string hash /////////////////////////////////////////
			if(!Do_annotation_only_bool)
			{	
		    	//tmpData_log_phase2_ofs << "start to load SJs in alignments" << endl;
	    		//cout << "start to load SJs in alignments" << endl;
				alignInferJunctionHashInfo_tmpData->convert2SJhashInfo(SJ_tmpData, indexInfo);
				junctionNum_in_alignInferJuncHashInfo_tmpData = alignInferJunctionHashInfo_tmpData->returnAlignInferInfoVecSize();
			}
			if(annotation_provided_bool)
			{	
				// loading SJs in annotation file (if provided)
		    	//tmpData_log_phase2_ofs << "start to load SJs in annotation" << endl;
	    		//cout << "start to load SJs in annotation" << endl;
				ifstream annotatedSJ_ifs_tmpData(annotation_file_path.c_str());
				while(!annotatedSJ_ifs_tmpData.eof())
				{
					getline(annotatedSJ_ifs_tmpData, entryString_tmpData);
					if(entryString_tmpData == "")
						break;
					junctionNum_in_annotation_tmpData ++;
					//entryString = entry;
					tabLocation1_tmpData = entryString_tmpData.find('\t', 0);
					tabLocation2_tmpData = entryString_tmpData.find('\t', tabLocation1_tmpData + 1);
					tabLocation3_tmpData = entryString_tmpData.find('\t', tabLocation2_tmpData + 1);
					chrIntString_tmpData = entryString_tmpData.substr(0, tabLocation1_tmpData);
					spliceStartPosString_tmpData = entryString_tmpData.substr(tabLocation1_tmpData + 1, tabLocation2_tmpData - tabLocation1_tmpData - 1);
					if(tabLocation3_tmpData == string::npos)
						spliceEndPosString_tmpData = entryString_tmpData.substr(tabLocation2_tmpData + 1);
					else
						spliceEndPosString_tmpData = entryString_tmpData.substr(tabLocation2_tmpData + 1, tabLocation3_tmpData - tabLocation2_tmpData - 1);
					chrInt_tmpData = indexInfo->convertStringToInt(chrIntString_tmpData);
					if(chrInt_tmpData >= 0)
					{	
						spliceStartPos_tmpData = atoi(spliceStartPosString_tmpData.c_str());
						spliceEndPos_tmpData = atoi(spliceEndPosString_tmpData.c_str());	
						SJ_tmpData->insert2AreaAndStringHash(chrInt_tmpData, spliceStartPos_tmpData, spliceEndPos_tmpData, indexInfo);
					}
				}
				annotatedSJ_ifs_tmpData.close();
			}
			junctionNum_tmpData = junctionNum_in_alignInferJuncHashInfo_tmpData + junctionNum_in_annotation_tmpData;
			if(junctionNum_tmpData == 0)
				spliceJunctionHashExists_tmpData = false;
			else
				spliceJunctionHashExists_tmpData = true;
			tmpData_log_phase2_ofs << "finish building spliceJunction Hash" << endl;
			tmpData_log_phase2_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum_tmpData << endl;
			tmpData_log_phase2_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			tmpTimeStr = asctime(local);
			tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
			tmpData_log_phase2_ofs << endl << tmpTimeStr << "... end of building separate spliceJunction Hash" << endl;
		}
		////////////////////////////////               #END# sam 2 junc         /////////////////////////////////////////////
		////////////////////////////////               #END# sam 2 junc         /////////////////////////////////////////////
		////////////////////////////////               #END# sam 2 junc         /////////////////////////////////////////////
		//////////////////////////////   #BEGIN#    Do mapping on unfixed head/tail Reads     ///////////////////////////////
		//////////////////////////////   #BEGIN#    Do mapping on unfixed head/tail Reads     ///////////////////////////////
		//////////////////////////////   #BEGIN#    Do mapping on unfixed head/tail Reads     ///////////////////////////////
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... fixing unfixed-head/tail reads starts for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase2_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads starts for dataSet " << tmpData + 1 << endl;
		//tmpData_log_phase2_ofs << endl << "**********************************************************************************" << endl;

	  	string tmpData_outputDirStr_logs_phase2_fixHeadTail = tmpData_outputDirStr + "/logs/phase2_log/fixHeadTail";
	   	string tmpData_inputLogStr_phase2_fixHeadTail = tmpData_outputDirStr_logs_phase2_fixHeadTail + "/input.log";
	  	ofstream tmpData_input_log_ofs_phase2_fixHeadTail(tmpData_inputLogStr_phase2_fixHeadTail.c_str());
		string tmpData_outputLogStr_phase2_fixHeadTail = tmpData_outputDirStr_logs_phase2_fixHeadTail + "/output.log";
	   	ofstream tmpData_output_log_ofs_phase2_fixHeadTail(tmpData_outputLogStr_phase2_fixHeadTail.c_str());
	   	string tmpData_mappingLogStr_phase2_fixHeadTail = tmpData_outputDirStr_logs_phase2_fixHeadTail + "/mapping.log";
		ofstream tmpData_mapping_log_ofs_phase2_fixHeadTail(tmpData_mappingLogStr_phase2_fixHeadTail.c_str());

		string tmpData_OutputSamFile_fixHeadTail_complete_pair = tmpData_outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_complete_pair_ofs(tmpData_OutputSamFile_fixHeadTail_complete_pair.c_str());
		string tmpData_OutputSamFile_fixHeadTail_incomplete_pair = tmpData_outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_incomplete_pair_ofs(tmpData_OutputSamFile_fixHeadTail_incomplete_pair.c_str());
		string tmpData_OutputSamFile_fixHeadTail_complete_unpair = tmpData_outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_complete_unpair_ofs(tmpData_OutputSamFile_fixHeadTail_complete_unpair.c_str());
		string tmpData_OutputSamFile_fixHeadTail_incomplete_unpair = tmpData_outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_incomplete_unpair_ofs(tmpData_OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	
		string tmpData_OutputSamFile_fixHeadTail_pair_lowScore = tmpData_outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_pair_lowScore_ofs(tmpData_OutputSamFile_fixHeadTail_pair_lowScore.c_str());	
		string tmpData_OutputSamFile_fixHeadTail_complete_SE = tmpData_outputDirStr + "/phase2_output/fixHeadTail_complete_SE.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_complete_SE_ofs(tmpData_OutputSamFile_fixHeadTail_complete_SE.c_str());
		string tmpData_OutputSamFile_fixHeadTail_incomplete_SE = tmpData_outputDirStr + "/phase2_output/fixHeadTail_incomplete_SE.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_incomplete_SE_ofs(tmpData_OutputSamFile_fixHeadTail_incomplete_SE.c_str());
		string tmpData_OutputSamFile_fixHeadTail_lowScore_SE = tmpData_outputDirStr + "/phase2_output/fixHeadTail_lowScore_SE.sam";
		ofstream tmpData_OutputSamFile_fixHeadTail_lowScore_SE_ofs(tmpData_OutputSamFile_fixHeadTail_lowScore_SE.c_str());
		
		string headTailSoftClippingFile;
		if(SE_or_PE_bool)
			headTailSoftClippingFile = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete_SE.alignInfo";
		else
			headTailSoftClippingFile = tmpData_outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo";
		AlignInfoInput_Array_Queue* tmpData_fixHeadTailAlignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixHeadTail_Array_Queue* tmpData_fixHeadTailResultQueue = new Result_FixHeadTail_Array_Queue();
		ifstream tmpData_inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());
		//int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool tmpData_endOfFile_bool_unfixedHeadTail = false;
		bool tmpData_endOfProcessing_bool_unfixedHeadTail = false;
		int tmpData_tmpInputReadNumInBatchArray_fixHeadTail = normalRecordNum_fixHeadTail;//ReadNumInReadArray_FixHeadTail;
		int tmpData_tmpInputTimeWeight_fixHeadTail = InputTimeWeight_FixHeadTail;
		int tmpData_tmpOutputTimeWeight_fixHeadTail = OutputTimeWeight_FixHeadTail;

		if(separateThreadForIO_bool)
		{	
			omp_set_num_threads(2);
			omp_set_nested(1);
			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					io_stage_fixHeadTail_separateThreadForIO(
						tmpData_inputUnfixedHeadTailRecord_ifs, tmpData_fixHeadTailAlignInfoInputQueue, tmpData_fixHeadTailResultQueue, 
						tmpData_endOfFile_bool_unfixedHeadTail, tmpData_endOfProcessing_bool_unfixedHeadTail, 
						tmpData_tmpInputReadNumInBatchArray_fixHeadTail, tmpData_tmpInputTimeWeight_fixHeadTail,
						tmpData_tmpOutputTimeWeight_fixHeadTail, tmpData_log_phase2_ofs, tmpData_OutputSamFile_fixHeadTail_complete_pair_ofs,
						tmpData_OutputSamFile_fixHeadTail_incomplete_pair_ofs, tmpData_OutputSamFile_fixHeadTail_complete_unpair_ofs,
						tmpData_OutputSamFile_fixHeadTail_incomplete_unpair_ofs, tmpData_OutputSamFile_fixHeadTail_pair_lowScore_ofs,
						tmpData_input_log_ofs_phase2_fixHeadTail, tmpData_output_log_ofs_phase2_fixHeadTail);
					#pragma omp section
					process_stage_fixHeadTail_separateThreadForIO(
						tmpData_fixHeadTailAlignInfoInputQueue, tmpData_fixHeadTailResultQueue, tmpData_endOfFile_bool_unfixedHeadTail, 
						tmpData_endOfProcessing_bool_unfixedHeadTail, threads_num-1, fasta_or_fastq_bool, tmpBatchManagerInfo.statsInfoVec[tmpData], 
						secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
						////////////////////////////////////////////
						shared_SJ, SJ_tmpData,
						////////////////////////////////////////////
						Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, 
						annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
						//////////////////////////////////////////////////////////////////////////
						spliceJunctionHashExists_shared, spliceJunctionHashExists_tmpData,
						//////////////////////////////////////////////////////////////////////////
						Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
						Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain,
						tmpData_mapping_log_ofs_phase2_fixHeadTail, SE_or_PE_bool,
						outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
						// provided SNPs
						tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, 
						tmpData_indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool,
						// learned SNPs
						tmpData_sa_SNP_learned, tmpData_lcpCompress_SNP_learned, tmpData_childTab_SNP_learned, tmpData_chrom_SNP_learned, tmpData_verifyChild_SNP_learned, 
						tmpData_indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, tmpData_snpLearned_success_bool);
				}
			}
		}
		else
		{
			io_process_fixHeadTail_allThreadsSharedByBothStage(
				tmpData_inputUnfixedHeadTailRecord_ifs, tmpData_fixHeadTailAlignInfoInputQueue, tmpData_fixHeadTailResultQueue, 
				tmpData_tmpInputReadNumInBatchArray_fixHeadTail, tmpData_log_phase2_ofs,
				tmpData_OutputSamFile_fixHeadTail_complete_pair_ofs, tmpData_OutputSamFile_fixHeadTail_incomplete_pair_ofs,
				tmpData_OutputSamFile_fixHeadTail_complete_unpair_ofs, tmpData_OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
				tmpData_OutputSamFile_fixHeadTail_pair_lowScore_ofs, tmpData_input_log_ofs_phase2_fixHeadTail,
				tmpData_output_log_ofs_phase2_fixHeadTail, threads_num, fasta_or_fastq_bool, tmpBatchManagerInfo.statsInfoVec[tmpData],
				secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab,
				secondLevelDetChild, indexInfo, 
				////////////////////////////////////////////
				shared_SJ, SJ_tmpData,
				////////////////////////////////////////////
				Do_extendHeadTail_fixHeadTail, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, 
				checkQualSeqForShortAnchorSeqToTargetMap,
				//////////////////////////////////////////////////////////////////////////
				spliceJunctionHashExists_shared, spliceJunctionHashExists_tmpData,
				//////////////////////////////////////////////////////////////////////////
				Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain,
				tmpData_mapping_log_ofs_phase2_fixHeadTail, SE_or_PE_bool,
				outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
				// provided SNPs
				tmpData_sa_SNP, tmpData_lcpCompress_SNP, tmpData_childTab_SNP, tmpData_chrom_SNP, tmpData_verifyChild_SNP, 
				tmpData_indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool, 
				// learned SNPs
				tmpData_sa_SNP_learned, tmpData_lcpCompress_SNP_learned, tmpData_childTab_SNP_learned, tmpData_chrom_SNP_learned, tmpData_verifyChild_SNP_learned, 
				tmpData_indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, tmpData_snpLearned_success_bool);
		}

		delete alignInferJunctionHashInfo_tmpData; alignInferJunctionHashInfo_tmpData = NULL;
		delete SJ_tmpData; SJ_tmpData = NULL;
		tmpData_fixHeadTailAlignInfoInputQueue->free(); delete tmpData_fixHeadTailAlignInfoInputQueue; tmpData_fixHeadTailAlignInfoInputQueue = NULL;
		tmpData_fixHeadTailResultQueue->free(); delete tmpData_fixHeadTailResultQueue; tmpData_fixHeadTailResultQueue = NULL;
		tmpData_inputUnfixedHeadTailRecord_ifs.close();		
		tmpData_input_log_ofs_phase2_fixHeadTail.close();
		tmpData_output_log_ofs_phase2_fixHeadTail.close();
		tmpData_mapping_log_ofs_phase2_fixHeadTail.close();
		tmpData_OutputSamFile_fixHeadTail_complete_pair_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_complete_SE_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_incomplete_SE_ofs.close();
		tmpData_OutputSamFile_fixHeadTail_lowScore_SE_ofs.close();

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... fixing unfixed-head/tail reads ends for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase2_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads ends for dataSet " << tmpData + 1 << endl;

		int collectedInPhase1_readTotalNum = tmpBatchManagerInfo.returnTotalReadNum_withIndexInVec(tmpData);
		if(SE_or_PE_bool)
		{
			//tmpBatchManagerInfo.statsInfoVec[tmpData]->getPhase1Stats_SE();
			tmpBatchManagerInfo.statsInfoVec[tmpData]->getFixHeadTailStats_SE();
			tmpBatchManagerInfo.statsInfoVec[tmpData]->outputAllStats_SE_fixHeadTail(tmpData_stats_phase2_ofs, collectedInPhase1_readTotalNum);
			tmpBatchManagerInfo.statsInfoVec[tmpData]->outputFinalStats_SE(tmpData_stats_phase2_ofs, collectedInPhase1_readTotalNum);
		}	
		else
		{	
			tmpBatchManagerInfo.statsInfoVec[tmpData]->getPhase1Stats();
			tmpBatchManagerInfo.statsInfoVec[tmpData]->getFixUnpairedStats();
			tmpBatchManagerInfo.statsInfoVec[tmpData]->getFixHeadTailStats();
			//tmpBatchManagerInfo.statsInfoVec[tmpData]->outputAllStats(log_ofs, readTotalNum);
			tmpBatchManagerInfo.statsInfoVec[tmpData]->outputAllStats(tmpData_stats_phase2_ofs, Do_Phase1_Only, collectedInPhase1_readTotalNum);
			tmpBatchManagerInfo.statsInfoVec[tmpData]->outputFinalStats(tmpData_stats_phase2_ofs, Do_Phase1_Only, collectedInPhase1_readTotalNum);	
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		//cout << endl << tmpTimeStr << "... start to prepare for final output files for dataSet " << tmpData + 1 << endl;
		tmpData_log_phase2_ofs << endl << tmpTimeStr << "... start to prepare for final output files for dataSet " << tmpData + 1 << endl;
   
		if(tmpData_snpLearned_success_bool)
		{
			free(tmpData_sa_SNP_learned);
			free(tmpData_chrom_SNP_learned);
		    free(tmpData_lcpCompress_SNP_learned);
	    	free(tmpData_childTab_SNP_learned);
	    	free(tmpData_verifyChild_SNP_learned);
    	}
    	tmpData_sa_SNP_learned = NULL;
    	tmpData_chrom_SNP_learned = NULL;
    	tmpData_lcpCompress_SNP_learned = NULL;
    	tmpData_childTab_SNP_learned = NULL;
    	tmpData_verifyChild_SNP_learned = NULL;
    	delete tmpData_indexInfo_SNP_learned;
    	tmpData_indexInfo_SNP_learned = NULL;

    	if(segMap2SNPmer_phase2_bool && (!segMap2SNPmer_phase1_bool))
    	{
    		free(tmpData_sa_SNP);
			free(tmpData_chrom_SNP);
		    free(tmpData_lcpCompress_SNP);
		    free(tmpData_childTab_SNP);
	    	free(tmpData_verifyChild_SNP);
    	}
    	tmpData_sa_SNP = NULL;
    	tmpData_chrom_SNP = NULL;
    	tmpData_lcpCompress_SNP = NULL;
    	tmpData_childTab_SNP = NULL;
    	tmpData_verifyChild_SNP = NULL;
    	delete tmpData_indexInfo_SNP;
    	tmpData_indexInfo_SNP = NULL;

    	tmpData_log_phase2_ofs.close();
    	tmpData_settings_phase2_ofs.close();
    	tmpData_stats_phase2_ofs.close();
	
		string finalOutputSam = tmpData_outputDirStr + "/output.sam";
		string cat_cmd;
		string tmpPhase1_tmpHeadSectionInfo = tmpData_outputDirStr + "/headSectionInfo";
		string tmpPhase1_tmpAlignCompleteRead = tmpData_outputDirStr + "/phase1_output/completePair/completePair.sam";
		string tmpPhase1_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = tmpData_outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam";
		string tmpPhase1_tmpAlignBothEndsUnmapped_lowScore = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
		string tmpPhase1_tmpAlignBothEndsUnmapped = tmpData_outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
		cat_cmd = "cat " + tmpPhase1_tmpHeadSectionInfo
			+ " " + tmpPhase1_tmpAlignCompleteRead 
			+ " " + tmpData_OutputSamFile_oneEndMapped 
			+ " " + tmpData_OutputSamFile_fixHeadTail_complete_pair
			+ " " + tmpData_OutputSamFile_fixHeadTail_incomplete_pair 
			+ " " + tmpData_OutputSamFile_fixHeadTail_complete_unpair 
			+ " " + tmpData_OutputSamFile_fixHeadTail_incomplete_unpair
			+ " " + tmpData_OutputSamFile_oneEndMapped_unpair 
			+ " " + tmpPhase1_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile 
			+ " " + tmpPhase1_tmpAlignBothEndsUnmapped_lowScore
			+ " " + tmpData_OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore 
			+ " " + tmpData_OutputSamFile_fixHeadTail_pair_lowScore
			+ " " + tmpPhase1_tmpAlignBothEndsUnmapped + " > " + finalOutputSam;
		system(cat_cmd.c_str());

		if(reportJunc_bool)
		{
			string sam2junc_folder = tmpData_outputDirStr + "/sam2junc";
			string mkdir_sam2junc_cmd = "mkdir -p " + sam2junc_folder;
			system(mkdir_sam2junc_cmd.c_str());
			string reportJunc_file = tmpData_outputDirStr + "/output.junc";
			string reportJunc_cmd = "./reportJunc_embeded " + indexStr + " " + int_to_str(threads_num) + " "
				+ sam2junc_folder + " " + finalOutputSam  + " " + reportJunc_file;
			system(reportJunc_cmd.c_str());
		}

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... end of phase2 mapping for dataset " << tmpData + 1 << endl;
		log_batch_ofs << endl << tmpTimeStr << "... end of phase2 mapping for dataset " << tmpData + 1 << endl;
		runtime_batch_ofs << endl << tmpTimeStr << "... end of phase2 mapping for dataset " << tmpData + 1 << endl;
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... end of phase2 for all datasets ..." << endl;
	log_batch_ofs << endl << tmpTimeStr << "... end of phase2 for all datasets ..." << endl;
	runtime_batch_ofs << endl << tmpTimeStr << "... end of phase2 for all datasets ..." << endl;

	annotation_ifs.close();
	delete alignInferJunctionHashInfo_shared; alignInferJunctionHashInfo_shared = NULL;
	delete shared_SJ; shared_SJ = NULL;
	delete annotationInfo; annotationInfo = NULL;
	delete indexInfo; indexInfo = NULL;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << endl << tmpTimeStr << "... all jobs done ..." << endl << endl ;  
	log_batch_ofs << endl << endl << tmpTimeStr << "... all jobs done ..." << endl << endl ;  
	runtime_batch_ofs << endl << endl << tmpTimeStr << "... all jobs done ..." << endl << endl ;  

    log_batch_ofs.close();
   	settings_batch_ofs.close();
    process_batch_ofs.close();
   	runtime_batch_ofs.close();

    return 0;
} //end main