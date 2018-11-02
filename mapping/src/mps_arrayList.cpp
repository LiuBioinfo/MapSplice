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
#include <cstdlib>
//#include "switch.h"
#include "otherProjects/indexing/utils/build_snpMerIndex.h"
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

using namespace std;  

clock_t loadGlobalIndex_start, loadGlobalIndex_end, loadGlobalIndex_cost = 0,
		loadLocalIndex_start, loadLocalIndex_end, loadLocalIndex_cost = 0,
		phase1map_start, phase1map_end, phase1map_cost = 0,
		oneEndMap_start, oneEndMap_end, oneEndMap_cost = 0,
		headTailMap_start, headTailMap_end, headTailMap_cost = 0,
		buildSNPcontext_start, buildSNPcontext_end, buildSNPcontext_cost = 0,
		buildSpliceContext_start, buildSpliceContext_end, buildSpliceContext_cost = 0,
		samReport_start, samReport_end, samReport_cost = 0,
		juncReport_start, juncReport_end, juncReport_cost = 0;

void createFolders(string& outputDirStr)
{
    string outputDirStr_logs = outputDirStr + "/logs";
    string outputDirStr_logs_phase1 = outputDirStr_logs + "/phase1_log";
    string outputDirStr_logs_phase2 = outputDirStr_logs + "/phase2_log";
    string outputDirStr_logs_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2 + "/fixOneEndUnmapped";
    string outputDirStr_logs_phase2_fixHeadTail = outputDirStr_logs_phase2 + "/fixHeadTail";
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	const int mkdir_err = system(mkdirOutputCommand.c_str());	
   	string mkdirOutputCommand_log = "mkdir -p " + outputDirStr_logs;
   	const int mkdir_log_err = system(mkdirOutputCommand_log.c_str());
   	string mkdirOutputCommand_log_phase1 = "mkdir -p " + outputDirStr_logs_phase1;
   	const int mkdir_log_phase1_err = system(mkdirOutputCommand_log_phase1.c_str());   	
   	string mkdirOutputCommand_log_phase2 = "mkdir -p " + outputDirStr_logs_phase2;
   	const int mkdir_log_phase2_err = system(mkdirOutputCommand_log_phase2.c_str());
   	string mkdirOutputCommand_log_phase2_fixOneEndUnmapped = "mkdir -p " + outputDirStr_logs_phase2_fixOneEndUnmapped;
   	const int mkdir_log_phase2_fixOneEndUnmapped_err = system(mkdirOutputCommand_log_phase2_fixOneEndUnmapped.c_str());
   	string mkdirOutputCommand_log_phase2_fixHeadTail = "mkdir -p " + outputDirStr_logs_phase2_fixHeadTail;
   	const int mkdir_log_phase2_fixHeadTail_err = system(mkdirOutputCommand_log_phase2_fixHeadTail.c_str());
	string mkdirOutputCommand_phase1 = "mkdir -p " + outputDirStr + "/phase1_output";
	const int mkdir_phase1_err = system(mkdirOutputCommand_phase1.c_str());
	string mkdirOutputCommand_repeatRegionFile = mkdirOutputCommand_phase1 + "/repeat_region";	
   	const int mkdir_repeatRegionFile_err = system(mkdirOutputCommand_repeatRegionFile.c_str());
	string mkdirOutputCommand_tmpAlignOneEndUnmapped = mkdirOutputCommand_phase1 + "/oneEndUnmapped";
	const int mkdir_tmpAlignOneEndUnmapped_err = system(mkdirOutputCommand_tmpAlignOneEndUnmapped.c_str());   	
	string mkdirOutputCommand_tmpAlignIncompletePair = mkdirOutputCommand_phase1 + "/incomplete";
	const int mkdir_tmpAlignIncompletePair_err = system(mkdirOutputCommand_tmpAlignIncompletePair.c_str());
	string mkdirOutputCommand_tmpAlignCompleteRead = mkdirOutputCommand_phase1 + "/completePair";
	const int mkdir_tmpAlignCompleteRead_err = system(mkdirOutputCommand_tmpAlignCompleteRead.c_str());	
	string mkdirOutputCommand_tmpAlignBothEndsUnmapped = mkdirOutputCommand_phase1 + "/bothEndsUnmapped";
	const int mkdir_tmpAlignBothEndsUnmapped_err = system(mkdirOutputCommand_tmpAlignBothEndsUnmapped.c_str());
   	string mkdir_SNPmer_learned = "mkdir -p " + outputDirStr + "/SNPmer_learned";
   	const int mkdir_SNPmer_learned_err = system(mkdir_SNPmer_learned.c_str());
	string mkdirOutputCommand_phase2 = "mkdir -p " + outputDirStr + "/phase2_output";
	const int mkdir_phase2_err = system(mkdirOutputCommand_phase2.c_str());
	string mkdir_sam2junc_cmd = "mkdir -p " + outputDirStr + "/sam2junc";
	const int mkdir_sam2junc_err = system(mkdir_sam2junc_cmd.c_str());	
   	if((-1 == mkdir_err)||(-1 == mkdir_log_err)||(-1 == mkdir_log_phase1_err)||(-1 == mkdir_log_phase2_err)
   		||(-1 == mkdir_log_phase2_fixOneEndUnmapped_err)||(-1 == mkdir_log_phase2_fixHeadTail_err)||(-1 == mkdir_phase1_err)
   		||(-1 == mkdir_repeatRegionFile_err)||(-1 == mkdir_tmpAlignCompleteRead_err)||(-1 == mkdir_tmpAlignOneEndUnmapped_err)
   		||(-1 == mkdir_tmpAlignBothEndsUnmapped_err)||(-1 == mkdir_tmpAlignIncompletePair_err)||(-1 == mkdir_SNPmer_learned_err)
   		||(-1 == mkdir_phase2_err)||(-1 == mkdir_sam2junc_err))
   	{
   		cout << "error in folder creating ..." << endl;
   		cout << "mkdir_err: " << mkdir_err << "\nmkdir_log_err: " << mkdir_log_err << "\nmkdir_log_phase1_err: " << mkdir_log_phase1_err
   			<< "\nmkdir_log_phase2_err: " << mkdir_log_phase2_err << "\nmkdir_log_phase2_fixOneEndUnmapped_err: " << mkdir_log_phase2_fixOneEndUnmapped_err
   			<< "\nmkdir_log_phase2_fixHeadTail_err: " << mkdir_log_phase2_fixHeadTail_err << endl;
   		cout << "mkdir_phase1_err: " << mkdir_phase1_err << "\nmkdir_repeatRegionFile_err: " << mkdir_repeatRegionFile_err 
   			<< "\nmkdir_tmpAlignCompleteRead_err: " << mkdir_tmpAlignCompleteRead_err
   			<< "\nmkdir_tmpAlignOneEndUnmapped_err: " << mkdir_tmpAlignOneEndUnmapped_err 
   			<< "\nmkdir_tmpAlignBothEndsUnmapped_err: " << mkdir_tmpAlignBothEndsUnmapped_err
   			<< "\nmkdir_tmpAlignIncompletePair_err: " << mkdir_tmpAlignIncompletePair_err << endl;   
   		cout << "mkdir_SNPmer_learned_err: " << mkdir_SNPmer_learned_err << endl;   
   		cout << "mkdir_phase2_err: " << mkdir_phase2_err << endl;
   	}
}

int main(int argc, char**argv)
{
	time_t nowtime;
	struct tm *local;
	string tmpTimeStr;

	loadGlobalIndex_start = clock();

	int normalRecordNum_1stMapping = 500000;
	int normalRecordNum_fixOneEndUnmapped = 500000;
	int normalRecordNum_fixHeadTail = 500000;

    bool checkQualSeqForShortAnchorSeqToTargetMap = false;//true;
    //cout << "Attention! checkQualSeqForShortAnchorSeqToTargetMap true or not: " << checkQualSeqForShortAnchorSeqToTargetMap << endl;
    bool checkQualSeqForReadSegSeq = false;//true;
	//cout << "Attention! checkQualSeqForReadSegSeq true or not: " << checkQualSeqForReadSegSeq << endl;
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

    //bool avoid_creating_folder_bool = optionInfo->return_avoid_creating_folder_bool();
    bool stopLearningSNP_bool = optionInfo->return_stopLearningSNP_bool();
  //FIXME: turn off it for now  
  stopLearningSNP_bool = true;

	bool segMap2SNPmer_phase2_bool_learned = (!stopLearningSNP_bool);//true;
	bool snpLearned_success_bool = segMap2SNPmer_phase2_bool_learned;
	int SNPmerLength = 201;
	int SNPlocInSyntheticSNPseq = (SNPmerLength - 1)/2;
	int SNPlocInSyntheticSNPseq_learned = (SNPmerLength - 1)/2;

	bool stopLearningSplicing_bool = optionInfo->return_stopLearningSplicing_bool();

    bool sharedThreadForIO_bool = optionInfo->sharedThreadForIO_bool;
    bool separateThreadForIO_bool = (!sharedThreadForIO_bool);

    bool extractUnmapAlignment2ReadFile_bool = optionInfo->return_extractUnmapAlignment2ReadFile_bool();
    string inputSamFilePathStr = optionInfo->returnInputSamFilePath();
    if(extractUnmapAlignment2ReadFile_bool)
    	unmappedAlignemnt2ReadFile(inputSamFilePathStr);
    bool SE_or_PE_bool = optionInfo->returnSEorPE_bool();
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   initiate log files    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
    string outputDirStr = optionInfo->outputFolder_path; //argv[3];
    createFolders(outputDirStr);
    string outputDirStr_logs = outputDirStr + "/logs";
    string outputDirStr_logs_phase1 = outputDirStr_logs + "/phase1_log";
    string outputDirStr_logs_phase2 = outputDirStr_logs + "/phase2_log";
    string outputDirStr_logs_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2 + "/fixOneEndUnmapped";
    string outputDirStr_logs_phase2_fixHeadTail = outputDirStr_logs_phase2 + "/fixHeadTail";

   	string settingsLogStr = outputDirStr_logs + "/settings.log";
   	ofstream settings_log_ofs(settingsLogStr.c_str());
   	string progressLogStr = outputDirStr_logs + "/process.log";
   	ofstream log_ofs(progressLogStr.c_str());
   	string runtimeLogStr = outputDirStr_logs + "/runtime.log";
   	ofstream runtime_log_ofs(runtimeLogStr.c_str());
   	string statsStr = outputDirStr + "/stats.txt";
   	ofstream stats_ofs(statsStr.c_str());
   	// string runTimeStr = outputDirStr + "/runtime.txt";
   	// ofstream runtime_overall_ofs(runTimeStr.c_str());
   	string inputLogStr_phase1 = outputDirStr_logs_phase1 + "/input.log";
   	ofstream input_log_ofs_phase1(inputLogStr_phase1.c_str());
   	string outputLogStr_phase1 = outputDirStr_logs_phase1 + "/output.log";
   	ofstream output_log_ofs_phase1(outputLogStr_phase1.c_str());
   	string mappingLogStr_phase1 = outputDirStr_logs_phase1 + "/mapping.log";
   	ofstream mapping_log_ofs_phase1(mappingLogStr_phase1.c_str());   	
   	string inputLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/input.log";
   	ofstream input_log_ofs_phase2_fixOneEndUnmapped(inputLogStr_phase2_fixOneEndUnmapped.c_str());
   	string outputLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/output.log";
   	ofstream output_log_ofs_phase2_fixOneEndUnmapped(outputLogStr_phase2_fixOneEndUnmapped.c_str());
   	string mappingLogStr_phase2_fixOneEndUnmapped = outputDirStr_logs_phase2_fixOneEndUnmapped + "/mapping.log";
   	ofstream mapping_log_ofs_phase2_fixOneEndUnmapped(mappingLogStr_phase2_fixOneEndUnmapped.c_str());   
   	string inputLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/input.log";
   	ofstream input_log_ofs_phase2_fixHeadTail(inputLogStr_phase2_fixHeadTail.c_str());
   	string outputLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/output.log";
   	ofstream output_log_ofs_phase2_fixHeadTail(outputLogStr_phase2_fixHeadTail.c_str());
   	string mappingLogStr_phase2_fixHeadTail = outputDirStr_logs_phase2_fixHeadTail + "/mapping.log";
   	ofstream mapping_log_ofs_phase2_fixHeadTail(mappingLogStr_phase2_fixHeadTail.c_str());
   	optionInfo->outputOptStr(settings_log_ofs);
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   switches of seperate processes    ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
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
    cout << "fusion_search_bool: " << fusion_search_bool << endl;
    cout << "report_junc_bool: " << reportJunc_bool << endl;
  	settings_log_ofs << "backSplice_search_bool: " << backSplice_search_bool << endl;
    settings_log_ofs << "fusion_search_bool: " << fusion_search_bool << endl;
    settings_log_ofs << "report_junc_bool: " << reportJunc_bool << endl;

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
	bool DoRemappingOnUnmapEndReadsBool = false;
	DoRemappingOnUnmapEndReadsBool = true;
	bool DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	DoRemappingOnUnfixedHeadTailAlignmentBool = true;
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
		settings_log_ofs << "Do_Phase1 only!" << endl;
		DoSam2JuncBool = false;
		load2ndLevelIndexBool = false;
		load2ndLevelIndexBool_compressedSize = false;
		DoRemappingOnUnmapEndReadsBool = false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = false;
	}	
	else
	{
		settings_log_ofs << "Do_Phase1_Phase2! " << endl;
		DoSam2JuncBool = true;//false;
		load2ndLevelIndexBool = true;//false;
		load2ndLevelIndexBool_compressedSize = true;//false;
		DoRemappingOnUnmapEndReadsBool = true;//false;
		DoRemappingOnUnfixedHeadTailAlignmentBool = true;//false;
	}
	int readTotalNum = 0;
	optionInfo->outputSwitchInfo(Do_Phase1_Only, outputAlignInfoAndSamForAllPairedAlignmentBool,
		removeAllIntermediateFilesBool, Do_cirRNA, outputDirectlyBool_Phase1Only, 
		normalRecordNum_1stMapping, normalRecordNum_fixOneEndUnmapped,
		normalRecordNum_fixHeadTail, Do_extendHeadTail_phase1, 
		Do_extendHeadTail_fixOneEndUnmapped, Do_extendHeadTail_fixHeadTail, 
		Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
		Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain, settings_log_ofs);
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... MPS starts ......" << endl;  
	log_ofs << endl << tmpTimeStr << "... MPS starts ......" << endl; 
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
	if(checkQualSeqForShortAnchorSeqToTargetMap && fasta_or_fastq_bool)
	{
		cout << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		log_ofs << "if checkQualSeqForShortAnchorSeqToTargetMap, fasta_or_fastq_bool must be true" << endl;
		exit(1);
	}
	/////////////////////////////////////          LOAD INDEX         ////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start to load whole genome index ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... start to load whole genome index ......" << endl; 

    //cout << "start to load preIndex ..." << endl;
    log_ofs << "start to load preIndex ..." << endl;
    string indexStr = optionInfo->global_index_file_path_prefix; //argv[6];
    string preIndexArrayPreStr = indexStr;
    string secondLevelIndexStr = optionInfo->local_index_file_path_prefix; //argv[7];
    preIndexArrayPreStr.append("/");
    indexStr.append("/");
    secondLevelIndexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));
 	log_ofs << "finish loading preIndex ..." << endl;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);
	settings_log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	settings_log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;
 	
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

	log_ofs << "start to load SA, lcpCompress, childTab, detChild ....." << endl;
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	log_ofs << "All index files loaded" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... all index loaded ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... all index loaded ......" << endl;
	runtime_log_ofs << endl << tmpTimeStr << "... all index loaded ......" << endl;	

	//////////////////////////////////////////////////
	HeaderSection_Info* headerSectionInfo = new HeaderSection_Info(indexInfo);	
	/////////////////////////////////////   start to load annotation  /////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start to load annotation file (SJs)......" << endl; 
	log_ofs << endl << tmpTimeStr << "... start to load annotation file (SJs) ......" << endl; 	
	string annotation_file_path = optionInfo->annotation_file_path; // junction files
	ifstream annotation_ifs(annotation_file_path.c_str());
	Annotation_Info* annotationInfo = new Annotation_Info();
	if(annotation_provided_bool)
		annotationInfo->initiateAndReadAnnotationFile(indexInfo, annotation_ifs);
	/////////////////////////////////////   finish loading annotation  /////////////////////////////////////	
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	/////^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   load alignInferJunctionHashInfo ************************************************///
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	bool SJalignInferHash_provided_bool = optionInfo->spliceJunctionAlignInferHash_provided_bool;
	//cout << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "start to initaite alignInferJunctionHashInfo " << endl;
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	if(SJalignInferHash_provided_bool)
	{	
		//cout << "start to insert SJ into SJmap" << endl;
		string tmpInputJuncFile = optionInfo->spliceJunctionAlignInferHash_file_path;
		alignInferJunctionHashInfo->insertJuncFromJuncFile_onlyChrNamePos(tmpInputJuncFile, indexInfo);
		//cout << "start to output SJ map" << endl;
	}
	//////////////////////////////////////////////////       finish LOADing INDEX      ////////////////////////////////////////////////////////////////
	/////////////////////////////////////////    initiating   learning SNPs from phase1 (by default)        /////////////////////////////////////////////////////
	//cout << "start to initiate learnedCandiSNPhashInfo ....." << endl;
	LearnedCandiSNPhash_Info_Vec mismatchHashInfoVec;
	int mismatchHashInfoVecSize = threads_num;
	mismatchHashInfoVec.initiateLearnedCandiSNPhashInfoVec(mismatchHashInfoVecSize, chromNum);
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
	//////////////////////////////////////////////////////////////////////
	// check if updateChrSeqWithSNPonly_avoidSegMap2SNPmer or not .... ///
	//////////////////////////////////////////////////////////////////////
	bool updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool = optionInfo->return_updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool();
	bool updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool = optionInfo->return_updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool();
	if(updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool)
	{
		segMap2SNPmer_phase1_bool = false;
		segMap2SNPmer_phase2_bool = false;		
		segMap2SNPmer_phase2_bool_learned = false;
		snpLearned_success_bool = false;
	}

	//cout << "segMap2SNPmer_phase_defined_bool: " << segMap2SNPmer_phase_defined_bool << endl;
	//cout << "segMap2SNPmer_phase1_bool: " << segMap2SNPmer_phase1_bool << endl;
	//cout << "segMap2SNPmer_phase2_bool: " << segMap2SNPmer_phase2_bool << endl;
	//cout << "segMap2SNPmer_phase2_bool_learned: " << segMap2SNPmer_phase2_bool_learned << endl;
	//cout << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
	log_ofs << "segMap2SNPmer_phase_defined_bool: " << segMap2SNPmer_phase_defined_bool << endl;
	log_ofs << "segMap2SNPmer_phase1_bool: " << segMap2SNPmer_phase1_bool << endl;
	log_ofs << "segMap2SNPmer_phase2_bool: " << segMap2SNPmer_phase2_bool << endl;
	log_ofs << "segMap2SNPmer_phase2_bool_learned: " << segMap2SNPmer_phase2_bool_learned << endl;
	log_ofs << "snpLearned_success_bool: " << snpLearned_success_bool << endl; 
	Index_Info* indexInfo_SNP = new Index_Info();
	char *chrom_SNP;
    unsigned int *sa_SNP; 
    BYTE *lcpCompress_SNP;
    unsigned int *childTab_SNP;
    BYTE *verifyChild_SNP;
	if(segMap2SNPmer_phase1_bool)
	{	
		log_ofs << "start to load provided SNPmer index" << endl;
		string SNPfilePath = optionInfo->SNPfilePath;
		string SNP_seq_index_path = optionInfo->SNP_seq_index_path;
		//cout << "SNPfilePath: " << endl << SNPfilePath << endl;
		//log_ofs << "SNPfilePath: " << endl << SNPfilePath << endl;
		// Xinan: for now, do not replace ref bases with alternate bases. In the future, 
		// MPS3 should incorporate SNPs when doing sequence matching and canonical splice site detection
		indexInfo->insertSNP2chromStr(SNPfilePath, log_ofs);
		//cout << "start to load indexes" << endl;
		string indexStr_SNP = SNP_seq_index_path;
		//cout << "indexStr_SNP: " << endl << indexStr_SNP << endl;
		//log_ofs << "indexStr_SNP: " << endl << indexStr_SNP << endl;
		indexStr_SNP.append("/");
		string SA_file_SNP = indexStr_SNP; SA_file_SNP.append("_SA"); ifstream SA_file_ifs_SNP(SA_file_SNP.c_str(),ios::binary); 
		string lcpCompress_file_SNP = indexStr_SNP; lcpCompress_file_SNP.append("_lcpCompress"); ifstream lcpCompress_file_ifs_SNP(lcpCompress_file_SNP.c_str(),ios::binary);
		string childTab_file_SNP = indexStr_SNP; childTab_file_SNP.append("_childTab"); ifstream childTab_file_ifs_SNP(childTab_file_SNP.c_str(),ios::binary);
		string verifyChild_file_SNP = indexStr_SNP; verifyChild_file_SNP.append("_detChild"); ifstream verifyChild_file_ifs_SNP(verifyChild_file_SNP.c_str(),ios::binary);	
		string chrom_bit_file_SNP = indexStr_SNP; chrom_bit_file_SNP.append("_chrom"); ifstream chrom_bit_file_ifs_SNP(chrom_bit_file_SNP.c_str(),ios::binary);
		string parameter_file_SNP = indexStr_SNP; parameter_file_SNP.append("_parameter"); ifstream parameter_file_ifs_SNP(parameter_file_SNP.c_str(),ios::binary);

		indexInfo_SNP->initiate(parameter_file_ifs_SNP, log_ofs);
		chrom_SNP = (char*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs_SNP.read((char*)chrom_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
		indexInfo_SNP->readGenome(chrom_SNP);
		indexInfo_SNP->initiate();	
		indexInfo_SNP->initiateChrNameIndexArray(1000);
	    sa_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs_SNP.read((char*)sa_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
		lcpCompress_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs_SNP.read((char*)lcpCompress_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));	 
		childTab_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs_SNP.read((char*)childTab_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
		verifyChild_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs_SNP.read((char*)verifyChild_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
		//cout << "SyntheticSNPtransSeq index files loaded" << endl;
		log_ofs << "end of loading provided SNPmer index" << endl;
	}
	else if(updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool && updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool)
	{
		cout << "updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool: bothPhase" << endl;
		log_ofs << "updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool: bothPhase" << endl;
		string SNPfilePath = optionInfo->SNPfilePath;
		indexInfo->insertSNP2chromStr(SNPfilePath, log_ofs);
	}
	else
	{}
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       Personal Genome           ////////////////////////////////////////////////////////////////
	Stats_Info* statsInfo = new Stats_Info();
	if(SE_or_PE_bool)
		statsInfo->initiate_stats_info_SE(threads_num);
	else
		statsInfo->initiate_stats_info_PE(threads_num);
	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////       1st Mapping Process           ////////////////////////////////////////////////////////////////
	/* align main*/
   	string tmpHeadSectionInfo = outputDirStr + "/headSectionInfo";
   	ofstream tmpHeadSectionInfo_ofs(tmpHeadSectionInfo.c_str());    	
   	tmpHeadSectionInfo_ofs << headerSectionInfo->returnHeaderSectionInfoStr() << endl;
	tmpHeadSectionInfo_ofs.close();
	string repeatRegionFile = outputDirStr + "/phase1_output/repeat_region/repeatRegion";
	ofstream repeatRegionFile_ofs(repeatRegionFile.c_str());
	string tmpAlignCompleteRead_SE = outputDirStr + "/phase1_output/completePair/complete_SE.sam";
	ofstream tmpAlignCompleteRead_SE_ofs(tmpAlignCompleteRead_SE.c_str());
	string tmpAlignCompleteRead = outputDirStr + "/phase1_output/completePair/completePair.sam";
	ofstream tmpAlignCompleteRead_ofs(tmpAlignCompleteRead.c_str());
	string tmpAlignOneEndUnmapped = outputDirStr + "/phase1_output/oneEndUnmapped/oneEndUnmapped";
	if(Do_Phase1_Only)
		tmpAlignOneEndUnmapped += ".sam";	
	else
		tmpAlignOneEndUnmapped += ".alignInfo";
	ofstream tmpAlignOneEndUnmapped_ofs(tmpAlignOneEndUnmapped.c_str());
	string tmpAlignUnmapped_mappedToRepeatRegionFile_SE = outputDirStr + "/phase1_output/repeat_region/unmapped_mappedToRepeatRegion_SE.sam";
	ofstream tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs(tmpAlignUnmapped_mappedToRepeatRegionFile_SE.c_str());
	string tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = outputDirStr + "/phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam";
	ofstream tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs(tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
	string tmpAlignUnmapped_SE = outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_SE.sam";
	ofstream tmpAlignUnmapped_SE_ofs(tmpAlignUnmapped_SE.c_str());
	string tmpAlignBothEndsUnmapped = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";	
	ofstream tmpAlignBothEndsUnmapped_ofs(tmpAlignBothEndsUnmapped.c_str());
	string tmpAlignUnmapped_lowScore_SE = outputDirStr + "/phase1_output/bothEndsUnmapped/unmapped_lowScore_SE.sam";
	ofstream tmpAlignUnmapped_lowScore_SE_ofs(tmpAlignUnmapped_lowScore_SE.c_str());
	string tmpAlignBothEndsUnmapped_lowScore = outputDirStr + "/phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
	ofstream tmpAlignBothEndsUnmapped_lowScore_ofs(tmpAlignBothEndsUnmapped_lowScore.c_str());
	string tmpAlignIncomplete_SE = outputDirStr + "/phase1_output/incomplete/incomplete_SE.alignInfo";
	ofstream tmpAlignIncomplete_SE_ofs(tmpAlignIncomplete_SE);
	string tmpAlignIncompletePair = outputDirStr + "/phase1_output/incomplete/incomplete.alignInfo";
	ofstream tmpAlignIncompletePair_ofs(tmpAlignIncompletePair.c_str());
	string tmpAlignIncomplete_SE_SAM = outputDirStr + "/phase1_output/incomplete/incomplete_SE.sam";
	ofstream tmpAlignIncomplete_SE_SAM_ofs(tmpAlignIncomplete_SE_SAM.c_str());
	string tmpAlignIncompletePair_SAM = outputDirStr + "/phase1_output/incomplete/incompletePair.sam"; 
	ofstream tmpAlignIncompletePair_SAM_ofs(tmpAlignIncompletePair_SAM.c_str());	
	string tmpIntermediateJunctionFile = outputDirStr + "/phase2_output/inter.junc";

	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
	if(SE_or_PE_bool)
		inputRead_PE_ifs.close();

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE, line2_afterProcess, line2_PE_afterProcess;
	vector< RepeatRegion_Info* > repeatRegionInfoVec;
	for(int tmp = 0; tmp < threads_num; tmp++)
	{
		RepeatRegion_Info* repeatRegionInfo = new RepeatRegion_Info();
		repeatRegionInfoVec.push_back(repeatRegionInfo);
	}
	loadGlobalIndex_end = clock();
	phase1map_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... 1st mapping process starts ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... 1st mapping process starts ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... 1st mapping process starts ......" << endl; 

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();
	Read_Array_Queue* readArrayQueue = new Read_Array_Queue();
	Result_Array_Queue* resultArrayQueue = new Result_Array_Queue();
	int tmpInputReadNumInBatchArray_phase1 = normalRecordNum_1stMapping;//ReadNumInReadArray_Phase1;
	int tmpInputTimeWeight_phase1 = InputTimeWeight_Phase1;
	int tmpOutputTimeWeigth_phase1 = OutputTimeWeight_Phase1; 
	bool endOfFile_bool = false;
	bool endOfProcessing_bool = false;

	if(separateThreadForIO_bool)
	{	
		omp_set_num_threads(2);
		omp_set_nested(1);
	#pragma omp parallel
		{
	#pragma omp sections
			{
	#pragma omp section
				io_stage_phase1_separateThreadForIO(inputRead_ifs, inputRead_PE_ifs, readArrayQueue, resultArrayQueue, endOfFile_bool, endOfProcessing_bool, 
					tmpInputReadNumInBatchArray_phase1, tmpInputTimeWeight_phase1, tmpOutputTimeWeigth_phase1, log_ofs, readPreProcessInfo,
					tmpAlignCompleteRead_ofs, tmpAlignIncompletePair_ofs, tmpAlignOneEndUnmapped_ofs, tmpAlignBothEndsUnmapped_ofs,
					tmpAlignBothEndsUnmapped_lowScore_ofs, tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
					tmpAlignIncompletePair_SAM_ofs, input_log_ofs_phase1, output_log_ofs_phase1, fasta_or_fastq_bool, SE_or_PE_bool, readTotalNum);
	#pragma omp section
				process_stage_phase1_separateThreadForIO(readArrayQueue, resultArrayQueue, endOfFile_bool, endOfProcessing_bool, threads_num-1, 
					sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, preIndexIntervalStartArray,
					preIndexIntervalEndArray, repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
					Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only, Do_Phase1_Only,	statsInfo, fasta_or_fastq_bool, 
					mapping_log_ofs_phase1, checkQualSeqForReadSegSeq, SE_or_PE_bool, sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
					indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase1_bool, mismatchHashInfoVec, segMap2SNPmer_phase2_bool_learned);
			}
		}
	}
	else
	{
		io_process_phase1_allThreadsSharedByBothStage(inputRead_ifs, inputRead_PE_ifs, readArrayQueue, resultArrayQueue,
			tmpInputReadNumInBatchArray_phase1, log_ofs, readPreProcessInfo, tmpAlignCompleteRead_ofs, 
			tmpAlignIncompletePair_ofs, tmpAlignOneEndUnmapped_ofs, tmpAlignBothEndsUnmapped_ofs, tmpAlignBothEndsUnmapped_lowScore_ofs, 
			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs, tmpAlignIncompletePair_SAM_ofs, input_log_ofs_phase1, output_log_ofs_phase1, 
			fasta_or_fastq_bool, SE_or_PE_bool, threads_num, sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, preIndexMapLengthArray, 
			preIndexIntervalStartArray, preIndexIntervalEndArray, repeatRegionInfoVec, Do_cirRNA, Do_extendHeadTail_phase1, annotation_provided_bool, 
			Do_annotation_only_bool, annotationInfo, outputDirectlyBool_Phase1Only, Do_Phase1_Only,	statsInfo, mapping_log_ofs_phase1, 
			checkQualSeqForReadSegSeq, sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, 
			SNPlocInSyntheticSNPseq, segMap2SNPmer_phase1_bool, mismatchHashInfoVec, segMap2SNPmer_phase2_bool_learned, readTotalNum);
	}
	//log_ofs << "perfectMatch_pair #: " << perfectMatch_pair << endl;
	repeatRegionFile_ofs << "Repeat Region Info: size = " << repeatRegionInfoVec.size() << endl;
	for(int tmpThread = 0; tmpThread < threads_num; tmpThread++)
		repeatRegionInfoVec[tmpThread]->outputRepeatRegion(tmpThread+1, indexInfo, sa, 100, repeatRegionFile_ofs);

	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	tmpAlignCompleteRead_ofs.close();
	tmpAlignOneEndUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_ofs.close();
	tmpAlignBothEndsUnmapped_lowScore_ofs.close();
	tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs.close();
	tmpAlignIncompletePair_SAM_ofs.close();
	if(Do_Phase1_Only)
		tmpAlignIncompletePair_ofs.close();
	tmpAlignCompleteRead_SE_ofs.close();
	tmpAlignUnmapped_SE_ofs.close();
	tmpAlignUnmapped_lowScore_SE_ofs.close();
	tmpAlignUnmapped_mappedToRepeatRegionFile_SE_ofs.close();
	tmpAlignIncomplete_SE_SAM_ofs.close();
	if(Do_Phase1_Only)
		tmpAlignIncomplete_SE_ofs.close();

	free(preIndexMapLengthArray); free(preIndexIntervalStartArray); free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);
	free(verifyChild);
	free(chrom);
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... 1st mapping process ends ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... 1st mapping process ends ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... 1st mapping process ends ......" << endl; 
	log_ofs << endl << "**********************************" << endl << "**********************************";
	runtime_log_ofs << endl << "**********************************" << endl << "**********************************";
	phase1map_end = clock();
	///////////////////////////////////////////      merge learned candi SNPs      ////////////////////////////////////////////
	///////////////////////////////////////////      merge learned candi SNPs      ////////////////////////////////////////////
	///////////////////////////////////////////      merge learned candi SNPs      ////////////////////////////////////////////
	buildSNPcontext_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start to build SNPmer index for learned SNPs from phase1 ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... start to build SNPmer index for learned SNPs from phase1 ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... start to build SNPmer index for learned SNPs from phase1 ......" << endl; 

	LearnedCandiSNPhash_Info mismatchHashInfo_merged;
	mismatchHashInfo_merged.initiate(chromNum);	
	// merge 2 one
	mismatchHashInfoVec.mergeLearnedCandiSNPhashInfoVec2one(mismatchHashInfo_merged, indexInfo);
	mismatchHashInfoVec.freeMemory();
	string dir_SNPmer_learned = outputDirStr + "/SNPmer_learned";
   	string filteredLearnedSNP_file = dir_SNPmer_learned + "/filteredLearnedSNP.txt";
   	int called_snp_num = 0;
	int filterLearnedSNP_supNumMin = FILTER_LEARNSNP_SUPNUMMIN;
	double filterLearnedSNP_ratioMin = FILTER_LEARNSNP_RATIOMIN;
	mismatchHashInfo_merged.output_filteredLearnedSNP(filteredLearnedSNP_file, called_snp_num, 
		filterLearnedSNP_supNumMin, filterLearnedSNP_ratioMin, indexInfo);
	//cout << "filterLearnedSNP_supNumMin: " << filterLearnedSNP_supNumMin << endl;
	//cout << "filterLearnedSNP_ratioMin: " << filterLearnedSNP_ratioMin << endl;
	log_ofs << "filterLearnedSNP_supNumMin: " << filterLearnedSNP_supNumMin << endl;
	log_ofs << "filterLearnedSNP_ratioMin: " << filterLearnedSNP_ratioMin << endl;	
	//cout << "called_snp_num: " << called_snp_num << endl;
	log_ofs << "called_snp_num: " << called_snp_num << endl;
	if(called_snp_num == 0)
		snpLearned_success_bool = false;
   	// pull out SNPmers ... 
	string filteredLearnedSNP_SNPmer_fa_file = dir_SNPmer_learned + "/filteredLearnedSNP_SNPmer.fa";
	if(snpLearned_success_bool)
		indexInfo->insertSNP2chromStr_outputSNPmer(filteredLearnedSNP_file, filteredLearnedSNP_SNPmer_fa_file, log_ofs, SNPmerLength);
	// indexing SNPmers
   	string dir_SNPmer_index_learned = dir_SNPmer_learned + "/index";

   	// separate SNPmer index building program
   	if(snpLearned_success_bool)
   	{	
   		string cmd_buildSNPmerIndex = "buildSNPmerIndex_embeded_mps_nonBatch " + filteredLearnedSNP_SNPmer_fa_file + " " + dir_SNPmer_index_learned;
   		const int buildSNPmerIndex_err = system(cmd_buildSNPmerIndex.c_str());
   		if(-1 == buildSNPmerIndex_err)
   		{
	   		log_ofs << "error in buildSNPmerIndex ..." << endl;
	   		log_ofs << "buildSNPmerIndex_err: " << cmd_buildSNPmerIndex << endl;
	   		cout << "error in buildSNPmerIndex ..." << endl;
	   		cout << "buildSNPmerIndex_err: " << cmd_buildSNPmerIndex << endl;	   		
	   		//exit(1);   			
   		}
   		string dir_SNPmer_index_learned_success_file = dir_SNPmer_index_learned + "/results.txt";
   		string buildSNPmerIndex_result_str;
   		ifstream dir_SNPmer_index_learned_success_ifs(dir_SNPmer_index_learned_success_file.c_str());
   		getline(dir_SNPmer_index_learned_success_ifs, buildSNPmerIndex_result_str);
   		dir_SNPmer_index_learned_success_ifs.close();
   		if(buildSNPmerIndex_result_str == "success")
   			snpLearned_success_bool = true;
   		else //if(buildSNPmerIndex_result_str == "fail")
   			snpLearned_success_bool = false;
   		cout << "buildSNPmerIndex_result_str: " << buildSNPmerIndex_result_str << endl;
   		cout << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
   	}
 	//Build_snpMerIndex buildSnpMerIndexInfo;
 	//if(snpLearned_success_bool)
	//	snpLearned_success_bool = buildSnpMerIndexInfo.build_snpMerIndex(filteredLearnedSNP_SNPmer_fa_file, dir_SNPmer_index_learned); 
	//cout << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
	log_ofs << "snpLearned_success_bool: " << snpLearned_success_bool << endl;
	// loading SNPmers
	Index_Info* indexInfo_SNP_learned = new Index_Info();
	char *chrom_SNP_learned;
    unsigned int *sa_SNP_learned; 
    BYTE *lcpCompress_SNP_learned;
    unsigned int *childTab_SNP_learned;
    BYTE *verifyChild_SNP_learned;

    if(snpLearned_success_bool)
    {
	    string indexStr_SNP_learned = dir_SNPmer_index_learned + "/";
		string SA_file_SNP_learned = indexStr_SNP_learned + "_SA"; 
		ifstream SA_file_ifs_SNP_learned(SA_file_SNP_learned.c_str(),ios::binary); 
		string lcpCompress_file_SNP_learned = indexStr_SNP_learned + "_lcpCompress"; 
		ifstream lcpCompress_file_ifs_SNP_learned(lcpCompress_file_SNP_learned.c_str(),ios::binary);
		string childTab_file_SNP_learned = indexStr_SNP_learned + "_childTab"; 
		ifstream childTab_file_ifs_SNP_learned(childTab_file_SNP_learned.c_str(),ios::binary);
		string verifyChild_file_SNP_learned = indexStr_SNP_learned + "_detChild"; 
		ifstream verifyChild_file_ifs_SNP_learned(verifyChild_file_SNP_learned.c_str(),ios::binary);	
		string chrom_bit_file_SNP_learned = indexStr_SNP_learned + "_chrom"; 
		ifstream chrom_bit_file_ifs_SNP_learned(chrom_bit_file_SNP_learned.c_str(),ios::binary);
		string parameter_file_SNP_learned = indexStr_SNP_learned + "_parameter"; 
		ifstream parameter_file_ifs_SNP_learned(parameter_file_SNP_learned.c_str(),ios::binary);

		indexInfo_SNP_learned->initiate(parameter_file_ifs_SNP_learned, log_ofs);
		chrom_SNP_learned = (char*)malloc((indexInfo_SNP_learned->returnIndexSize()) * sizeof(char)); 
		chrom_bit_file_ifs_SNP_learned.read((char*)chrom_SNP_learned, (indexInfo_SNP_learned->returnIndexSize()) * sizeof(char)); 
		indexInfo_SNP_learned->readGenome(chrom_SNP_learned);
		indexInfo_SNP_learned->initiate();	
		indexInfo_SNP_learned->initiateChrNameIndexArray(1000);
		sa_SNP_learned = (unsigned int*)malloc((indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int)); 
		SA_file_ifs_SNP_learned.read((char*)sa_SNP_learned, (indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int));	 
		lcpCompress_SNP_learned = (BYTE*)malloc((indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE)); 
		lcpCompress_file_ifs_SNP_learned.read((char*)lcpCompress_SNP_learned, (indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE));	 
		childTab_SNP_learned = (unsigned int*)malloc((indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int)); 
		childTab_file_ifs_SNP_learned.read((char*)childTab_SNP_learned, (indexInfo_SNP_learned->returnIndexSize()) * sizeof(unsigned int));	 
		verifyChild_SNP_learned = (BYTE*)malloc((indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE)); 
		verifyChild_file_ifs_SNP_learned.read((char*)verifyChild_SNP_learned, (indexInfo_SNP_learned->returnIndexSize()) * sizeof(BYTE));

		SA_file_ifs_SNP_learned.close();
		lcpCompress_file_ifs_SNP_learned.close();
		childTab_file_ifs_SNP_learned.close();
		verifyChild_file_ifs_SNP_learned.close();
		chrom_bit_file_ifs_SNP_learned.close();
		parameter_file_ifs_SNP_learned.close();	
		cout << "Self learned SNPmer index files loaded" << endl;
		log_ofs << "Self learned SNPmer index files loaded" << endl;
	}

	///////////////////////////////////////////   build SNPmer index for phase2    ////////////////////////////////////////////
	if(segMap2SNPmer_phase2_bool && (!segMap2SNPmer_phase1_bool))
	{	
		string SNPfilePath = optionInfo->SNPfilePath;
		string SNP_seq_index_path = optionInfo->SNP_seq_index_path;		
		cout << "SNPfilePath: " << endl << SNPfilePath << endl;
		log_ofs << "SNPfilePath: " << endl << SNPfilePath << endl;
		// Xinan: for now, do not replace ref bases with alternate bases. In the future, 
		// MPS3 should incorporate SNPs when doing sequence matching and canonical splice site detection
		indexInfo->insertSNP2chromStr(SNPfilePath, log_ofs);
		cout << "start to load indexes" << endl;
		string indexStr_SNP = SNP_seq_index_path;
		cout << "indexStr_SNP: " << endl << indexStr_SNP << endl;
		log_ofs << "indexStr_SNP: " << endl << indexStr_SNP << endl;
		indexStr_SNP.append("/");
		string SA_file_SNP = indexStr_SNP; SA_file_SNP.append("_SA"); ifstream SA_file_ifs_SNP(SA_file_SNP.c_str(),ios::binary); 
		string lcpCompress_file_SNP = indexStr_SNP; lcpCompress_file_SNP.append("_lcpCompress"); ifstream lcpCompress_file_ifs_SNP(lcpCompress_file_SNP.c_str(),ios::binary);
		string childTab_file_SNP = indexStr_SNP; childTab_file_SNP.append("_childTab"); ifstream childTab_file_ifs_SNP(childTab_file_SNP.c_str(),ios::binary);
		string verifyChild_file_SNP = indexStr_SNP; verifyChild_file_SNP.append("_detChild"); ifstream verifyChild_file_ifs_SNP(verifyChild_file_SNP.c_str(),ios::binary);	
		string chrom_bit_file_SNP = indexStr_SNP; chrom_bit_file_SNP.append("_chrom"); ifstream chrom_bit_file_ifs_SNP(chrom_bit_file_SNP.c_str(),ios::binary);
		string parameter_file_SNP = indexStr_SNP; parameter_file_SNP.append("_parameter"); ifstream parameter_file_ifs_SNP(parameter_file_SNP.c_str(),ios::binary);

		indexInfo_SNP->initiate(parameter_file_ifs_SNP, log_ofs);
		chrom_SNP = (char*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs_SNP.read((char*)chrom_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(char)); 
		indexInfo_SNP->readGenome(chrom_SNP);
		indexInfo_SNP->initiate();	
		indexInfo_SNP->initiateChrNameIndexArray(1000);
	    sa_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); SA_file_ifs_SNP.read((char*)sa_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
		lcpCompress_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); lcpCompress_file_ifs_SNP.read((char*)lcpCompress_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));	 
		childTab_SNP = (unsigned int*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int)); childTab_file_ifs_SNP.read((char*)childTab_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(unsigned int));	 
		verifyChild_SNP = (BYTE*)malloc((indexInfo_SNP->returnIndexSize()) * sizeof(BYTE)); verifyChild_file_ifs_SNP.read((char*)verifyChild_SNP, (indexInfo_SNP->returnIndexSize()) * sizeof(BYTE));
		cout << "SyntheticSNPtransSeq index files loaded" << endl;
	}
	else if(updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool && (!updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool))
	{
		cout << "updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool: phase2only" << endl;
		string SNPfilePath = optionInfo->SNPfilePath;
		cout << "SNPfilePath: " << endl << SNPfilePath << endl;
		indexInfo->insertSNP2chromStr(SNPfilePath, log_ofs);
	}
	else
	{

	}
	buildSNPcontext_end = clock();
	loadLocalIndex_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... end of building SNPmer index for learned SNPs from phase1 ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... end of building SNPmer index for learned SNPs from phase1 ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... end of building SNPmer index for learned SNPs from phase1 ......" << endl; 
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	cout << endl << tmpTimeStr << "... start to load 2nd level index ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... load 2nd level index starts ......" << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;
	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	if(load2ndLevelIndexBool)
	{
		//log_ofs << "start to load second-level index ..." << endl;
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
        {  
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
          log_ofs << "invalid index: " << endl << inputIndexFileStr << endl;
        }
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
  			{
        	indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				  log_ofs << "invalid index: " << endl << inputIndexFileStr << endl;
        }
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
			//log_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
		}
		//log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
		//loadIndex_end = clock(); 
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... load 2nd level index ends ......" << endl;		
	log_ofs << endl << tmpTimeStr << "... load 2nd level index ends ......" << endl;
	loadLocalIndex_end = clock();
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On one end unmapped Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	oneEndMap_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts ......" << endl;	 
	log_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts ......" << endl;
	runtime_log_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads starts ......" << endl; 

	string OutputSamFile_oneEndMapped = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam";
	ofstream OutputSamFile_oneEndMapped_ofs(OutputSamFile_oneEndMapped.c_str());
	string OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = outputDirStr + "/phase2_output/oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
	ofstream OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore.c_str());
	string OutputSamFile_oneEndMapped_unpair = outputDirStr + "/phase2_output/oneEndUnmapped.unpaired.sam";
	ofstream OutputSamFile_oneEndMapped_unpair_ofs(OutputSamFile_oneEndMapped_unpair.c_str());
	string OutputSamFile_oneEndMapped_alignInfo = outputDirStr + "/phase2_output/oneEndUnmapped.pairedComplete.sam_alignInfo";
	ofstream OutputSamFile_oneEndMapped_alignInfo_ofs(OutputSamFile_oneEndMapped_alignInfo.c_str());	

	int tmpRecordNum_oneEndUnmapped = 0;
	if(SE_or_PE_bool)
		DoRemappingOnUnmapEndReadsBool = false;

	if(DoRemappingOnUnmapEndReadsBool)
	{
		// nowtime = time(NULL);
		// local = localtime(&nowtime);
		// tmpTimeStr = asctime(local);
		// tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		// log_ofs << endl << tmpTimeStr << "... start doing remapping on unmapped end reads" << endl;
		// runtime_log_ofs << endl << tmpTimeStr << "... start doing remapping on unmapped end reads" << endl;
		// cout << endl << tmpTimeStr << "... start doing remapping on unmapped end reads" << endl;
		AlignInfoInput_Array_Queue* alignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixOneEndUnmapped_Array_Queue* fixOneEndUnmappedResultQueue = new Result_FixOneEndUnmapped_Array_Queue();

		string oneEndMappedFileStr = tmpAlignOneEndUnmapped;
		ifstream inputRecord_ifs(oneEndMappedFileStr.c_str());		
		bool endOfFile_bool = false;
		bool endOfProcessing_bool = false;
		int tmpInputReadNumInBatchArray_fixOneEndUnmapped = normalRecordNum_fixOneEndUnmapped;//ReadNumInReadArray_FixOneEndUnmapped;
		int tmpInputTimeWeight_fixOneEndUnmapped = InputTimeWeight_FixOneEndUnmapped;
		int tmpOutputTimeWeight_fixOneEndUnmapped = OutputTimeWeight_FixOneEndUnmapped;

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
							inputRecord_ifs, alignInfoInputQueue, fixOneEndUnmappedResultQueue, endOfFile_bool, endOfProcessing_bool,
							tmpInputReadNumInBatchArray_fixOneEndUnmapped, tmpInputTimeWeight_fixOneEndUnmapped,
							tmpOutputTimeWeight_fixOneEndUnmapped, log_ofs, OutputSamFile_oneEndMapped_ofs,
							tmpAlignIncompletePair_ofs, OutputSamFile_oneEndMapped_unpair_ofs, 
							OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs,
							input_log_ofs_phase2_fixOneEndUnmapped, output_log_ofs_phase2_fixOneEndUnmapped);
		#pragma omp section
						process_stage_fixOneEndUnmapped_separateThreadForIO(
							alignInfoInputQueue, fixOneEndUnmappedResultQueue, endOfFile_bool, endOfProcessing_bool, threads_num-1,
							fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab,
							secondLevelDetChild, indexInfo, Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, 
							Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq,
							mapping_log_ofs_phase2_fixOneEndUnmapped, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
							// provided SNP
							sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
							indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool, 
							// learned SNP
							sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
							indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, snpLearned_success_bool);
					}
				}
		}
		else
		{
			io_process_fixOneEndUnmapped_allThreadsSharedByBothStage(
				inputRecord_ifs, alignInfoInputQueue, fixOneEndUnmappedResultQueue,	tmpInputReadNumInBatchArray_fixOneEndUnmapped,
				log_ofs, OutputSamFile_oneEndMapped_ofs, tmpAlignIncompletePair_ofs, OutputSamFile_oneEndMapped_unpair_ofs,
				OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs, input_log_ofs_phase2_fixOneEndUnmapped,
				output_log_ofs_phase2_fixOneEndUnmapped, threads_num, fasta_or_fastq_bool, statsInfo, secondLevelChrom,
				secondLevelSa, secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, 
				Do_extendHeadTail_fixOneEndUnmapped, annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
				MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, mapping_log_ofs_phase2_fixOneEndUnmapped, SE_or_PE_bool,
				outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
				// provided SNP
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
				indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool,
				// learned SNP
				sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
				indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, snpLearned_success_bool);
		}
		alignInfoInputQueue->free();
		delete alignInfoInputQueue;
		fixOneEndUnmappedResultQueue->free();
		delete fixOneEndUnmappedResultQueue;
		inputRecord_ifs.close();
	}

	OutputSamFile_oneEndMapped_ofs.close();
	OutputSamFile_oneEndMapped_unpair_ofs.close();
	tmpAlignIncompletePair_ofs.close();
	OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs.close();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... fixing oneEndUnmapped reads ends ......" << endl;  
	log_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads ends ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... fixing oneEndUnmapped reads ends ......" << endl; 
	log_ofs << endl << "**********************************************************************************" << endl;
	runtime_log_ofs << endl << "**********************************************************************************" << endl;
	oneEndMap_end = clock();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Sam 2 Junc   ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	buildSpliceContext_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... generating SJhashInfo starts ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... generating SJhashInfo starts ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... generating SJhashInfo starts ......" << endl; 

	string juncfile = tmpIntermediateJunctionFile;
	string juncfile_alignInferHash = juncfile + ".alignInferHash";
	if(stopLearningSplicing_bool)
		DoSam2JuncBool = false;
	if(DoSam2JuncBool)
	{
		// generate SJ from already mapped reads
		vector<string> tmpAlignmentFileVec;
		if(SE_or_PE_bool)
		{
			tmpAlignmentFileVec.push_back(tmpAlignCompleteRead_SE);
			tmpAlignmentFileVec.push_back(tmpAlignIncomplete_SE_SAM);
		}
		else
		{
			tmpAlignmentFileVec.push_back(tmpAlignCompleteRead);
			tmpAlignmentFileVec.push_back(tmpAlignIncompletePair_SAM);
			tmpAlignmentFileVec.push_back(OutputSamFile_oneEndMapped);
		}

		AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec 
			= new AlignInferJunctionHash_Info_Vec();
		int alignInferJuncHashInfoVecSize = threads_num;
		alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(
			alignInferJuncHashInfoVecSize, chromNum);		
		alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(
			tmpAlignmentFileVec, indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(
			alignInferJunctionHashInfo, indexInfo);
		alignInferJunctionHashInfo->outputAlignInferInfoHashInfo_chrNamePos_supportNum(
			indexInfo, juncfile_alignInferHash);
		alignInferJunctionHashInfoVec->freeMemory();
		delete alignInferJunctionHashInfoVec;
	}

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... generating SJhashInfo ends ......" << endl; 
	log_ofs << endl << tmpTimeStr << "... generating SJhashInfo ends ......" << endl; 
	runtime_log_ofs << endl << tmpTimeStr << "... generating SJhashInfo ends ......" << endl; 
	buildSpliceContext_end = clock();
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	headTailMap_start = clock();
	string OutputSamFile_fixHeadTail_complete_pair = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_pair_ofs(OutputSamFile_fixHeadTail_complete_pair.c_str());
	string OutputSamFile_fixHeadTail_incomplete_pair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());
	string OutputSamFile_fixHeadTail_complete_unpair = outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());
	string OutputSamFile_fixHeadTail_incomplete_unpair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	
	string OutputSamFile_fixHeadTail_pair_lowScore = outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
	ofstream OutputSamFile_fixHeadTail_pair_lowScore_ofs(OutputSamFile_fixHeadTail_pair_lowScore.c_str());	
	string OutputSamFile_fixHeadTail_complete_SE = outputDirStr + "/phase2_output/fixHeadTail_complete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_complete_SE_ofs(OutputSamFile_fixHeadTail_complete_SE.c_str());
	string OutputSamFile_fixHeadTail_incomplete_SE = outputDirStr + "/phase2_output/fixHeadTail_incomplete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_SE_ofs(OutputSamFile_fixHeadTail_incomplete_SE.c_str());
	string OutputSamFile_fixHeadTail_lowScore_SE = outputDirStr + "/phase2_output/fixHeadTail_lowScore_SE.sam";
	ofstream OutputSamFile_fixHeadTail_lowScore_SE_ofs(OutputSamFile_fixHeadTail_lowScore_SE.c_str());

	// start to read splice junction
	int junctionNum = 0;
	int junctionNum_in_alignInferJuncHashInfo = 0;
	int junctionNum_in_annotation = 0;

	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	string InputSpliceJunction = juncfile_alignInferHash;
	// settings_log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	// log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	// cout << "InputSpliceJunction: " << InputSpliceJunction << endl;
	///////////////////////////////////////////////////////////////////
	if(DoRemappingOnUnfixedHeadTailAlignmentBool)
	{
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... start to build spliceJunction Hash" << endl;
		log_ofs << endl << tmpTimeStr << "... start to build spliceJunction Hash" << endl;
		runtime_log_ofs << endl << tmpTimeStr << "... start to build spliceJunction Hash" << endl;

    	bool spliceJunctionHashExists = true;
		string entryString;
		int tabLocation1, tabLocation2, tabLocation3, tabLocation4, tabLocation5;
		//char entry[500];
		int chrInt;
		int spliceStartPos;
		int spliceEndPos;
		string chrIntString;
		string spliceStartPosString;
		string spliceEndPosString;
		//////////////////////  string hash /////////////////////////////////////////
		if(!Do_annotation_only_bool)
		{	
	    	log_ofs << "start to load SJs in alignments" << endl;
    		//cout << "start to load SJs in alignments" << endl;
			if(!stopLearningSplicing_bool)
			{	
				alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
				junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
			}
			else
				junctionNum_in_alignInferJuncHashInfo = 0;
		}
		if(annotation_provided_bool)
		{	
			// loading SJs in annotation file (if provided)
	    	log_ofs << "start to load SJs in annotation" << endl;
    		//cout << "start to load SJs in annotation" << endl;
			ifstream annotatedSJ_ifs(annotation_file_path.c_str());
			while(!annotatedSJ_ifs.eof())
			{
				getline(annotatedSJ_ifs, entryString);
				if(entryString == "")
					break;
				junctionNum_in_annotation ++;
				//entryString = entry;
				tabLocation1 = entryString.find('\t', 0);
				tabLocation2 = entryString.find('\t', tabLocation1+1);
				tabLocation3 = entryString.find('\t', tabLocation2+1);
				chrIntString = entryString.substr(0, tabLocation1);
				spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
				if(tabLocation3 == string::npos)
					spliceEndPosString = entryString.substr(tabLocation2+1);
				else
					spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
				chrInt = indexInfo->convertStringToInt(chrIntString);
				if(chrInt >= 0)
				{	
					spliceStartPos = atoi(spliceStartPosString.c_str());
					spliceEndPos = atoi(spliceEndPosString.c_str());	
					SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
				}
			}
			annotatedSJ_ifs.close();
		}
		junctionNum = junctionNum_in_alignInferJuncHashInfo + junctionNum_in_annotation;
		if((junctionNum == 0))
			spliceJunctionHashExists = false;

		settings_log_ofs << "finish building spliceJunction Hash" << endl;
		settings_log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum << endl;
		settings_log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		
		log_ofs << "After inserting SJs generated from alignments and annotation, junctionNum = " << junctionNum << endl;
		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... end of building spliceJunction Hash" << endl;
		log_ofs << endl << tmpTimeStr << "... end of building spliceJunction Hash" << endl;
		runtime_log_ofs << endl << tmpTimeStr << "... end of building spliceJunction Hash" << endl;

		nowtime = time(NULL);
		local = localtime(&nowtime);
		tmpTimeStr = asctime(local);
		tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
		cout << endl << tmpTimeStr << "... fixing unfixed-head/tail reads starts ......" << endl ;  
		log_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads starts ......" << endl ; 
		runtime_log_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads starts ......" << endl; 

		string headTailSoftClippingFile;
		if(SE_or_PE_bool)
			headTailSoftClippingFile = tmpAlignIncomplete_SE;
		else
			headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";
		AlignInfoInput_Array_Queue* fixHeadTailAlignInfoInputQueue = new AlignInfoInput_Array_Queue();
		Result_FixHeadTail_Array_Queue* fixHeadTailResultQueue = new Result_FixHeadTail_Array_Queue();
		if(outputDirectlyBool_Phase1Only)
			headTailSoftClippingFile += ".all";
		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());
		//int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool endOfFile_bool = false;
		bool endOfProcessing_bool = false;

		int tmpInputReadNumInBatchArray_fixHeadTail = normalRecordNum_fixHeadTail;//ReadNumInReadArray_FixHeadTail;
		int tmpInputTimeWeight_fixHeadTail = InputTimeWeight_FixHeadTail;
		int tmpOutputTimeWeight_fixHeadTail = OutputTimeWeight_FixHeadTail;
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
						inputUnfixedHeadTailRecord_ifs, fixHeadTailAlignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool,
						endOfProcessing_bool, tmpInputReadNumInBatchArray_fixHeadTail, tmpInputTimeWeight_fixHeadTail,
						tmpOutputTimeWeight_fixHeadTail, log_ofs, OutputSamFile_fixHeadTail_complete_pair_ofs,
						OutputSamFile_fixHeadTail_incomplete_pair_ofs, OutputSamFile_fixHeadTail_complete_unpair_ofs,
						OutputSamFile_fixHeadTail_incomplete_unpair_ofs, OutputSamFile_fixHeadTail_pair_lowScore_ofs,
						input_log_ofs_phase2_fixHeadTail, output_log_ofs_phase2_fixHeadTail);
				#pragma omp section
					process_stage_fixHeadTail_separateThreadForIO(
						fixHeadTailAlignInfoInputQueue, fixHeadTailResultQueue, endOfFile_bool, endOfProcessing_bool,
						threads_num-1, fasta_or_fastq_bool, statsInfo, secondLevelChrom, secondLevelSa,
						secondLevelLcpCompress, secondLevelChildTab, secondLevelDetChild, indexInfo, SJ, 
						Do_extendHeadTail_fixHeadTail, annotation_provided_bool, Do_annotation_only_bool, 
						annotationInfo, checkQualSeqForReadSegSeq, checkQualSeqForShortAnchorSeqToTargetMap,
						spliceJunctionHashExists, Do_fixHeadTail_remapping, Do_fixHeadTail_greedyMapping,
						Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain,
						mapping_log_ofs_phase2_fixHeadTail, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
						// provided SNPs
						sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
						indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool, 
						// learned SNPs
						sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
						indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, snpLearned_success_bool);
				}
			}
		}
		else
		{
			io_process_fixHeadTail_allThreadsSharedByBothStage(
				inputUnfixedHeadTailRecord_ifs, fixHeadTailAlignInfoInputQueue, fixHeadTailResultQueue, 
				tmpInputReadNumInBatchArray_fixHeadTail, log_ofs,
				OutputSamFile_fixHeadTail_complete_pair_ofs, OutputSamFile_fixHeadTail_incomplete_pair_ofs,
				OutputSamFile_fixHeadTail_complete_unpair_ofs, OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
				OutputSamFile_fixHeadTail_pair_lowScore_ofs, input_log_ofs_phase2_fixHeadTail,
				output_log_ofs_phase2_fixHeadTail, threads_num, fasta_or_fastq_bool, statsInfo,
				secondLevelChrom, secondLevelSa, secondLevelLcpCompress, secondLevelChildTab,
				secondLevelDetChild, indexInfo, SJ, Do_extendHeadTail_fixHeadTail, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, checkQualSeqForReadSegSeq, 
				checkQualSeqForShortAnchorSeqToTargetMap, spliceJunctionHashExists, Do_fixHeadTail_remapping, 
				Do_fixHeadTail_greedyMapping, Do_fixHeadTail_remappingAndTargetMapping, Do_fixHeadTail_remappingAgain,
				mapping_log_ofs_phase2_fixHeadTail, SE_or_PE_bool, outputUnpairedSAM_bool, outputUnpairedSAM_bothEndsUniqueMappedOnly_bool,
				// provided SNPs
				sa_SNP, lcpCompress_SNP, childTab_SNP, chrom_SNP, verifyChild_SNP, 
				indexInfo_SNP, SNPlocInSyntheticSNPseq, segMap2SNPmer_phase2_bool, 
				// learned SNPs
				sa_SNP_learned, lcpCompress_SNP_learned, childTab_SNP_learned, chrom_SNP_learned, verifyChild_SNP_learned, 
				indexInfo_SNP_learned, SNPlocInSyntheticSNPseq_learned, segMap2SNPmer_phase2_bool_learned, snpLearned_success_bool);
		}
		fixHeadTailAlignInfoInputQueue->free();
		delete fixHeadTailAlignInfoInputQueue;
		fixHeadTailResultQueue->free();
		delete fixHeadTailResultQueue;
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete(SJ);
	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();
	OutputSamFile_fixHeadTail_complete_SE_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_SE_ofs.close();
	OutputSamFile_fixHeadTail_lowScore_SE_ofs.close();
	//delete learnedCandiSNPhashInfo_merged;
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... fixing unfixed-head/tail reads ends ......" << endl;  
	log_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads ends ......" << endl;  
	runtime_log_ofs << endl << tmpTimeStr << "... fixing unfixed-head/tail reads ends ......" << endl;  
	headTailMap_end = clock();
	if(SE_or_PE_bool)
	{
		//statsInfo->getPhase1Stats_SE();
		statsInfo->getFixHeadTailStats_SE();
		statsInfo->outputAllStats_SE_fixHeadTail(stats_ofs, readTotalNum);
		statsInfo->outputFinalStats_SE(stats_ofs, readTotalNum);
	}	
	else
	{	
		statsInfo->getPhase1Stats();
		statsInfo->getFixUnpairedStats();
		statsInfo->getFixHeadTailStats();
		//statsInfo->outputAllStats(log_ofs, readTotalNum);
		statsInfo->outputAllStats(stats_ofs, Do_Phase1_Only, readTotalNum);
		statsInfo->outputFinalStats(stats_ofs, Do_Phase1_Only, readTotalNum);	
	}

	samReport_start = clock();
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... start to prepare for final output files ......" << endl ;  
	log_ofs << endl << tmpTimeStr << "... start to prepare for final output files ......" << endl ;  
	runtime_log_ofs << endl << tmpTimeStr << "... start to prepare for final output files ......" << endl ;  

	string finalOutputSam = outputDirStr + "/output.sam";
  int cat_err;
	if(Do_Phase1_Only)
	{
		string cat_cmd;
		if(SE_or_PE_bool)
		{
			cat_cmd = "cat " + tmpHeadSectionInfo + " " + tmpAlignCompleteRead_SE + " " + tmpAlignIncomplete_SE_SAM + " " 
				+ tmpAlignUnmapped_mappedToRepeatRegionFile_SE + " " + tmpAlignUnmapped_lowScore_SE + " " + tmpAlignUnmapped_SE + " > " + finalOutputSam;			
		}
		else	
		{
			cat_cmd = "cat " + tmpHeadSectionInfo + " " + tmpAlignCompleteRead + " " + tmpAlignIncompletePair_SAM  + " " + tmpAlignOneEndUnmapped
				+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + " " + tmpAlignBothEndsUnmapped_lowScore
				+ " " + tmpAlignBothEndsUnmapped + " > " + finalOutputSam;
		}
		cat_err = system(cat_cmd.c_str());
		if(-1 == cat_err)
		{
	   		cout << "error in catting all files to a final output.sam ..." << endl;
	   		cout << "cat_err: " << cat_err << endl;
	   		//exit(1);
		}		
	}
	else
	{
		//cout << "start to merge all alignment files into one...." << endl;
		log_ofs << "start to merge all alignment files into one...." << endl;
		string cat_cmd;
		if(!SE_or_PE_bool)
		{	
      // cout << "tmpHeadSectionInfo: " << tmpHeadSectionInfo << endl;
      // cout << "tmpAlignCompleteRead: " << tmpAlignCompleteRead << endl;
      // cout << "OutputSamFile_oneEndMapped: " << OutputSamFile_oneEndMapped << endl;
      // cout << "OutputSamFile_fixHeadTail_complete_pair: " << OutputSamFile_fixHeadTail_complete_pair << endl;
      // cout << "OutputSamFile_fixHeadTail_incomplete_pair: " << OutputSamFile_fixHeadTail_incomplete_pair << endl;
      // cout << "OutputSamFile_fixHeadTail_complete_unpair: " << OutputSamFile_fixHeadTail_complete_unpair << endl;
      // cout << "OutputSamFile_fixHeadTail_incomplete_unpair: " << OutputSamFile_fixHeadTail_incomplete_unpair << endl;
      // cout << "OutputSamFile_oneEndMapped_unpair: " << OutputSamFile_oneEndMapped_unpair << endl;
      // cout << "tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile: " << tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile << endl;
      // cout << "tmpAlignBothEndsUnmapped_lowScore: " << tmpAlignBothEndsUnmapped_lowScore << endl;
      // cout << "OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore: " << OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore << endl;
      // cout << "OutputSamFile_fixHeadTail_pair_lowScore: " << OutputSamFile_fixHeadTail_pair_lowScore << endl;
      // cout << "tmpAlignBothEndsUnmapped: " << tmpAlignBothEndsUnmapped << endl;
      // cout << "finalOutputSam: " << finalOutputSam << endl;
			cat_cmd = "cat " + tmpHeadSectionInfo + " " 
        + tmpAlignCompleteRead + " " 
        + OutputSamFile_oneEndMapped + " " 
        + OutputSamFile_fixHeadTail_complete_pair + " " 
        + OutputSamFile_fixHeadTail_incomplete_pair + " " 
        + OutputSamFile_fixHeadTail_complete_unpair + " " 
        + OutputSamFile_fixHeadTail_incomplete_unpair + " " 
        + OutputSamFile_oneEndMapped_unpair + " " 
        + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + " " 
        + tmpAlignBothEndsUnmapped_lowScore + " " 
        + OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore + " " 
        + OutputSamFile_fixHeadTail_pair_lowScore + " " 
        + tmpAlignBothEndsUnmapped + " > " + finalOutputSam;
		}
		else
		{
      // cout << "tmpHeadSectionInfo: " << tmpHeadSectionInfo << endl;
      // cout << "tmpAlignCompleteRead_SE: " << tmpAlignCompleteRead_SE << endl;
      // cout << "OutputSamFile_fixHeadTail_complete_SE: " << OutputSamFile_fixHeadTail_complete_SE << endl;
      // cout << "OutputSamFile_fixHeadTail_incomplete_SE: " << OutputSamFile_fixHeadTail_incomplete_SE << endl;
      // cout << "tmpAlignUnmapped_mappedToRepeatRegionFile_SE: " << tmpAlignUnmapped_mappedToRepeatRegionFile_SE << endl;
      // cout << "tmpAlignUnmapped_lowScore_SE: " << tmpAlignUnmapped_lowScore_SE << endl;
      // cout << "OutputSamFile_fixHeadTail_lowScore_SE: " << OutputSamFile_fixHeadTail_lowScore_SE << endl;
      // cout << "tmpAlignUnmapped_SE: " << tmpAlignUnmapped_SE << endl;
      // cout << "finalOutputSam: " << finalOutputSam << endl;
			cat_cmd = "cat " + tmpHeadSectionInfo + " " 
        + tmpAlignCompleteRead_SE + " " 
        + OutputSamFile_fixHeadTail_complete_SE + " " 
        + OutputSamFile_fixHeadTail_incomplete_SE + " " 
        + tmpAlignUnmapped_mappedToRepeatRegionFile_SE + " " 
        + tmpAlignUnmapped_lowScore_SE + " " 
        + OutputSamFile_fixHeadTail_lowScore_SE + " " 
        + tmpAlignUnmapped_SE + " > " + finalOutputSam;
		}
		cat_err = system(cat_cmd.c_str());
		if(-1 == cat_err)
		{
	   		log_ofs << "error in catting all files to a final output.sam ..." << endl;
	   		log_ofs << "cat_err: " << cat_err << endl;
	   		log_ofs << "re-reading all files and printing to a final output.sam" << endl;
        vector<string> fileVec;
        if(!SE_or_PE_bool) // PE
        {
          fileVec.push_back(tmpHeadSectionInfo);
          fileVec.push_back(tmpAlignCompleteRead);
          fileVec.push_back(OutputSamFile_oneEndMapped);
          fileVec.push_back(OutputSamFile_fixHeadTail_complete_pair);
          fileVec.push_back(OutputSamFile_fixHeadTail_incomplete_pair);
          fileVec.push_back(OutputSamFile_fixHeadTail_complete_unpair);
          fileVec.push_back(OutputSamFile_fixHeadTail_incomplete_unpair);
          fileVec.push_back(OutputSamFile_oneEndMapped_unpair);
          fileVec.push_back(tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile);
          fileVec.push_back(tmpAlignBothEndsUnmapped_lowScore);
          fileVec.push_back(OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore);
          fileVec.push_back(OutputSamFile_fixHeadTail_pair_lowScore);
          fileVec.push_back(tmpAlignBothEndsUnmapped);
        }
        else // SE
        {
          fileVec.push_back(tmpHeadSectionInfo);
          fileVec.push_back(tmpAlignCompleteRead_SE);
          fileVec.push_back(OutputSamFile_fixHeadTail_complete_SE);
          fileVec.push_back(OutputSamFile_fixHeadTail_incomplete_SE);
          fileVec.push_back(tmpAlignUnmapped_mappedToRepeatRegionFile_SE);
          fileVec.push_back(tmpAlignUnmapped_lowScore_SE);
          fileVec.push_back(OutputSamFile_fixHeadTail_lowScore_SE);
          fileVec.push_back(tmpAlignUnmapped_SE);
        }
        merge2finalOutput(fileVec, finalOutputSam);
		    log_ofs << "end of generating the final output.sam" << endl;
    }
	}

  annotation_ifs.close();
  delete alignInferJunctionHashInfo;
  delete annotationInfo;  
  delete indexInfo;

	samReport_end = clock();
  // if(cat_err == -1)
  // { 
  // 	juncReport_start = clock();
  // 	if(reportJunc_bool)
  // 	{
  // 		string sam2junc_folder = outputDirStr + "/sam2junc";
  // 		//string mkdir_sam2junc_cmd = "mkdir -p " + sam2junc_folder;
  // 		//system(mkdir_sam2junc_cmd.c_str());
  // 		string reportJunc_file = outputDirStr + "/output.junc";
  // 		string reportJunc_cmd = "reportJunc_embeded " + indexStr + " " + int_to_str(threads_num) + " "
  // 			+ sam2junc_folder + " " + finalOutputSam + " " + reportJunc_file;
  // 		const int reportJunc_err = system(reportJunc_cmd.c_str());
  //     if(-1 == reportJunc_err)
  //     {
  //       cout << "error in report junc ..." << endl;
  //       cout << "reportJunc_err: " << reportJunc_err << endl;
  //     }
  // 	}
  // 	juncReport_end = clock();
  // 	//string rm_headerSection = "rm " + outputDirStr + "/headSectionInfo";
  // 	//system(rm_headerSection.c_str());
  // 	if(removeAllIntermediateFilesBool)
  // 	{
  // 		string rm_logFolder_cmd = "mkdir -p " + outputDirStr + "/logs";
  // 		string rm_phase1output_cmd = "mkdir -p " + outputDirStr + "/phase1_output";
  // 		string rm_phase2output_cmd = "mkdir -p " + outputDirStr + "/phase2_output";
  // 		string rm_sam2junc_cmd = "mkdir -p " + outputDirStr + "/sam2junc";
  // 		string rm_SNPmerLearned_cmd = "mkdir -p " + outputDirStr + "/SNPmer_learned";
  // 		const int rm_logFolder_err = system(rm_logFolder_cmd.c_str());
  // 		const int rm_phase1output_err = system(rm_phase1output_cmd.c_str());
  // 		const int rm_phase2output_err = system(rm_phase2output_cmd.c_str());
  // 		const int rm_sam2junc_err = system(rm_sam2junc_cmd.c_str());
  // 		const int rm_SNPmerLearned_err = system(rm_SNPmerLearned_cmd.c_str());
  //     if((-1 == rm_logFolder_err)||(-1 == rm_phase1output_err)||(-1 == rm_phase2output_err)
  //       ||(-1 == rm_sam2junc_err)||(-1 == rm_SNPmerLearned_err))
  //       cout << "error in remove all intermediate files" << endl;
  // 	}

  //   string mv_headerSectionFile_cmd = "mv " + outputDirStr + "/headSectionInfo " + outputDirStr + "/logs/";
  //   string mv_stats_cmd = "mv " + outputDirStr + "/stats.txt " + outputDirStr + "/logs/"; 
  //   const int mv_headerSectionFile_err = system(mv_headerSectionFile_cmd.c_str());
  //   const int mv_stats_err = system(mv_stats_cmd.c_str());
  //   if((-1 == mv_headerSectionFile_err)||(-1 == mv_stats_err))
  //     cout << "error in mv_headerSectionFile_cmd or mv_stats_cmd" << endl;

  //   string mv_rawSam_cmd = "mv " + outputDirStr + "/output.sam " + outputDirStr + "/sam2junc/output.raw.sam";
  //   string removeFalseJuncSam_cmd = "removeFalseJuncSam " + indexStr + " " + outputDirStr 
  //     + "/sam2junc/output.alignInferJunc_classified_filterOut.txt " + outputDirStr + "/sam2junc/output.raw.sam "
  //     + outputDirStr + "/output.sam " + int_to_str(threads_num);
  //   const int mv_rawSam_err = system(mv_rawSam_cmd.c_str());
  //   const int removeFalseJuncSam_err = system(removeFalseJuncSam_cmd.c_str());
  //   if((-1 == mv_rawSam_err)||(-1 == removeFalseJuncSam_err))
  //     cout << "error in mv_rawSam_cmd or removeFalseJuncSam_cmd" << endl; 
  // }

  #ifdef MPS_FUSION_POST_NEW
  nowtime = time(NULL);
  local = localtime(&nowtime);
  tmpTimeStr = asctime(local);
  tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
  cout << endl << tmpTimeStr << "... fusion post analysis starts ......" << endl ;  
  log_ofs << endl << tmpTimeStr << "... fusion post analysis starts ......" << endl ; 
  runtime_log_ofs << endl << tmpTimeStr << "... fusion post analysis starts ......" << endl;
  int fusion_post_sup_min = optionInfo->return_fusion_post_sup_min();
  string fusion_post_sup_min_str = int_to_str(fusion_post_sup_min);
  string fusion_post_thread_num_str = int_to_str(threads_num);
  string fusion_post_gtf_path = optionInfo->return_fusion_post_formatted_gtf_path();
  string fusion_post_paralog_gene_path = optionInfo->return_fusion_post_paralog_gene_path();
  string fusion_post_output_path = outputDirStr + "/fusion_post/";
  string fusion_post_cmd = "./mps3_fusion_post_new " + indexStr 
    + " " + outputDirStr 
    + "/ " + fusion_post_gtf_path 
    + " " + fusion_post_thread_num_str 
    + " " + fusion_post_sup_min_str 
    + " " + fusion_post_paralog_gene_path 
    + " " + fusion_post_output_path;
  int fusion_post_err = system(fusion_post_cmd.c_str());
  if(-1 == fusion_post_err)
    cout << "error in fusion_post_cmd" << endl;
  nowtime = time(NULL);
  local = localtime(&nowtime);
  tmpTimeStr = asctime(local);
  tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
  cout << endl << tmpTimeStr << "... fusion post analysis ends ......" << endl ;  
  log_ofs << endl << tmpTimeStr << "... fusion post analysis ends ......" << endl ; 
  runtime_log_ofs << endl << tmpTimeStr << "... fusion post analysis ends ......" << endl;
  #endif


	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "... all jobs done ......" << endl << endl; 
	log_ofs << endl << tmpTimeStr << "... all jobs done ......" << endl << endl;
	runtime_log_ofs << endl << tmpTimeStr << "... all jobs done ......" << endl << endl << endl;

	log_ofs.close();
	runtime_log_ofs.close();
    return 0;
} //end main
