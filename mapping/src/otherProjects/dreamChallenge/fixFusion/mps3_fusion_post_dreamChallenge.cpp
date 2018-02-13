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

using namespace std;

time_t nowtime;
struct tm *local;	
string tmpTimeStr;

int main(int argc, char** argv)
{
	if(argc != 8)
	{
		cout << "Executable inputIndexFolderPath inputOriMPS3resultsFolder inputFormattedGTFpath threadsNum supNumMin paralogGeneFile outputFolder" << endl;
		exit(1);
	}
	////////////////////////////////////////////////////////////
	////////////////////   PARAMETERS  /////////////////////////
	////////////////////////////////////////////////////////////
	string offsetWhenComparingWithGTF = "5";
	string checkStrandedOrNotWhenComparingWithGTF = "Y";

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "creating folder ......" << endl;
	string outputFolderStr = argv[7];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	string inputIndexFolderPath = argv[1];
	string inputOriMPS3resultsFolder = argv[2];
	string inputFormattedGtf = argv[3];	
	string threadsNumStr = argv[4];
	//int threadsNum = atoi(threadsNumStr.c_str());
	string supNumMinStr = argv[5];
	int supNumMin = atoi(supNumMinStr.c_str());
	string inputReformattedParalogGeneGroupFile = argv[6];
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////    STEP 1:    SAM 2 JUNC        ////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 1:    SAM 2 JUNC ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 1:    SAM 2 JUNC ......" << endl;

	string sam2juncForFusionFilter_folderPath = outputFolderStr + "sam2juncForFusionFilter";
	string oriOutputSam = inputOriMPS3resultsFolder + "/phase1_output/completePair/completePair.sam";
	string sam2junc_forFusionFilter_cmd = "sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify " 
		+ inputIndexFolderPath + " " + threadsNumStr + " " + sam2juncForFusionFilter_folderPath + " " + oriOutputSam;
	cout << "sam2junc cmd: " << sam2junc_forFusionFilter_cmd << endl << endl;
	log_ofs << "sam2junc cmd: " << sam2junc_forFusionFilter_cmd << endl << endl;
	system(sam2junc_forFusionFilter_cmd.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 1:    SAM 2 JUNC ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 1:    SAM 2 JUNC ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////    STEP 2:   GLOBAL DETECTION    ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 2:   GLOBAL DETECTION ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 2:   GLOBAL DETECTION......" << endl;

	string inputIncompletePairSamFile = inputOriMPS3resultsFolder + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	string globalDetection_folderPath = outputFolderStr + "globalDetection";
	string globalDetection_cmd = "detection_globalMap " + inputIndexFolderPath + " " + inputFormattedGtf 
		+ " " + threadsNumStr + " " + inputIncompletePairSamFile + " " + globalDetection_folderPath;
	cout << "globalDetection cmd: " << globalDetection_cmd << endl << endl;
	log_ofs << "globalDetection cmd: " << globalDetection_cmd << endl << endl;
	system(globalDetection_cmd.c_str());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 2:   GLOBAL DETECTION ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 2:   GLOBAL DETECTION ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////   STEP 3: GENERATE FUSION JUNC FROM BREAK POINT FILE   ////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 3: GENERATE FUSION JUNC FROM BREAK POINT FILE ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 3: GENERATE FUSION JUNC FROM BREAK POINT FILE ......" << endl;

	string inputAdjustedFusionBreakPointFilePath = outputFolderStr + "globalDetection/fusion.breakPoint";
	string fusionJunc_file = inputAdjustedFusionBreakPointFilePath + ".junc.raw";
	string generateFusionJuncFromBreakPointFile_cmd = "generateFusionJuncInfoFromAdjustedBreakPointFile "
		+ inputIndexFolderPath + " " + inputAdjustedFusionBreakPointFilePath + " " + fusionJunc_file;
	cout << "generateFusionJuncFromBreakPointFile cmd: " << generateFusionJuncFromBreakPointFile_cmd << endl << endl;
	log_ofs << "generateFusionJuncFromBreakPointFile cmd: " << generateFusionJuncFromBreakPointFile_cmd << endl << endl;
	system(generateFusionJuncFromBreakPointFile_cmd.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 3: GENERATE FUSION JUNC FROM BREAK POINT FILE ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 3: GENERATE FUSION JUNC FROM BREAK POINT FILE ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////  STEP 4:   SEPARATE FUSION JUNC AS STRANDED OR NOT    /////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 4:   SEPARATE FUSION JUNC AS STRANDED OR NOT ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 4:   SEPARATE FUSION JUNC AS STRANDED OR NOT ......" << endl;

	string fusionJunc_folder = outputFolderStr + "/fusionJunc";
	string mkdir_fusionJunc_folder = "mkdir " + fusionJunc_folder;
	system(mkdir_fusionJunc_folder.c_str());

	string fusionJunc_separate_prefix = fusionJunc_folder + "/raw";
	string separateFusionJunc_strandedOrNot_cmd = "separateFusionJuncStrandedOrNot " + fusionJunc_file 
		+ " " + fusionJunc_separate_prefix;
	cout << "separateFusionJunc_strandedOrNot cmd: " << separateFusionJunc_strandedOrNot_cmd << endl << endl;
	log_ofs << "separateFusionJunc_strandedOrNot cmd: " << separateFusionJunc_strandedOrNot_cmd << endl << endl;
	system(separateFusionJunc_strandedOrNot_cmd.c_str());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 4:   SEPARATE FUSION JUNC AS STRANDED OR NOT......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 4:   SEPARATE FUSION JUNC AS STRANDED OR NOT......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////  STEP 5:  FILTER STRANDED FUSION JUNC       ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 5:  FILTER STRANDED FUSION JUNC ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 5:  FILTER STRANDED FUSION JUNC  ......" << endl;

	string sam2juncForFusionFilter_file = sam2juncForFusionFilter_folderPath + "/output.alignInferJunc.txt";
	string fusionJunc_stranded = fusionJunc_separate_prefix + ".stranded";
	string fusionJunc_stranded_filter_folder = fusionJunc_stranded + ".filter";
	string filterFusionJunc_stranded_cmd = "filterFusionJunc_anchorSeqSimilarity " + inputIndexFolderPath + " "
		+ fusionJunc_stranded + " " + fusionJunc_stranded_filter_folder + " " + sam2juncForFusionFilter_file;
	cout << "filterFusionJunc_stranded cmd: " << filterFusionJunc_stranded_cmd << endl << endl;
	log_ofs << "filterFusionJunc_stranded cmd: " << filterFusionJunc_stranded_cmd << endl << endl;
	system(filterFusionJunc_stranded_cmd.c_str());			

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 5:  FILTER STRANDED FUSION JUNC  ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 5:  FILTER STRANDED FUSION JUNC  ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////   STEP 6:  FILTER NON-STRANDED FUSION JUNC      ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 6:  FILTER NON-STRANDED FUSION JUNC......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 6:  FILTER NON-STRANDED FUSION JUNC ......" << endl;

	string fusionJunc_nonStranded = fusionJunc_separate_prefix + ".nonStranded";
	string fusionJunc_nonStranded_filter_folder = fusionJunc_nonStranded + ".filter";
	string filterFusionJunc_nonStranded_cmd = "filterNonStrandedFusionJunc_anchorSeqSimilarity " 
		+ inputIndexFolderPath + " " + fusionJunc_nonStranded + " " + fusionJunc_nonStranded_filter_folder
		+ " " + sam2juncForFusionFilter_file;
	cout << "filterFusionJunc_nonStranded cmd: " << filterFusionJunc_nonStranded_cmd << endl << endl;
	log_ofs << "filterFusionJunc_nonStranded cmd: " << filterFusionJunc_nonStranded_cmd << endl << endl;
	system(filterFusionJunc_nonStranded_cmd.c_str());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 6:  FILTER NON-STRANDED FUSION JUNC......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 6:  FILTER NON-STRANDED FUSION JUNC......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////         STEP 7:  CAT PASSED JUNC            ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 7:  CAT PASSED JUNC......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 7:  CAT PASSED JUNC ......" << endl;

	string passed_fusionJunc_stranded = fusionJunc_stranded_filter_folder + "/fusionJunc_classified_pass.txt";
	string passed_fusionJunc_nonStranded = fusionJunc_nonStranded_filter_folder + "/fusionJunc_classified_pass.txt";
	string passed_fusionJunc_bothStrandedOrNot = fusionJunc_folder + "/fusionJunc_pass_bothStrandedOrNot.txt";
	string cat_passed_fusionJunc_bothStrandedOrNot_cmd = "cat " + passed_fusionJunc_stranded + " "
		+ passed_fusionJunc_nonStranded + " > " + passed_fusionJunc_bothStrandedOrNot;
	cout << "cat_passed_fusionJunc_bothStrandedOrNot cmd: " << cat_passed_fusionJunc_bothStrandedOrNot_cmd << endl << endl;
	log_ofs << "cat_passed_fusionJunc_bothStrandedOrNot cmd: " << cat_passed_fusionJunc_bothStrandedOrNot_cmd << endl << endl;
	system(cat_passed_fusionJunc_bothStrandedOrNot_cmd.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 7:  CAT PASSED JUNC ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 7:  CAT PASSED JUNC ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////       STEP 8:  CMP PASSED JUNC 2 GTF        ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 8:  CMP PASSED JUNC 2 GTF ......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 8:  CMP PASSED JUNC 2 GTF ......" << endl;

	string passedJunc_cmp2gtf_folder = outputFolderStr + "/fusionJunc_pass_cmp2ann";
	string passedJunc_cmp2gtf_cmd = "compareFusionJuncResultsWithGTFannotation " + inputIndexFolderPath
		+ " " + inputFormattedGtf + " " + passed_fusionJunc_bothStrandedOrNot + " " + passedJunc_cmp2gtf_folder + " " + offsetWhenComparingWithGTF 
		+ " " + checkStrandedOrNotWhenComparingWithGTF;
	cout << "passedJunc_cmp2gtf cmd: " << passedJunc_cmp2gtf_cmd << endl << endl;
	log_ofs << "passedJunc_cmp2gtf cmd: " << passedJunc_cmp2gtf_cmd << endl << endl;
	system(passedJunc_cmp2gtf_cmd.c_str());		

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 8:  CMP PASSED JUNC 2 GTF ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 8:  CMP PASSED JUNC 2 GTF ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////      STEP 9:  FILTERING_SUPNUM_REPEAT       ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "start to do STEP 9:  FILTERING_SUPNUM_REPEAT......" << endl;
	log_ofs << endl << tmpTimeStr << "start to do STEP 9:  FILTERING_SUPNUM_REPEAT ......" << endl;

	string filtering_supNum_repeat_cmd = "filtering_supNum_repeat " + supNumMinStr + " " + inputReformattedParalogGeneGroupFile + " "
		+ passedJunc_cmp2gtf_folder + "/fusionJunc_geneInfo_interGene.txt " + passedJunc_cmp2gtf_folder + "/fusionJunc_geneInfo_interGene";
	cout << "filtering_supNum_repeat cmd: " << filtering_supNum_repeat_cmd << endl << endl;
	log_ofs << "filtering_supNum_repeat cmd: " << filtering_supNum_repeat_cmd << endl << endl;
	system(filtering_supNum_repeat_cmd.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	tmpTimeStr = asctime(local);
	tmpTimeStr = "[" + tmpTimeStr.substr(0, tmpTimeStr.length()-1) + "] ";
	cout << endl << tmpTimeStr << "end of STEP 9:  FILTERING_SUPNUM_REPEAT ......" << endl;
	log_ofs << endl << tmpTimeStr << "end of STEP 9:  FILTERING_SUPNUM_REPEAT ......" << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////      STEP 9:  FILTERING_SUPNUM_REPEAT       ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////	
	cout << endl << tmpTimeStr << "all jobs done ......" << endl;
	log_ofs << endl << tmpTimeStr << "all jobs done ......" << endl;	
	log_ofs.close();
	return 0;	
}