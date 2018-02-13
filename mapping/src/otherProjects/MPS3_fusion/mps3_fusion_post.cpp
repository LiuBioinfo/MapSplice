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

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputResultsFolderPath threadsNum inputGtf" << endl;
		exit(1);
	}

	time_t nowtime;
	struct tm *local;

	string inputIndexFolderPath = argv[1];	
	string inputResultsFolderPath = argv[2];
	string threadsNumStr = argv[3];
	int threadsNum = atoi(threadsNumStr.c_str());
	string inputGtf = argv[4];
	string inputResultsFolderPath_outputSam = inputResultsFolderPath + "/output.sam";
	string inputResultsFolderPath_phase2 = inputResultsFolderPath + "/phase2_output/";
	string inputResultsFolderPath_phase2_incompletePairSam = inputResultsFolderPath_phase2 + "/fixHeadTail_incomplete_pair.sam";
	string inputResultsFolderPath_phase2_globalMap = inputResultsFolderPath_phase2 + "/globalMap";
	string inputResultsFolderPath_phase2_sam2juncForFusionFilter = inputResultsFolderPath_phase2 + "/sam2juncForFusionFilter";
	string inputResultsFolderPath_phase2_fusionDetection = inputResultsFolderPath_phase2_globalMap + "/fusionDetection";

	string sam2junc_forFusionFilter_cmd = "./sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify " + inputIndexFolderPath
		+ " " + threadsNumStr + " " + inputResultsFolderPath_phase2_sam2juncForFusionFilter + " " + inputResultsFolderPath_outputSam;
	cout << "sam2junc_forFusionFilter_cmd: " << endl << sam2junc_forFusionFilter_cmd << endl;
	system(sam2junc_forFusionFilter_cmd.c_str());

	string globalMap_cmd = "./globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion " + inputIndexFolderPath + " "
		+ inputResultsFolderPath_phase2_incompletePairSam + " " + inputResultsFolderPath_phase2_globalMap + " " + threadsNumStr
		+ " fastq 300000 Y";
	cout << "globalMap_cmd: " << endl << globalMap_cmd << endl;
	system(globalMap_cmd.c_str());

	string generateFusionJunc_filter_cmp2gtf_cmd = "./generateFusionJunc_filter_compare2gtf " + inputIndexFolderPath + " "
		+ inputGtf + " " + inputResultsFolderPath_phase2_globalMap + "/fusionReadsWithBreakPoint_adjusted_withMapRange.txt "
		+ inputResultsFolderPath_phase2_sam2juncForFusionFilter + "/output.alignInferJunc.txt 5 N " + inputResultsFolderPath_phase2_fusionDetection;
	cout << "generateFusionJunc_filter_cmp2gtf_cmd: " << endl << generateFusionJunc_filter_cmp2gtf_cmd << endl;
	system(generateFusionJunc_filter_cmp2gtf_cmd.c_str());

	string separateFusionJuncStrandedOrNot_cmd = "./separateFusionJuncStrandedOrNot " + inputResultsFolderPath_phase2_fusionDetection
		+ "/fusionJunc_pass_compare2gtf/fusionJunc_geneInfo_interGene.txt " + inputResultsFolderPath_phase2_fusionDetection
		+ "/fusionJunc_pass_compare2gtf/fusionJunc_geneInfo_interGene";
	cout << "separateFusionJuncStrandedOrNot_cmd: " << separateFusionJuncStrandedOrNot_cmd << endl;
	system(separateFusionJuncStrandedOrNot_cmd.c_str());

	// string remapAndCountEncompRead_cmd = "./remapAndCountEncompRead " + inputIndexFolderPath + " " + inputResultsFolderPath_phase2_fusionDetection 
	// 	+ "/fusionJunc_ori.stranded.filterBasedOnAnchorSeqSimi/fusionJunc_classified_pass.txt " + inputResultsFolderPath_phase2 + " " 
	// 	+ inputResultsFolderPath_phase2_globalMap + "/nonFusion.sam " + inputResultsFolderPath_phase2_fusionDetection 
	// 	+ "fusionJunc_ori.stranded.filterBasedOnAnchorSeqSimi/remapAndCountEncompRead " + threadsNumStr + " fastq stranded";
	// cout << "remapAndCountEncompRead_cmd: " << endl << remapAndCountEncompRead_cmd << endl;
	// system(remapAndCountEncompRead_cmd.c_str());
	return 0;
}