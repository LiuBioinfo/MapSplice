#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFolderList outputScriptPath" << endl;
		exit(1);
	}
	string inputFolderListFile = argv[1];
	ifstream folderList_ifs(inputFolderListFile.c_str());
	vector<string> folderPathVec;
	while(!folderList_ifs.eof())
	{
		string tmpStr;
		getline(folderList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		folderPathVec.push_back(tmpStr);
	}
	folderList_ifs.close();

	string outputScriptPath = argv[2];
	ofstream script_ofs(outputScriptPath.c_str());
	for(int tmp = 0; tmp < folderPathVec.size(); tmp++)
	{
		string tmpFolderPath = folderPathVec[tmp];
		tmpFolderPath += "/";
		string tmpHeadSectionInfo = tmpFolderPath + "headSectionInfo";
		string tmpAlignCompleteRead = tmpFolderPath + "phase1_output/completePair/completePair.sam";
		string OutputSamFile_oneEndMapped = tmpFolderPath + "phase2_output/oneEndUnmapped.pairedComplete.sam";
		string OutputSamFile_fixHeadTail_complete_pair = tmpFolderPath + "phase2_output/fixHeadTail_complete_pair.sam";
		string OutputSamFile_fixHeadTail_incomplete_pair = tmpFolderPath + "phase2_output/fixHeadTail_incomplete_pair.sam";
		string OutputSamFile_fixHeadTail_complete_unpair = tmpFolderPath + "phase2_output/fixHeadTail_complete_unpair.sam"; 
		string OutputSamFile_fixHeadTail_incomplete_unpair = tmpFolderPath + "phase2_output/fixHeadTail_incomplete_unpair.sam";
		string OutputSamFile_oneEndMapped_unpair = tmpFolderPath + "phase2_output/oneEndUnmapped.unpaired.sam";
		string tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = tmpFolderPath + "phase1_output/repeat_region/bothEndsUnmapped_mappedToRepeatRegion.sam";
		string tmpAlignBothEndsUnmapped_lowScore = tmpFolderPath + "phase1_output/bothEndsUnmapped/bothEndsUnmapped_lowScore.sam";
		string OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore = tmpFolderPath + "phase2_output/oneEndUnmapped.bothEndsUnmapped_lowScore.sam";
		string OutputSamFile_fixHeadTail_pair_lowScore = tmpFolderPath + "phase2_output/fixHeadTail_pair_lowScore.sam";
		string tmpAlignBothEndsUnmapped  = tmpFolderPath + "phase1_output/bothEndsUnmapped/bothEndsUnmapped.sam";
		string finalOutputSam = tmpFolderPath + "output.sam";
		string cat_cmd = "cat " + tmpHeadSectionInfo + " " + tmpAlignCompleteRead 
			+ " " + OutputSamFile_oneEndMapped + " " + OutputSamFile_fixHeadTail_complete_pair
			+ " " + OutputSamFile_fixHeadTail_incomplete_pair + " " + OutputSamFile_fixHeadTail_complete_unpair 
			+ " " + OutputSamFile_fixHeadTail_incomplete_unpair + " " + OutputSamFile_oneEndMapped_unpair 
			+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + " " + tmpAlignBothEndsUnmapped_lowScore
			+ " " + OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore 
			+ " " + OutputSamFile_fixHeadTail_pair_lowScore + " " + tmpAlignBothEndsUnmapped + " > " + finalOutputSam;
		script_ofs << cat_cmd << endl;
	}
	script_ofs.close();
	return 0;
}