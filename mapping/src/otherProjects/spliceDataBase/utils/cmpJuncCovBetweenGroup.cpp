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
#include <hash_map>
#include <map>
#include <set>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

void generateTwoSpliceDataBaseFormatFileVec_groupedFilePathFile(
	vector<string>& groupNameVec, vector< vector<string> >& spliceDataBaseFormatFileVecVec_groupWise,
	vector<string>& spliceDataBaseFormatFileVecVec_total, string& groupedFilePathFile)
{
	cout << "generateTwoSpliceDataBaseFormatFileVec_groupedFilePathFile starts ......" << endl;
	ifstream groupedFilePath_ifs(groupedFilePathFile.c_str());
	vector<string> strVec;
	while(!groupedFilePath_ifs.eof())
	{
		string tmpStr;
		getline(groupedFilePath_ifs, tmpStr);
		if(tmpStr == "")
			break;
		if(tmpStr.substr(0,1) == "@")
			groupNameVec.push_back(tmpStr);
		else
			spliceDataBaseFormatFileVecVec_total.push_back(tmpStr);
		strVec.push_back(tmpStr);
	}
	groupedFilePath_ifs.close();
	int groupNum = groupNameVec.size();
	cout << "groupNum: " << groupNum << endl;
	for(int tmpGroup = 0; tmpGroup < groupNum; tmpGroup ++)
	{
		vector<string> tmpGroupFileVec;
		spliceDataBaseFormatFileVecVec_groupWise.push_back(tmpGroupFileVec);
	}
	int tmpGroupIndex = -1;
	for(int tmp = 0; tmp < strVec.size(); tmp++)
	{
		string tmpStr = strVec[tmp];
		cout << "tmpStr: " << tmpStr << endl;
		if(tmpStr.substr(0,1) == "@")
			tmpGroupIndex ++;
		else
			spliceDataBaseFormatFileVecVec_groupWise[tmpGroupIndex].push_back(tmpStr);
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputGroupedFilePathFilePath outputFolder" << endl;
		exit(1);
	}
    string outputDirStr = argv[3];
   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());	
   	string settingsLogStr = outputDirStr + "/log.txt";
   	ofstream log_ofs(settingsLogStr.c_str());

	string indexFolderPath = argv[1];
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	Index_Info* indexInfo = new Index_Info(indexParameterFileStr);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	nowtime = time(NULL);
	local = localtime(&nowtime);	
	cout << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "... end of initiating indexInfo" << endl;

	string groupedFilePathFilePath = argv[2];
	vector<string> groupNameVec;
	vector< vector<string> > spliceDataBaseFormatFileVecVec_groupWise;
	vector<string> spliceDataBaseFormatFileVecVec_total;
	generateTwoSpliceDataBaseFormatFileVec_groupedFilePathFile(
		groupNameVec, spliceDataBaseFormatFileVecVec_groupWise, 
		spliceDataBaseFormatFileVecVec_total, groupedFilePathFilePath);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initaite merged alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_merged = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_merged->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_merged->insertJuncFromSpliceDataBaseFormatFileVec_chrNamePosOnly(
		spliceDataBaseFormatFileVecVec_total, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initaite alignInferJuncHashInfoVec" << endl;
	vector< vector<AlignInferJunctionHash_Info*> > alignInferJuncHashInfoVecVec;
	int groupNum = groupNameVec.size();
	for(int tmpGroupIndex = 0; tmpGroupIndex < groupNum; tmpGroupIndex ++)
	{
		cout << "tmpGroupIndex: " << tmpGroupIndex << endl;
		cout << "tmpGroupName: " << groupNameVec[tmpGroupIndex] << endl;
		vector<AlignInferJunctionHash_Info*> tmpAlignInferJuncHashInfoVec;
		int juncHashNumInTmpGroup = spliceDataBaseFormatFileVecVec_groupWise[tmpGroupIndex].size();
		for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncHashNumInTmpGroup; tmpJuncHashIndex ++)
		{
			string tmpSpliceDataBaseFormatFilePath = (spliceDataBaseFormatFileVecVec_groupWise[tmpGroupIndex])[tmpJuncHashIndex];
			cout << "tmpJuncHashIndex: " << tmpJuncHashIndex << endl;
			cout << "tmpSpliceDataBaseFormatFilePath: " << tmpSpliceDataBaseFormatFilePath << endl;
			AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
			tmpJuncHash->initiateAlignInferJunctionHashInfo(chromNum);
			tmpJuncHash->insertJuncFromSpliceDataBaseFormatFile_chrNamePos_supportNum_anchorSize_branchRatio_hapFlkStr(
				tmpSpliceDataBaseFormatFilePath, indexInfo);
			tmpAlignInferJuncHashInfoVec.push_back(tmpJuncHash);
		}
		alignInferJuncHashInfoVecVec.push_back(tmpAlignInferJuncHashInfoVec);
	}

	string output_groupWiseJuncSupNum_file = outputDirStr + "/groupWiseJuncSupNum.txt";
	ofstream groupWiseJuncSupNum_ofs(output_groupWiseJuncSupNum_file.c_str()); 
	int totalJuncNum = alignInferJunctionHashInfo_merged->returnAlignInferInfoVecSize();
	for(int tmpJunc = 0; tmpJunc < totalJuncNum; tmpJunc ++)
	{
		int tmpJunc_chrNameInt = alignInferJunctionHashInfo_merged->returnAlignInferInfo_chrNameInt(tmpJunc);
		string tmpJunc_chrName = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
		int tmpJunc_donerPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_donerEndPos(tmpJunc);
		int tmpJunc_acceptorPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_acceptorStartPos(tmpJunc);
		int tmpJunc_totalSupNum = 0;
		vector<int> tmpJunc_groupWiseSupNumVec;
		for(int tmpGroupIndex = 0; tmpGroupIndex < groupNum; tmpGroupIndex ++)
		{
			int tmpJunc_tmpGroup_supNumSum = 0;
			int juncHashNumInTmpGroup = spliceDataBaseFormatFileVecVec_groupWise[tmpGroupIndex].size();
			for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncHashNumInTmpGroup; tmpJuncHashIndex ++)
			{				
				int tmpJunc_tmpJuncHash_supNum = (alignInferJuncHashInfoVecVec[tmpGroupIndex])[tmpJuncHashIndex]->searchAndReturnAlignInferJuncHashSupNum(
					tmpJunc_chrNameInt, tmpJunc_donerPos, tmpJunc_acceptorPos);
				tmpJunc_totalSupNum += tmpJunc_tmpJuncHash_supNum;
				tmpJunc_tmpGroup_supNumSum += tmpJunc_tmpJuncHash_supNum;
			}
			tmpJunc_groupWiseSupNumVec.push_back(tmpJunc_tmpGroup_supNumSum);		
		}
		groupWiseJuncSupNum_ofs << tmpJunc_chrName << "\t" << tmpJunc_donerPos << "\t" << tmpJunc_acceptorPos << "\t"
			<< tmpJunc_totalSupNum;
		for(int tmpGroupIndex = 0; tmpGroupIndex < groupNum; tmpGroupIndex ++)
			groupWiseJuncSupNum_ofs << "\t" << tmpJunc_groupWiseSupNumVec[tmpGroupIndex];
		groupWiseJuncSupNum_ofs << endl;
	}
	groupWiseJuncSupNum_ofs.close();
	for(int tmpGroupIndex = 0; tmpGroupIndex < groupNum; tmpGroupIndex ++)
	{
		int juncHashNumInTmpGroup = spliceDataBaseFormatFileVecVec_groupWise[tmpGroupIndex].size();
		for(int tmpJuncHashIndex = 0; tmpJuncHashIndex < juncHashNumInTmpGroup; tmpJuncHashIndex ++)
			delete alignInferJuncHashInfoVecVec[tmpGroupIndex][tmpJuncHashIndex];
	}
	delete alignInferJunctionHashInfo_merged;
	log_ofs.close();
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	return 0;
}