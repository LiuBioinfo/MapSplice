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
//#include "../../../general/alignInferJunctionHash_info.h"
#include "../../../general/alignInferJunctionHash_info_vec.h"
#include "../../../general/baseContiguousReadCount_info_vec.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if((argc < 5)||(argc > 7))
	{
		cout << "Executable <IndexPath> <threadNum> <inputSAM> <outputJuncFile> (<SNP_hap1_or_singlSNPfile>) (<SNP_hap2>)" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to initiate indexInfo" << endl;
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	//ifstream parameter_ifs(indexParameterFileStr.c_str());
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

	cout << "start to load single SNP file or two SNP files from both haplotypes" << endl;
	cout << "start initiating indexInfo_withSNP" << endl;
	Index_Info* indexInfo_withSNP = new Index_Info(indexParameterFileStr);
	cout << "start initiating indexInfo_hap1" << endl;
	Index_Info* indexInfo_hap1 = new Index_Info(indexParameterFileStr);
	cout << "start initiating indexInfo_hap2" << endl; 
	Index_Info* indexInfo_hap2 = new Index_Info(indexParameterFileStr);
	if(argc == 6) // load single SNP file
	{
		string singleSNPfile = argv[5];
		indexInfo_withSNP->readGenome(chrom);
		indexInfo_withSNP->initiate();
		indexInfo_withSNP->initiateChrNameIndexArray(1000);
		indexInfo_withSNP->insertSNP2chromStr(singleSNPfile);
	}
	else if(argc == 7) // load two SNP files from both haplotypes 
	{
		string hap1SNPfile = argv[5];
		indexInfo_hap1->readGenome(chrom);
		indexInfo_hap1->initiate();
		indexInfo_hap1->initiateChrNameIndexArray(1000);
		indexInfo_hap1->insertSNP2chromStr(hap1SNPfile);

		string hap2SNPfile = argv[6];
		indexInfo_hap2->readGenome(chrom);
		indexInfo_hap2->initiate();
		indexInfo_hap2->initiateChrNameIndexArray(1000);
		indexInfo_hap2->insertSNP2chromStr(hap2SNPfile);
	}
	else
	{}

	string threadNumStr = argv[2];
	int threadNum = atoi(threadNumStr.c_str());
	cout << "threadNum: " << threadNum << endl;
	int alignInferJuncHashInfoVecSize = threadNum;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initaite merged alignInferJunctionHashInfo " << endl;
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_merged = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_merged->initiateAlignInferJunctionHashInfo(chromNum);	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initaite alignInferJunctionHashInfo vec " << endl;	
	AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec = new AlignInferJunctionHash_Info_Vec();
	alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(alignInferJuncHashInfoVecSize, chromNum);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to insert SJ into SJmap from SAM file" << endl;
	// extract SJs from alignment file
	string inputSAMfile = argv[3];
	alignInferJunctionHashInfoVec->insertJuncFromAlignmentFile_chrNamePos_supportNum_parallel(
		inputSAMfile, indexInfo, alignInferJuncHashInfoVecSize);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to merge alignInferJuncHashInfo in vec" << endl;
	alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(
		alignInferJunctionHashInfo_merged, indexInfo);
	alignInferJunctionHashInfoVec->freeMemory();
	delete alignInferJunctionHashInfoVec;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to initiate alignInferJunctionHashInfo ...." << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	alignInferJunctionHashInfo_merged->convert2SJhashInfo(SJhashInfo, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to generate spliceSiteContiguousReadCount ...." << endl;
	BaseContiguousReadCount_Info_Vec* baseContiguousReadCountInfoVec = new BaseContiguousReadCount_Info_Vec(); 
	baseContiguousReadCountInfoVec->initiateBaseContiguousReadCountInfoVec(threadNum, chromNum);
	baseContiguousReadCountInfoVec->initiateBase_spliceSiteFromAlignInferJuncHash(alignInferJunctionHashInfo_merged);
	baseContiguousReadCountInfoVec->generatePosAreaMap_withStoredTotalBase();
	baseContiguousReadCountInfoVec->update_alignmentFile_parallel(inputSAMfile, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to merge spliceSiteContiguousReadCount in vec ...." << endl;
	BaseContiguousReadCount_Info* baseContiguousReadCountInfo_merged = new BaseContiguousReadCount_Info();
	baseContiguousReadCountInfo_merged->initaite_chromNum(chromNum);
	baseContiguousReadCountInfo_merged->initiateBase_spliceSiteFromAlignInferJuncHash(alignInferJunctionHashInfo_merged);
	baseContiguousReadCountInfo_merged->generatePosAreaMap_withStoredTotalBase();
	baseContiguousReadCountInfoVec->mergeBaseContiguousReadCountInfoVec2one(baseContiguousReadCountInfo_merged);
	baseContiguousReadCountInfoVec->freeMemory();
	delete baseContiguousReadCountInfoVec;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...start to output junc with branching ratio and flank strings " << endl;
	string outputJuncFile = argv[4];
	ofstream branchingRatioJunc_ofs(outputJuncFile.c_str());
	int juncNum = alignInferJunctionHashInfo_merged->returnAlignInferInfoVecSize();
	for(int tmpJunc = 0; tmpJunc < juncNum; tmpJunc ++)
	{
		//cout << "tmpJunc: " << tmpJunc + 1 << endl;
		int tmpJunc_chrNameInt = alignInferJunctionHashInfo_merged->returnAlignInferInfo_chrNameInt(tmpJunc);
		int tmpJunc_donerEndPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_donerEndPos(tmpJunc);
		int tmpJunc_acceptorStartPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_acceptorStartPos(tmpJunc);
		int tmpJunc_supportNum = alignInferJunctionHashInfo_merged->returnAlignInferInfo_supportNum(tmpJunc);
		//cout << "tmpJunc_chrNameInt: " << tmpJunc_chrNameInt << endl;
		//cout << "tmpJunc_donerEndPos: " << tmpJunc_donerEndPos << endl;
		//cout << "tmpJunc_acceptorStartPos: " << tmpJunc_acceptorStartPos << endl;
		double tmpJunc_donerBrranchingRatio = 1.0;
		double tmpJunc_acceptorBranchingRatio = 1.0;
		
		// count contiguous reads at splice sites (donerPos & acceptorPos-1)
		int tmpJunc_contiguousReadCount_donerSite = baseContiguousReadCountInfo_merged->return_contiguousReadCount(
			tmpJunc_chrNameInt, tmpJunc_donerEndPos);
		//cout << "tmpJunc_contiguousReadCount_donerSite: " << tmpJunc_contiguousReadCount_donerSite << endl;
		int tmpJunc_contiguousReadCount_acceptorSite = baseContiguousReadCountInfo_merged->return_contiguousReadCount(
			tmpJunc_chrNameInt, tmpJunc_acceptorStartPos - 1);
		//cout << "tmpJunc_contiguousReadCount_acceptorSite: " << tmpJunc_contiguousReadCount_acceptorSite << endl;
		// check branching ratio at doner site, step 1 -- count total branching out junction sup #; step 2 -- total sup # / junc sup #
		vector<int> tmpAlterAcceptorSpliceSitePosVec;
		SJhashInfo->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpAlterAcceptorSpliceSitePosVec);
		int tmpAlterAcceptorSpliceSitePosVecSize = tmpAlterAcceptorSpliceSitePosVec.size();
		//cout << "tmpAlterAcceptorSpliceSitePosVecSize: " << tmpAlterAcceptorSpliceSitePosVecSize << endl;
		if(tmpAlterAcceptorSpliceSitePosVecSize == 0)
		{
			cout << "error in tmpAlterAcceptorSpliceSitePosVecSize: " << tmpAlterAcceptorSpliceSitePosVecSize << endl;
			exit(1);
		}
		else if(tmpAlterAcceptorSpliceSitePosVecSize == 1) // no other junc sharing the same doner site
			tmpJunc_donerBrranchingRatio = (double)tmpJunc_supportNum / (double)(tmpJunc_supportNum + tmpJunc_contiguousReadCount_donerSite);//1.0;
		else
		{
			int totalOtherJuncSupNum = 0;
			for(int tmpAlterJunc = 0; tmpAlterJunc < tmpAlterAcceptorSpliceSitePosVecSize; tmpAlterJunc ++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmpAlterJunc];
				if(tmpAlterAcceptorSpliceSitePos == tmpJunc_acceptorStartPos)
					continue;
				int tmpAlterJuncSupNum = alignInferJunctionHashInfo_merged->searchAndReturnAlignInferJuncHashSupNum(tmpJunc_chrNameInt, 
					tmpJunc_donerEndPos, tmpAlterAcceptorSpliceSitePos);
				if(tmpAlterJuncSupNum < 1)
				{
					cout << "error in tmpAlterJuncSupNum: " << tmpAlterJuncSupNum << endl;
					cout << "tmpAlterJunc: " << tmpJunc_chrNameInt << "\t" << tmpJunc_donerEndPos 
						<< "\t" << tmpAlterAcceptorSpliceSitePos << endl;
					exit(1);
				}
				totalOtherJuncSupNum += tmpAlterJuncSupNum;
			}
			tmpJunc_donerBrranchingRatio = (double)tmpJunc_supportNum / (double)(totalOtherJuncSupNum + tmpJunc_supportNum + tmpJunc_contiguousReadCount_donerSite);
		}
		//cout << "tmpJunc_donerBrranchingRatio: " << tmpJunc_donerBrranchingRatio << endl;
		// check branching ratio at acceptor site, step 1 -- count total branching out junction sup #; step 2 -- total sup # / junc sup #
		vector<int> tmpAlterDonerSpliceSitePosVec;
		SJhashInfo->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpJunc_chrNameInt, tmpJunc_acceptorStartPos, tmpAlterDonerSpliceSitePosVec);
		int tmpAlterDonerSpliceSitePosVecSize = tmpAlterDonerSpliceSitePosVec.size();
		//cout << "tmpAlterDonerSpliceSitePosVecSize: " << tmpAlterDonerSpliceSitePosVecSize << endl;
		if(tmpAlterDonerSpliceSitePosVecSize == 0)
		{
			cout << "error in tmpAlterDonerSpliceSitePosVecSize: " << tmpAlterDonerSpliceSitePosVecSize << endl;
			exit(1);			
		}
		else if(tmpAlterDonerSpliceSitePosVecSize == 1)
			tmpJunc_acceptorBranchingRatio = 1.0;
		else
		{
			int totalOtherJuncSupNum = 0;
			for(int tmpAlterJunc = 0; tmpAlterJunc < tmpAlterDonerSpliceSitePosVecSize; tmpAlterJunc ++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmpAlterJunc];
				if(tmpAlterDonerSpliceSitePos == tmpJunc_donerEndPos)
					continue;
				int tmpAlterJuncSupNum = alignInferJunctionHashInfo_merged->searchAndReturnAlignInferJuncHashSupNum(tmpJunc_chrNameInt,
					tmpAlterDonerSpliceSitePos, tmpJunc_acceptorStartPos);
				if(tmpAlterJuncSupNum < 1)
				{
					cout << "error in tmpAlterJuncSupNum: " << tmpAlterJuncSupNum << endl;
					cout << "tmpAlterJunc: " << tmpJunc_chrNameInt << "\t" << tmpAlterDonerSpliceSitePos
						<< "\t" << tmpJunc_acceptorStartPos << endl;
					exit(1);
				}
				totalOtherJuncSupNum += tmpAlterJuncSupNum;
			}
			tmpJunc_acceptorBranchingRatio = (double)tmpJunc_supportNum / (double)(totalOtherJuncSupNum + tmpJunc_supportNum + tmpJunc_contiguousReadCount_acceptorSite);
		}
		//cout << "tmpJunc_acceptorBranchingRatio: " << tmpJunc_acceptorBranchingRatio << endl;
		// generate tmpJunc_flankString_ref, generate tmpJunc_flankString_hap1, generate tmpJunc_flankString_hap2
		string tmpJunc_flankString_ref = indexInfo->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
		//cout << "tmpJunc_flankString_ref: " << tmpJunc_flankString_ref << endl;
		string tmpJunc_flankString_hap1, tmpJunc_flankString_hap2;
		if(argc == 5) // no SNP provided
		{
			tmpJunc_flankString_hap1 = tmpJunc_flankString_ref;
			tmpJunc_flankString_hap2 = tmpJunc_flankString_ref;
		}
		else if(argc == 6) // single SNP file provided
		{
			tmpJunc_flankString_hap1 = indexInfo_withSNP->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
			tmpJunc_flankString_hap2 = tmpJunc_flankString_hap1;
		}
		else if(argc == 7) // two SNP files from both haplotypes provided
		{
			tmpJunc_flankString_hap1 = indexInfo_hap1->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
			tmpJunc_flankString_hap2 = indexInfo_hap2->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
		}
		else
		{
			cout << "error in argc: " << argc << endl;
			exit(1);
		}
		// cout << "tmpJunc_flankString_hap1: " << tmpJunc_flankString_hap1 << endl;
		// cout << "tmpJunc_flankString_hap2: " << tmpJunc_flankString_hap2 << endl;
		// output junctions with branching ratio and flank strings
		branchingRatioJunc_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t" << tmpJunc_donerEndPos << "\t" << tmpJunc_acceptorStartPos 
			<< "\t" << tmpJunc_supportNum << "\t" << tmpJunc_donerBrranchingRatio << "\t" << tmpJunc_acceptorBranchingRatio 
			<< "\t" << tmpJunc_flankString_ref << "\t" << tmpJunc_flankString_hap1 << "\t" << tmpJunc_flankString_hap2 << endl;;
	}
	branchingRatioJunc_ofs.close();
	delete alignInferJunctionHashInfo_merged;
	delete SJhashInfo;
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	//parameter_ifs.close();
	delete indexInfo_withSNP;
	delete indexInfo_hap1;
	delete indexInfo_hap2;
	//delete alignInferJunctionHashInfo_merged;
	delete baseContiguousReadCountInfo_merged;
	return 0;
}