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

int main(int argc, char** argv)
{
	if((argc < 4)||(argc > 6))
	{
		cout << "Executable <IndexInput> <inputJuncFileWithAtLeastChrNamePosSupNum> <outputJuncFileWithChrNamePosSupNumBranchingRatioAtBothSites> (<SNP_hap1_or_singlSNPfile>) (<SNP_hap2>)" << endl;
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

	string inputJuncFile = argv[2];
	cout << "start to initiate alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_raw = new AlignInferJunctionHash_Info();
	juncHash_raw->initiateAlignInferJunctionInfo(chromNum);
	juncHash_raw->insertJuncFromJuncFile_chrNamePos_supportNum_inclusive(inputJuncFile, indexInfo);
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	juncHash_raw->convert2SJhashInfo(SJhashInfo, indexInfo);
	cout << "start to calculate branching ratio at splice sites for each junc" << endl;

	cout << "start to load single SNP file or two SNP files from both haplotypes" << endl;
	cout << "start initiating indexInfo_withSNP" << endl;
	Index_Info* indexInfo_withSNP = new Index_Info(indexParameterFileStr);
	cout << "start initiating indexInfo_hap1" << endl;
	Index_Info* indexInfo_hap1 = new Index_Info(indexParameterFileStr);
	cout << "start initiating indexInfo_hap2" << endl; 
	Index_Info* indexInfo_hap2 = new Index_Info(indexParameterFileStr);
	if(argc == 5) // load single SNP file
	{
		string singleSNPfile = argv[4];
		indexInfo_withSNP->readGenome(chrom);
		indexInfo_withSNP->initiate();
		indexInfo_withSNP->initiateChrNameIndexArray(1000);
		indexInfo_withSNP->insertSNP2chromStr(singleSNPfile);
	}
	else if(argc == 6) // load two SNP files from both haplotypes 
	{
		string hap1SNPfile = argv[4];
		indexInfo_hap1->readGenome(chrom);
		indexInfo_hap1->initiate();
		indexInfo_hap1->initiateChrNameIndexArray(1000);
		indexInfo_hap1->insertSNP2chromStr(hap1SNPfile);

		string hap2SNPfile = argv[5];
		indexInfo_hap2->readGenome(chrom);
		indexInfo_hap2->initiate();
		indexInfo_hap2->initiateChrNameIndexArray(1000);
		indexInfo_hap2->insertSNP2chromStr(hap2SNPfile);
	}
	else
	{}
	cout << "start to output junc with branching ratio and flank strings " << endl;
	string outputJuncFile = argv[3];
	ofstream branchingRatioJunc_ofs(outputJuncFile.c_str());
	int juncNum = juncHash_raw->returnAlignInferInfoVecSize();
	for(int tmpJunc = 0; tmpJunc < juncNum; tmpJunc ++)
	{
		int tmpJunc_chrNameInt = juncHash_raw->returnAlignInferInfo_chrNameInt(tmpJunc);
		int tmpJunc_donerEndPos = juncHash_raw->returnAlignInferInfo_donerEndPos(tmpJunc);
		int tmpJunc_acceptorStartPos = juncHash_raw->returnAlignInferInfo_acceptorStartPos(tmpJunc);
		int tmpJunc_supportNum = juncHash_raw->returnAlignInferInfo_supportNum(tmpJunc);
		double tmpJunc_donerBrranchingRatio = 1.0;
		double tmpJunc_acceptorBranchingRatio = 1.0;
		// check branching ratio at doner site, step 1 -- count total branching out junction sup #; step 2 -- total sup # / junc sup #
		vector<int> tmpAlterAcceptorSpliceSitePosVec;
		SJhashInfo->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpAlterAcceptorSpliceSitePosVec);
		int tmpAlterAcceptorSpliceSitePosVecSize = tmpAlterAcceptorSpliceSitePosVec.size();
		if(tmpAlterAcceptorSpliceSitePosVecSize == 0)
		{
			cout << "error in tmpAlterAcceptorSpliceSitePosVecSize: " << tmpAlterAcceptorSpliceSitePosVecSize << endl;
			exit(1);
		}
		else if(tmpAlterAcceptorSpliceSitePosVecSize == 1) // no other junc sharing the same doner site
			tmpJunc_donerBrranchingRatio = 1.0;
		else
		{
			int totalOtherJuncSupNum = 0;
			for(int tmpAlterJunc = 0; tmpAlterJunc < tmpAlterAcceptorSpliceSitePosVecSize; tmpAlterJunc ++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmpAlterJunc];
				if(tmpAlterAcceptorSpliceSitePos == tmpJunc_acceptorStartPos)
					continue;
				int tmpAlterJuncSupNum = juncHash_raw->searchAndReturnAlignInferJuncHashSupNum(tmpJunc_chrNameInt, 
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
			tmpJunc_donerBrranchingRatio = (double)tmpJunc_supportNum / (double)(totalOtherJuncSupNum + tmpJunc_supportNum);
		}
		// check branching ratio at acceptor site, step 1 -- count total branching out junction sup #; step 2 -- total sup # / junc sup #
		vector<int> tmpAlterDonerSpliceSitePosVec;
		SJhashInfo->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpJunc_chrNameInt, tmpJunc_acceptorStartPos, tmpAlterDonerSpliceSitePosVec);
		int tmpAlterDonerSpliceSitePosVecSize = tmpAlterDonerSpliceSitePosVec.size();
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
				int tmpAlterJuncSupNum = juncHash_raw->searchAndReturnAlignInferJuncHashSupNum(tmpJunc_chrNameInt,
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
			tmpJunc_acceptorBranchingRatio = (double)tmpJunc_supportNum / (double)(totalOtherJuncSupNum + tmpJunc_supportNum);
		}

		string tmpJunc_flankString_ref = indexInfo->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
		string tmpJunc_flankString_hap1, tmpJunc_flankString_hap2;
		if(argc == 4)
		{
			tmpJunc_flankString_hap1 = tmpJunc_flankString_ref;
			tmpJunc_flankString_hap2 = tmpJunc_flankString_ref;
		}
		else if(argc == 5)
		{
			tmpJunc_flankString_hap1 = indexInfo_withSNP->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
			tmpJunc_flankString_hap2 = tmpJunc_flankString_hap1;
		}
		else if(argc == 6)
		{
			tmpJunc_flankString_hap1 = indexInfo_hap1->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
			tmpJunc_flankString_hap2 = indexInfo_hap2->returnFlankString(tmpJunc_chrNameInt, tmpJunc_donerEndPos, tmpJunc_acceptorStartPos);
		}
		else
		{
			cout << "error in argc: " << argc << endl;
			exit(1);
		}
		branchingRatioJunc_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t" << tmpJunc_donerEndPos << "\t" << tmpJunc_acceptorStartPos 
			<< "\t" << tmpJunc_supportNum << "\t" << tmpJunc_donerBrranchingRatio << "\t" << tmpJunc_acceptorBranchingRatio 
			<< "\t" << tmpJunc_flankString_ref << "\t" << tmpJunc_flankString_hap1 << "\t" << tmpJunc_flankString_hap2 << endl;;
	}
	branchingRatioJunc_ofs.close();
	delete juncHash_raw;
	delete SJhashInfo;
	delete indexInfo;
	free(chrom);
	chrom_bit_file_ifs.close();
	//parameter_ifs.close();
	delete indexInfo_withSNP;
	delete indexInfo_hap1;
	delete indexInfo_hap2;
	return 0;
}