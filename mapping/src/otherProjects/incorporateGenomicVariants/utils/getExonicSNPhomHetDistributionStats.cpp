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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../general/SNPhash_info.h"

using namespace std;

//typedef set<int> 

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputAnnotationPath inputFormattedSNPfile_prefix outputStatsFile" << endl;
		exit(1);
	}
	string inputFormattedSNPfile_prefix = argv[3];
	string outputStatsFile = argv[4];
	cout << "Command: \n" << argv[0] << endl << argv[1] << endl << argv[2] << endl << argv[3] << endl << argv[4] << endl;
	cout << "initiate indexInfo ..." << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;
	////////////////// transcript info loading .... //////////
	string transcript_file_path = argv[2];
	cout << "start to load transcriptInfo " << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type_GAF = "GAF";
	transcriptInfo->extractTranscript(transcript_file_ifs, indexInfo, transcript_type_GAF);
	///////////////////////  load SNPs //////////////////////////////////
	////// het_unique_1, het_unique_2, het_sharePos_diffBase_1, het_sharePos_diffBase_2, hom_shared, 
	SNPhash_Info tmpSNPhashInfo_het_unique_1;
	SNPhash_Info tmpSNPhashInfo_het_unique_2;
	SNPhash_Info tmpSNPhashInfo_het_sharedPos_diffBase_1;
	SNPhash_Info tmpSNPhashInfo_het_sharedPos_diffBase_2;
	SNPhash_Info tmpSNPhashInfo_hom_shared;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo_het_unique_1.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_het_unique_2.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_het_sharedPos_diffBase_1.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_het_sharedPos_diffBase_2.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_hom_shared.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);

	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	string inputFormattedSNPfile_het_unique_1 = inputFormattedSNPfile_prefix + "unique_1_SNP.txt";
	string inputFormattedSNPfile_het_unique_2 = inputFormattedSNPfile_prefix + "unique_2_SNP.txt";
	string inputFormattedSNPfile_het_sharedPos_diffBase_1 = inputFormattedSNPfile_prefix + "sharedPos_difBase_1_SNP.txt";
	string inputFormattedSNPfile_het_sharedPos_diffBase_2 = inputFormattedSNPfile_prefix + "sharedPos_difBase_2_SNP.txt";
	string inputFormattedSNPfile_hom_shared = inputFormattedSNPfile_prefix + "shared_SNP.txt";

	tmpSNPhashInfo_het_unique_1.generateSNPhash_formattedSNPfile(inputFormattedSNPfile_het_unique_1, indexInfo);
	tmpSNPhashInfo_het_unique_2.generateSNPhash_formattedSNPfile(inputFormattedSNPfile_het_unique_2, indexInfo);
	tmpSNPhashInfo_het_sharedPos_diffBase_1.generateSNPhash_formattedSNPfile(inputFormattedSNPfile_het_sharedPos_diffBase_1, indexInfo);
	tmpSNPhashInfo_het_sharedPos_diffBase_2.generateSNPhash_formattedSNPfile(inputFormattedSNPfile_het_sharedPos_diffBase_2, indexInfo);
	tmpSNPhashInfo_hom_shared.generateSNPhash_formattedSNPfile(inputFormattedSNPfile_hom_shared, indexInfo);

	/////////////////////// start to extract SNPs in transcriptome ///////////////////////////
 	vector< set<int> > SNPposSetVec_het_unique_1;
 	vector< set<int> > SNPposSetVec_het_unique_2;
 	vector< set<int> > SNPposSetVec_het_sharedPos_diffBase_1;
 	vector< set<int> > SNPposSetVec_het_sharedPos_diffBase_2;
 	vector< set<int> > SNPposSetVec_hom_shared;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpSNPposSet_het_unique_1;
		set<int> tmpSNPposSet_het_unique_2;
		set<int> tmpSNPposSet_het_sharedPos_diffBase_1;
		set<int> tmpSNPposSet_het_sharedPos_diffBase_2;
		set<int> tmpSNPposSet_hom_shared;
		SNPposSetVec_het_unique_1.push_back(tmpSNPposSet_het_unique_1);
		SNPposSetVec_het_unique_2.push_back(tmpSNPposSet_het_unique_2);
		SNPposSetVec_het_sharedPos_diffBase_1.push_back(tmpSNPposSet_het_sharedPos_diffBase_1);
		SNPposSetVec_het_sharedPos_diffBase_2.push_back(tmpSNPposSet_het_sharedPos_diffBase_2);
		SNPposSetVec_hom_shared.push_back(tmpSNPposSet_hom_shared);
 	}

 	cout << "start to extract SNP from gene ann" << endl;
	int transcriptNum = transcriptInfo->returnTranscriptNum();
	cout << "transcriptNum: " << transcriptNum << endl;
	for(int tmp = 0; tmp < transcriptNum; tmp++)
	{	
		string tmpTranscriptChromName = transcriptInfo->returnTranscriptChromName(tmp);
		int tmpTranscriptChromNameInt = indexInfo->convertStringToInt(tmpTranscriptChromName);
		vector<int> tmpTranscriptExonStartPosVec;
		vector<int> tmpTranscriptExonEndPosVec;
		transcriptInfo->copyExonPos2anotherVec(tmp,	tmpTranscriptExonStartPosVec, tmpTranscriptExonEndPosVec);
		int tmpExonNum = transcriptInfo->returnTranscriptExonNum(tmp);
		for(int tmpExon = 0; tmpExon < tmpExonNum; tmpExon++)
		{
			// generate SNP pos and base in each exon
			vector<int> candiSNPposVecInTmpExon_het_unique_1;
			vector<int> candiSNPposVecInTmpExon_het_unique_2;
			vector<int> candiSNPposVecInTmpExon_het_sharedPos_diffBase_1;
			vector<int> candiSNPposVecInTmpExon_het_sharedPos_diffBase_2;
			vector<int> candiSNPposVecInTmpExon_hom_shared;
			int tmpExonStartPos = tmpTranscriptExonStartPosVec[tmpExon];
			int tmpExonEndPos = tmpTranscriptExonEndPosVec[tmpExon];
			tmpSNPhashInfo_het_unique_1.returnSNPposVecWithinRegion(tmpTranscriptChromNameInt, 
				tmpExonStartPos, tmpExonEndPos, candiSNPposVecInTmpExon_het_unique_1);
			tmpSNPhashInfo_het_unique_2.returnSNPposVecWithinRegion(tmpTranscriptChromNameInt, 
				tmpExonStartPos, tmpExonEndPos, candiSNPposVecInTmpExon_het_unique_2);
			tmpSNPhashInfo_het_sharedPos_diffBase_1.returnSNPposVecWithinRegion(tmpTranscriptChromNameInt, 
				tmpExonStartPos, tmpExonEndPos, candiSNPposVecInTmpExon_het_sharedPos_diffBase_1);
			tmpSNPhashInfo_het_sharedPos_diffBase_2.returnSNPposVecWithinRegion(tmpTranscriptChromNameInt, 
				tmpExonStartPos, tmpExonEndPos, candiSNPposVecInTmpExon_het_sharedPos_diffBase_2);
			tmpSNPhashInfo_hom_shared.returnSNPposVecWithinRegion(tmpTranscriptChromNameInt, 
				tmpExonStartPos, tmpExonEndPos, candiSNPposVecInTmpExon_hom_shared);															
			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon_het_unique_1.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon_het_unique_1[tmpSNP];
				SNPposSetVec_het_unique_1[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}
			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon_het_unique_2.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon_het_unique_2[tmpSNP];
				SNPposSetVec_het_unique_2[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}

			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon_het_sharedPos_diffBase_1.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon_het_sharedPos_diffBase_1[tmpSNP];
				SNPposSetVec_het_sharedPos_diffBase_1[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}

			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon_het_sharedPos_diffBase_2.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon_het_sharedPos_diffBase_2[tmpSNP];
				SNPposSetVec_het_sharedPos_diffBase_2[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}

			for(int tmpSNP = 0; tmpSNP < candiSNPposVecInTmpExon_hom_shared.size(); tmpSNP++)
			{
				int tmpSNPpos = candiSNPposVecInTmpExon_hom_shared[tmpSNP];
				SNPposSetVec_hom_shared[tmpTranscriptChromNameInt].insert(tmpSNPpos);
			}										
		}	
	}

	cout << "start to get SNP_num in each category" << endl;
	ofstream stats_ofs(outputStatsFile.c_str());
	int exonic_SNPposNum_het_unique_1 = 0, 
		exonic_SNPposNum_het_unique_2 = 0,
		exonic_SNPposNum_het_sharedPos_diffBase_1 = 0, 
		exonic_SNPposNum_het_sharedPos_diffBase_2 = 0,
		exonic_SNPposNum_hom_shared = 0;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		exonic_SNPposNum_het_unique_1 += SNPposSetVec_het_unique_1[tmp].size();
		exonic_SNPposNum_het_unique_2 += SNPposSetVec_het_unique_2[tmp].size();
		exonic_SNPposNum_het_sharedPos_diffBase_1 += SNPposSetVec_het_sharedPos_diffBase_1[tmp].size(); 
		exonic_SNPposNum_het_sharedPos_diffBase_2 += SNPposSetVec_het_sharedPos_diffBase_2[tmp].size();
		exonic_SNPposNum_hom_shared += SNPposSetVec_hom_shared[tmp].size();
	}
	int SNPposNum_het_unique_1 = tmpSNPhashInfo_het_unique_1.returnSNPnum();
	int SNPposNum_het_unique_2 = tmpSNPhashInfo_het_unique_2.returnSNPnum();
	int	SNPposNum_het_sharedPos_diffBase_1 = tmpSNPhashInfo_het_sharedPos_diffBase_1.returnSNPnum(); 
	int SNPposNum_het_sharedPos_diffBase_2 = tmpSNPhashInfo_het_sharedPos_diffBase_2.returnSNPnum();
	int SNPposNum_hom_shared = tmpSNPhashInfo_hom_shared.returnSNPnum();

	stats_ofs << "SNPposNum_het_unique_1: " << endl << SNPposNum_het_unique_1 << endl;
	stats_ofs << "SNPposNum_het_unique_2: " << endl << SNPposNum_het_unique_2 << endl;
	stats_ofs << "SNPposNum_het_sharedPos_diffBase_1: " << endl << SNPposNum_het_sharedPos_diffBase_1 << endl;
	stats_ofs << "SNPposNum_het_sharedPos_diffBase_2: " << endl << SNPposNum_het_sharedPos_diffBase_2 << endl;
	stats_ofs << "SNPposNum_hom_shared: " << endl << SNPposNum_hom_shared << endl;
	stats_ofs << endl;
	stats_ofs << "exonic_SNPposNum_het_unique_1: " << endl << exonic_SNPposNum_het_unique_1 << endl;
	stats_ofs << "exonic_SNPposNum_het_unique_2: " << endl << exonic_SNPposNum_het_unique_2 << endl;
	stats_ofs << "exonic_SNPposNum_het_sharedPos_diffBase_1: " << endl << exonic_SNPposNum_het_sharedPos_diffBase_1 << endl;
	stats_ofs << "exonic_SNPposNum_het_sharedPos_diffBase_2: " << endl << exonic_SNPposNum_het_sharedPos_diffBase_2 << endl;
	stats_ofs << "exonic_SNPposNum_hom_shared: " << endl << exonic_SNPposNum_hom_shared << endl;	

	stats_ofs.close();
	cout << "all jobs done ..." << endl;
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}