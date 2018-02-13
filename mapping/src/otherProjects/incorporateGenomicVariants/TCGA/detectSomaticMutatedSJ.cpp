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
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"
#include "../general/SNPhash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolder inputSomaticSNPfile inputSJfile_normal inputSJfile_tumor outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[5];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "inputIndexFolder: " << argv[1] << endl;
	log_ofs << "inputSomaticSNPfile: " << argv[2] << endl;
	log_ofs << "inputSJfile_normal: " << argv[3] << endl;
	log_ofs << "inputSJfile_tumor: " << argv[4] << endl;
	log_ofs << "outputFolder: " << argv[5] << endl;

	cout << "start to initiate indexInfo for both sd and ps genome" << endl;
	log_ofs << "start to initiate indexInfo for both sd and ps genome" << endl;
	string indexFolderPath_normal = argv[1];
	string indexFolderPath_tumor = argv[1];
	indexFolderPath_normal += "/";
	indexFolderPath_tumor += "/";
	string chrom_bit_file_normal = indexFolderPath_normal; chrom_bit_file_normal.append("_chrom");
	string chrom_bit_file_tumor = indexFolderPath_tumor; chrom_bit_file_tumor.append("_chrom"); 
	ifstream chrom_bit_file_ifs_normal(chrom_bit_file_normal.c_str(),ios::binary);
	ifstream chrom_bit_file_ifs_tumor(chrom_bit_file_tumor.c_str(),ios::binary);	
	string indexParameterFileStr_normal = indexFolderPath_normal + "/_parameter";
	string indexParameterFileStr_tumor = indexFolderPath_tumor + "/_parameter";
	ifstream parameter_ifs_normal(indexParameterFileStr_normal.c_str());
	ifstream parameter_ifs_tumor(indexParameterFileStr_tumor.c_str());
	Index_Info* indexInfo_normal = new Index_Info(parameter_ifs_normal);
	Index_Info* indexInfo_tumor = new Index_Info(parameter_ifs_tumor);
	int chromNum_normal = indexInfo_normal->returnChromNum();
	int chromNum_tumor = indexInfo_tumor->returnChromNum();
	if(chromNum_normal != chromNum_tumor){
		cout << "different chrom # for sd and ps genome" << endl;
		exit(1);
	}
	char *chrom_normal; chrom_normal = (char*)malloc((indexInfo_normal->returnIndexSize()) * sizeof(char));
	char *chrom_tumor; chrom_tumor = (char*)malloc((indexInfo_tumor->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs_normal.read((char*)chrom_normal, (indexInfo_normal->returnIndexSize()) * sizeof(char));
	chrom_bit_file_ifs_tumor.read((char*)chrom_tumor, (indexInfo_tumor->returnIndexSize()) * sizeof(char));
	indexInfo_normal->readGenome(chrom_normal);
	indexInfo_tumor->readGenome(chrom_tumor);
	indexInfo_normal->initiate();
	indexInfo_tumor->initiate();
	indexInfo_normal->initiateChrNameIndexArray(1000);
	indexInfo_tumor->initiateChrNameIndexArray(1000);
	cout << "insert SNPs into indexInfo_tumor" << endl;
	log_ofs << "insert SNPs into indexInfo_tumor" << endl;
	string somaticSNPfile = argv[2];
	indexInfo_tumor->insertSNP2chromStr(somaticSNPfile, log_ofs);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;

	cout << "loading somatic SNPs ......" << endl;
	log_ofs << "loading somatic SNPs ......" << endl;
	SNPhash_Info* SNPhashInfo = new SNPhash_Info();
	SNPhashInfo->initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum_normal);
	string inputSomaticSNPfile = argv[2];
	SNPhashInfo->generateSNPhash_formattedSNPfile(inputSomaticSNPfile, indexInfo_normal);

	cout << "loading juncHash_normal ......" << endl;
	log_ofs << "loading juncHash_normal ......" << endl;
	AlignInferJunctionHash_Info* juncHash_normal = new AlignInferJunctionHash_Info();
	juncHash_normal->initiateAlignInferJunctionInfo(chromNum_normal);
	string inputSJfile_normal = argv[3];
	juncHash_normal->insertJuncFromJuncFile_chrNamePos_supportNum(inputSJfile_normal, indexInfo_normal);

	cout << "loading juncHash_tumor ......" << endl;
	log_ofs << "loading juncHash_tumor ......" << endl;
	AlignInferJunctionHash_Info* juncHash_tumor = new AlignInferJunctionHash_Info();
	juncHash_tumor->initiateAlignInferJunctionInfo(chromNum_tumor);
	string inputSJfile_tumor = argv[4];
	juncHash_tumor->insertJuncFromJuncFile_chrNamePos_supportNum(inputSJfile_tumor, indexInfo_tumor);		

	cout << "loading juncHash_merged ......" << endl;
	log_ofs << "loading juncHash_merged ......" << endl;	
	AlignInferJunctionHash_Info* juncHash_merged = new AlignInferJunctionHash_Info();
	juncHash_merged->initiateAlignInferJunctionInfo(chromNum_normal);
	vector<string> inputSJfileVec_merged;
	inputSJfileVec_merged.push_back(inputSJfile_normal);
	inputSJfileVec_merged.push_back(inputSJfile_tumor);
	juncHash_merged->insertJuncFromJuncFileVec_chrNamePosOnly(inputSJfileVec_merged, indexInfo_normal);

	cout << "start to output mutated splice junctions ......" << endl;
	log_ofs << "start to output mutated splice junctions ......" << endl;
	string somaticMutatedSJ_file = outputFolderStr + "/somaticMutatedSJ.txt";
	ofstream somaticMutatedSJ_ofs(somaticMutatedSJ_file.c_str());
	string somaticMutatedSJ_file_canonical2others = outputFolderStr + "/somaticMutatedSJ_canonical2others.txt";
	string somaticMutatedSJ_file_others2canonical = outputFolderStr + "/somaticMutatedSJ_others2canonical.txt";
	string somaticMutatedSJ_file_otherTypes = outputFolderStr + "/somaticMutatedSJ_otherTypes.txt";
	ofstream somaticMutatedSJ_canonical2others_ofs(somaticMutatedSJ_file_canonical2others.c_str());
	ofstream somaticMutatedSJ_others2canonical_ofs(somaticMutatedSJ_file_others2canonical.c_str());
	ofstream somaticMutatedSJ_otherTypes_ofs(somaticMutatedSJ_file_otherTypes.c_str());
	int juncNum_merged = juncHash_merged->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < juncNum_merged; tmp++)
	{
		int tmpChrNameInt = juncHash_merged->returnAlignInferInfo_chrNameInt(tmp); 
		int tmpDonorSite = juncHash_merged->returnAlignInferInfo_donerEndPos(tmp);
		int tmpAcceptorSite = juncHash_merged->returnAlignInferInfo_acceptorStartPos(tmp);
		string tmpFlankString_normal = indexInfo_normal->returnFlankString(
			tmpChrNameInt, tmpDonorSite, tmpAcceptorSite);
		string tmpFlankString_tumor = indexInfo_tumor->returnFlankString(
			tmpChrNameInt, tmpDonorSite, tmpAcceptorSite);
		if(tmpFlankString_normal != tmpFlankString_tumor)
		{
			int tmpSJ_supNum_normal = juncHash_normal->searchAndReturnAlignInferJuncHashSupNum(
				tmpChrNameInt, tmpDonorSite, tmpAcceptorSite);
			int tmpSJ_supNum_tumor = juncHash_tumor->searchAndReturnAlignInferJuncHashSupNum(
				tmpChrNameInt, tmpDonorSite, tmpAcceptorSite);
			somaticMutatedSJ_ofs << indexInfo_normal->returnChrNameStr(tmpChrNameInt) << "\t" 
				<< tmpDonorSite << "\t" << tmpAcceptorSite << "\t" << tmpFlankString_normal 
				<< "\t" << tmpFlankString_tumor << "\t" << tmpSJ_supNum_normal
				<< "\t" << tmpSJ_supNum_tumor << endl;
			if(((tmpFlankString_normal == "GTAG")||(tmpFlankString_normal == "CTAC"))
				&&(tmpFlankString_tumor != "GTAG")&&(tmpFlankString_tumor != "CTAC"))
				somaticMutatedSJ_canonical2others_ofs << indexInfo_normal->returnChrNameStr(tmpChrNameInt) << "\t" 
					<< tmpDonorSite << "\t" << tmpAcceptorSite << "\t" << tmpFlankString_normal 
					<< "\t" << tmpFlankString_tumor << "\t" << tmpSJ_supNum_normal
					<< "\t" << tmpSJ_supNum_tumor << endl;
			else if(((tmpFlankString_tumor == "GTAG")||(tmpFlankString_tumor == "CTAC"))
				&&(tmpFlankString_normal != "GTAG")&&(tmpFlankString_normal != "CTAC"))
				somaticMutatedSJ_others2canonical_ofs << indexInfo_normal->returnChrNameStr(tmpChrNameInt) << "\t" 
					<< tmpDonorSite << "\t" << tmpAcceptorSite << "\t" << tmpFlankString_normal 
					<< "\t" << tmpFlankString_tumor << "\t" << tmpSJ_supNum_normal
					<< "\t" << tmpSJ_supNum_tumor << endl;
			else
				somaticMutatedSJ_otherTypes_ofs << indexInfo_normal->returnChrNameStr(tmpChrNameInt) << "\t" 
					<< tmpDonorSite << "\t" << tmpAcceptorSite << "\t" << tmpFlankString_normal 
					<< "\t" << tmpFlankString_tumor << "\t" << tmpSJ_supNum_normal
					<< "\t" << tmpSJ_supNum_tumor << endl;
		}
	}
	somaticMutatedSJ_ofs.close();
	somaticMutatedSJ_canonical2others_ofs.close();
	somaticMutatedSJ_others2canonical_ofs.close();
	somaticMutatedSJ_otherTypes_ofs.close();
	delete SNPhashInfo;
	delete juncHash_normal;
	delete juncHash_tumor;
	free(chrom_normal);
	free(chrom_tumor);
	delete indexInfo_normal;
	delete indexInfo_tumor;
	return 0;
}