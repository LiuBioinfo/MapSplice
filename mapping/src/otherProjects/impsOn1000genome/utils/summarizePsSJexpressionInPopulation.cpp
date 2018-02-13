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
#include "../general/flkStrChangedSJ_info.h"

using namespace std;

bool isCanSJ(string& tmpFlkStr)	
{
	if((tmpFlkStr == "GTAG")||(tmpFlkStr == "CTAC")
		||(tmpFlkStr == "ATAC")||(tmpFlkStr == "GTAT")
		||(tmpFlkStr == "GCAG")||(tmpFlkStr == "CTGC"))
		return true;
	else
		return false;
}

int canSJhapNum(string& tmpFlkStr_hap1, string& tmpFlkStr_hap2)
{
	int tmpCanSJhapNum = 0;
	if(isCanSJ(tmpFlkStr_hap1))
		tmpCanSJhapNum ++;
	if(isCanSJ(tmpFlkStr_hap2))
		tmpCanSJhapNum ++;	
	return tmpCanSJhapNum;
}

void summarizeCanSJhapSampleSupNumPair(vector< pair<int,int> >& tmpCanSJhapNumAndSupNumPairVec, 
	int& tmpCanSJ_nonHap_sampleNum, int& tmpCanSJ_nonHap_supNum, 
	int& tmpCanSJ_existHap_sampleNum, int& tmpCanSJ_existHap_supNum, 
	int& tmpCanSJ_singleHap_sampleNum, int& tmpCanSJ_singleHap_supNum,
	int& tmpCanSJ_bothHap_sampleNum, int& tmpCanSJ_bothHap_supNum)
{
	tmpCanSJ_nonHap_sampleNum = 0; 
	tmpCanSJ_nonHap_supNum = 0;
	tmpCanSJ_existHap_sampleNum = 0;
	tmpCanSJ_existHap_supNum = 0;
	tmpCanSJ_singleHap_sampleNum = 0;
	tmpCanSJ_singleHap_supNum = 0;
	tmpCanSJ_bothHap_sampleNum = 0;
	tmpCanSJ_bothHap_supNum	= 0;
	for(int tmp = 0; tmp < tmpCanSJhapNumAndSupNumPairVec.size(); tmp++)
	{
		int tmpCanHapNum = tmpCanSJhapNumAndSupNumPairVec[tmp].first;
		int tmpSupNum = tmpCanSJhapNumAndSupNumPairVec[tmp].second;
		if(tmpCanHapNum == 0)
		{
			tmpCanSJ_nonHap_sampleNum ++;
			tmpCanSJ_nonHap_supNum += tmpSupNum;
		}
		else if(tmpCanHapNum == 1)
		{
			tmpCanSJ_singleHap_sampleNum ++;
			tmpCanSJ_singleHap_supNum += tmpSupNum;
			tmpCanSJ_existHap_sampleNum ++;
			tmpCanSJ_existHap_supNum += tmpSupNum;
		}
		else if(tmpCanHapNum == 2)
		{
			tmpCanSJ_bothHap_sampleNum ++;
			tmpCanSJ_bothHap_supNum += tmpSupNum;
			tmpCanSJ_existHap_sampleNum ++;
			tmpCanSJ_existHap_supNum += tmpSupNum;
		}
		else
		{
			cout << "invalid tmpCanHapNum: " << tmpCanHapNum << endl;
			exit(1);
		}
	}
}

bool checkFlkStrPatternChangeOrNot(int label_ref, vector<int> labelVec_hap)
{
	if(label_ref == 0)
	{
		for(int tmp = 0; tmp < labelVec_hap.size(); tmp++)
		{
			if(labelVec_hap[tmp] > 0)
				return true;
		}
		return false;
	}
	else if(label_ref == 1)
	{
		for(int tmp = 0; tmp < labelVec_hap.size(); tmp++)
		{
			if(labelVec_hap[tmp] == 0)
				return true;
		}
		return false;		
	}
	else
	{
		cout << "invalid label_ref: " << label_ref << endl;
		exit(1);
	}
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 index_info" << endl;
		cout << "#2 SNPfileList_hap1" << endl;
		cout << "#3 SNPfileList_hap2" << endl;
		cout << "#4 changedSJlist" << endl;
		cout << "#5 outputFilePrefix" << endl;
		exit(1);
	}
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
	//log_ofs << "end of initiating indexInfo" << endl;
	free(chrom);
	parameter_ifs.close();

	cout << "start to read SNPfileList_hap1/2 changedSJfileVec..." << endl;
	string SNPfileList_hap1 = argv[2];
	string SNPfileList_hap2 = argv[3];
	string changedSJfileList = argv[4];
	vector<string> SNPfileVec_hap1;
	vector<string> SNPfileVec_hap2;
	vector<string> changedSJfileVec;
	ifstream SNPlist_hap1_ifs(SNPfileList_hap1.c_str());
	ifstream SNPlist_hap2_ifs(SNPfileList_hap2.c_str());
	ifstream changedSJfileList_ifs(changedSJfileList.c_str());
	while((!SNPlist_hap1_ifs.eof())&&(!SNPlist_hap2_ifs.eof())&&(!changedSJfileList_ifs.eof()))
	{
		string tmpStr_SNP_hap1, tmpStr_SNP_hap2, tmpStr_SJ;
		getline(SNPlist_hap1_ifs, tmpStr_SNP_hap1);
		getline(SNPlist_hap2_ifs, tmpStr_SNP_hap2);
		getline(changedSJfileList_ifs, tmpStr_SJ);
		if((tmpStr_SNP_hap1 == "")||(tmpStr_SNP_hap2 == "")||(tmpStr_SJ == ""))
			break;
		SNPfileVec_hap1.push_back(tmpStr_SNP_hap1);
		SNPfileVec_hap2.push_back(tmpStr_SNP_hap2);
		changedSJfileVec.push_back(tmpStr_SJ);
	}
	SNPlist_hap1_ifs.close();
	SNPlist_hap2_ifs.close();
	changedSJfileList_ifs.close();

	int fileNum = changedSJfileVec.size();
	cout << "fileNum: " << fileNum << endl;
	vector<FlkStrChangedSJ_Hash_Info> flkStrChangedSJhashInfoVec;
	FlkStrChangedSJ_Hash_Info flkStrChangedSJhashInfo_merged;
	flkStrChangedSJhashInfo_merged.initiate_chrNum(chromNum);
	cout<< "start to load each flkStrChangedSJ_file" << endl;
	for(int tmp = 0; tmp < fileNum; tmp ++)
	{
		flkStrChangedSJhashInfo_merged.loadFlkStrChangedSJfile(changedSJfileVec[tmp], indexInfo);
		FlkStrChangedSJ_Hash_Info tmpFlkStrChangedSJhashInfo;
		tmpFlkStrChangedSJhashInfo.initiate_chrNum(chromNum);
		tmpFlkStrChangedSJhashInfo.loadFlkStrChangedSJfile(changedSJfileVec[tmp], indexInfo);
		flkStrChangedSJhashInfoVec.push_back(tmpFlkStrChangedSJhashInfo);
	}

	int totalFlkStrChangedSJnum = flkStrChangedSJhashInfo_merged.return_flkStrChangedSJnum();
	cout << "totalFlkStrChangedSJnum: " << totalFlkStrChangedSJnum << endl;
	
	vector<int> chrNameIntVec;           // juncVec <chrNameInt>
	vector<string> chrNameStrVec;        // juncVec <chrNameStr>
	vector< pair<int, int> > posPairVec; // juncVec <startPos, endPos>
	vector< vector<int> > supNumVecVec;  // juncVec < sampleVec <supNum> >
	vector<int> canJuncHapNumVec_ref;    // juncVec < sampleVec <0/1> > -- 0:noncanSJ; 1:canSJ	
	for(int tmp = 0; tmp < totalFlkStrChangedSJnum; tmp++)
	{
		int tmpChrNameInt = flkStrChangedSJhashInfo_merged.return_chrNameInt_vecIndex(tmp);
		string tmpChrNameStr = indexInfo->returnChrNameStr(tmpChrNameInt);
		int tmpDonerEndPos = flkStrChangedSJhashInfo_merged.return_donerEndPos_vecIndex(tmp);
		int tmpAcceptorStartPos = flkStrChangedSJhashInfo_merged.return_acceptorStartPos_vecIndex(tmp);
		string tmpFlkStr_ref = flkStrChangedSJhashInfo_merged.return_flkStr_ref_vecIndex(tmp);
		bool tmpFlkStr_ref_canOrNot_bool = isCanSJ(tmpFlkStr_ref);
		if(tmpFlkStr_ref_canOrNot_bool)
			canJuncHapNumVec_ref.push_back(1);
		else
			canJuncHapNumVec_ref.push_back(0);
		chrNameIntVec.push_back(tmpChrNameInt);
		chrNameStrVec.push_back(tmpChrNameStr);
		posPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
		vector<int> tmpSupNumVec;
		for(int tmp2 = 0; tmp2 < fileNum; tmp2++)
		{
			int tmpSupNum = flkStrChangedSJhashInfoVec[tmp2].search_and_return_supNum( 
				tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
			tmpSupNumVec.push_back(tmpSupNum);
		}
		supNumVecVec.push_back(tmpSupNumVec);
	}

	cout << "start to initiate canJuncHapNumVecVec_hap" << endl;
	vector< vector<int> > canJuncHapNumVecVec_hap; // juncVec < sampleVec <canJuncHapNum(0/1/2)> >
	for(int tmpJunc = 0; tmpJunc < totalFlkStrChangedSJnum; tmpJunc ++)
	{
		vector<int> tmpCanJuncHapNumVec;
		for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
			tmpCanJuncHapNumVec.push_back(0);
		canJuncHapNumVecVec_hap.push_back(tmpCanJuncHapNumVec);
	}	

	cout << "start to fill canJuncHapNumVecVec_hap" << endl;
	for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
	{
		// load personal genomes
		cout << "tmpSample: " << tmpSample << endl;
		string tmp_SNP_hap1_file = SNPfileVec_hap1[tmpSample];
		string tmp_SNP_hap2_file = SNPfileVec_hap2[tmpSample];
		Index_Info* indexInfo_hap1 = new Index_Info();
		Index_Info* indexInfo_hap2 = new Index_Info();
		indexInfo_hap1->cpIndex(indexInfo);
		indexInfo_hap2->cpIndex(indexInfo);
		indexInfo_hap1->insertSNP2chromStr(tmp_SNP_hap1_file);
		indexInfo_hap2->insertSNP2chromStr(tmp_SNP_hap2_file);
		for(int tmpJunc = 0; tmpJunc < totalFlkStrChangedSJnum; tmpJunc ++)
		{
			int tmpChrNameInt = chrNameIntVec[tmpJunc];
			int tmpDonerEndPos = posPairVec[tmpJunc].first;
			int tmpAcceptorStartPos = posPairVec[tmpJunc].second;
			string tmpFlkStr_hap1 = indexInfo_hap1->returnFlankString(
				tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
			string tmpFlkStr_hap2 = indexInfo_hap2->returnFlankString(
				tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
			int tmpCanSJhapNum = canSJhapNum(tmpFlkStr_hap1, tmpFlkStr_hap2);
			(canJuncHapNumVecVec_hap[tmpJunc])[tmpSample] = tmpCanSJhapNum;
		}
		// delete indexInfo for both hap
		delete indexInfo_hap1;
		delete indexInfo_hap2;
	}

	cout << "start to print each junc " << endl;
	string outputFilePrefix = argv[5]; 
	string outputFile_changed = outputFilePrefix + "patternChanged.txt";
	string outputFile_kept = outputFilePrefix + "patternKept.txt";
	ofstream changed_ofs(outputFile_changed.c_str());
	ofstream kept_ofs(outputFile_kept.c_str());
	changed_ofs << "chrName\tpos_1\tpos_2\tflkStr_ref";
	for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
		changed_ofs << "\tflkStr_hap_" << tmpSample + 1 << "\tsupNum_" << tmpSample + 1;
	changed_ofs << endl; 
	kept_ofs << "chrName\tpos_1\tpos_2\tflkStr_ref";
	for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
		kept_ofs << "\tflkStr_hap_" << tmpSample + 1 << "\tsupNum_" << tmpSample + 1;
	kept_ofs << endl;
	for(int tmpJunc = 0; tmpJunc < totalFlkStrChangedSJnum; tmpJunc ++)
	{
		bool SJflkStrPattern_changed_bool = checkFlkStrPatternChangeOrNot(
			canJuncHapNumVec_ref[tmpJunc], canJuncHapNumVecVec_hap[tmpJunc]);

		if(SJflkStrPattern_changed_bool)
		{
			changed_ofs << chrNameStrVec[tmpJunc] << "\t" << posPairVec[tmpJunc].first 
				<< "\t" << posPairVec[tmpJunc].second << "\t" << canJuncHapNumVec_ref[tmpJunc];
			for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
				changed_ofs << "\t" << (canJuncHapNumVecVec_hap[tmpJunc])[tmpSample]
					<< "\t" << (supNumVecVec[tmpJunc])[tmpSample];
			changed_ofs << endl;
		}
		else
		{
			kept_ofs << chrNameStrVec[tmpJunc] << "\t" << posPairVec[tmpJunc].first 
				<< "\t" << posPairVec[tmpJunc].second << "\t" << canJuncHapNumVec_ref[tmpJunc];
			for(int tmpSample = 0; tmpSample < fileNum; tmpSample ++)
				kept_ofs << "\t" << (canJuncHapNumVecVec_hap[tmpJunc])[tmpSample]
					<< "\t" << (supNumVecVec[tmpJunc])[tmpSample];
			kept_ofs << endl;
		}
	}
	changed_ofs.close();
	kept_ofs.close();
	delete indexInfo;
	return 0;
}