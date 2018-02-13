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

using namespace std;

void extractSJchrNamePosSupNumFromStr(string& tmpSJstr, string& tmpSJ_chrName, 
	int& tmpSJ_startPos, int& tmpSJ_endPos, int& tmpSJ_supNum, string& tmpOtherStr)
{
	int tabLoc_1 = tmpSJstr.find("\t");
	int tabLoc_2 = tmpSJstr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSJstr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpSJstr.find("\t", tabLoc_3 + 1);
	if(tabLoc_4 == string::npos)
	{
		cout << "incorrect SJ format" << endl;
		exit(1);
	}
	int tabLoc_5 = tmpSJstr.find("\t", tabLoc_4 + 1);
	string tmpSupNumStr;
	if(tabLoc_5 == string::npos)
	{	
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1);
		tmpOtherStr = "";
	}
	else
	{	
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);	
		tmpOtherStr = tmpSJstr.substr(tabLoc_5 + 1);
	}
	tmpSJ_chrName = tmpSJstr.substr(0,tabLoc_1);
	string tmpStartPosStr = tmpSJstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEndPosStr = tmpSJstr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);	
	tmpSJ_startPos = atoi(tmpStartPosStr.c_str());
	tmpSJ_endPos = atoi(tmpEndPosStr.c_str());
	tmpSJ_supNum = atoi(tmpSupNumStr.c_str());
}

bool flkStr_can_or_not(string& tmpFlkStr)
{
	if((tmpFlkStr == "GTAG")||(tmpFlkStr == "CTAC"))
		return true;
	else
		return false;
}

bool flkStr_sem_or_not(string& tmpFlkStr)
{
	if((tmpFlkStr == "ATAC")||(tmpFlkStr == "GTAT")||(tmpFlkStr == "CTGC")||(tmpFlkStr == "GCAG"))
		return true;
	else
		return false;
}

bool flkStr_non_or_not(string& tmpFlkStr)
{
	if((tmpFlkStr == "GTAG")||(tmpFlkStr == "CTAC")||(tmpFlkStr == "ATAC")||
		(tmpFlkStr == "GTAT")||(tmpFlkStr == "CTGC")||(tmpFlkStr == "GCAG"))
		return false;
	else
		return true;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexInfoDir" << endl;
		cout << "#2 SNP_hap1" << endl;
		cout << "#3 SNP_hap2" << endl;
		cout << "#4 inputSJfile" << endl;
		cout << "#5 outputDir" << endl;
		exit(1);
	}

	int SJ_size_min = 50;
	int SJ_size_max = 300000;

	string indexInfoDir = argv[1];
	string SNP_hap1 = argv[2];
	string SNP_hap2 = argv[3];
	string inputSJfile = argv[4];
	string outputFolderStr = argv[5];
	cout << "creating results folder ...." << endl;	
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());
	string stats_path = outputFolderStr + "stats.txt";
	ofstream stats_ofs(stats_path.c_str());

	cout << "start to initiate indexInfo for both sd and ps genome" << endl;
	log_ofs << "start to initiate indexInfo for both sd and ps genome" << endl;
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
	log_ofs << "end of initiating indexInfo" << endl;
	free(chrom);

	Index_Info* indexInfo_hap1 = new Index_Info();
	Index_Info* indexInfo_hap2 = new Index_Info();
	indexInfo_hap1->cpIndex(indexInfo);
	indexInfo_hap2->cpIndex(indexInfo);
	indexInfo_hap1->insertSNP2chromStr(SNP_hap1, log_ofs);
	indexInfo_hap2->insertSNP2chromStr(SNP_hap2, log_ofs); 
	
	// not changed
	// changed
	// hom:
	// 1. can->can; 2. can->sem; 3. can->non	
	// 4. sem->can; 5. sem->sem; 6. sem->non
	// 7. non->can; 8. non->sem; 9. non->non	

	// het (from can):
	// 1. can->can/can; 2. can->can/sem; 3. can->can/non
	// 4. can->sem/sem; 5. can->sem/non; 6. can->non/non
	// het (from sem):
	// 1. sem->can/can; 2. sem->can/sem; 3. sem->can/non
	// 4. sem->sem/sem; 5. sem->sem/non; 6. sem->non/non
	// het (from non):
	// 1. non->can/can; 2. non->can/sem; 3. non->can/non
	// 4. non->sem/sem; 5. non->sem/non; 6. non->non/non

	string output_all = outputFolderStr + "all.junc";
	string output_kept = outputFolderStr + "kept.junc";
	string output_changed = outputFolderStr + "changed.junc";
	// string output_changed_hom = outputFolderStr + "changed_hom.junc";
	// string output_changed_het = outputFolderStr + "changed_het.junc";
	// string output_can_2_can = outputFolderStr + "can_2_can.junc";
	// string output_can_2_sem = outputFolderStr + "can_2_sem.junc";
	// string output_can_2_non = outputFolderStr + "can_2_non.junc";
	// string output_sem_2_can = outputFolderStr + "sem_2_can.junc";
	// string output_sem_2_sem = outputFolderStr + "sem_2_sem.junc";
	// string output_sem_2_non = outputFolderStr + "sem_2_non.junc";
	// string output_non_2_can = outputFolderStr + "non_2_can.junc";
	// string output_non_2_sem = outputFolderStr + "non_2_sem.junc";
	// string output_non_2_non = outputFolderStr + "non_2_non.junc";	
	// string output_can_2_can_can = outputFolderStr + "can_2_can_can.junc";
	// string output_can_2_can_sem = outputFolderStr + "can_2_can_sem.junc";
	// string output_can_2_can_non = outputFolderStr + "can_2_can_non.junc";
	// string output_can_2_sem_sem = outputFolderStr + "can_2_sem_sem.junc";
	// string output_can_2_sem_non = outputFolderStr + "can_2_sem_non.junc";
	// string output_can_2_non_non = outputFolderStr + "can_2_non_non.junc";
	// string output_sem_2_can_can = outputFolderStr + "sem_2_can_can.junc";
	// string output_sem_2_can_sem = outputFolderStr + "sem_2_can_sem.junc";
	// string output_sem_2_can_non = outputFolderStr + "sem_2_can_non.junc";
	// string output_sem_2_sem_sem = outputFolderStr + "sem_2_sem_sem.junc";
	// string output_sem_2_sem_non = outputFolderStr + "sem_2_sem_non.junc";
	// string output_sem_2_non_non = outputFolderStr + "sem_2_non_non.junc";
	// string output_non_2_can_can = outputFolderStr + "non_2_can_can.junc";
	// string output_non_2_can_sem = outputFolderStr + "non_2_can_sem.junc";
	// string output_non_2_can_non = outputFolderStr + "non_2_can_non.junc";
	// string output_non_2_sem_sem = outputFolderStr + "non_2_sem_sem.junc";
	// string output_non_2_sem_non = outputFolderStr + "non_2_sem_non.junc";
	// string output_non_2_non_non = outputFolderStr + "non_2_non_non.junc";	

	ofstream all_ofs(output_all.c_str());
	ofstream kept_ofs(output_kept.c_str());
	ofstream changed_ofs(output_changed.c_str());
	// ofstream changed_hom_ofs(output_changed_hom.c_str());
	// ofstream changed_het_ofs(output_changed_het.c_str());
	// ofstream can_2_can_ofs(output_can_2_can.c_str());
	// ofstream can_2_sem_ofs(output_can_2_sem.c_str());
	// ofstream can_2_non_ofs(output_can_2_non.c_str());
	// ofstream sem_2_can_ofs(output_sem_2_can.c_str());
	// ofstream sem_2_sem_ofs(output_sem_2_sem.c_str());
	// ofstream sem_2_non_ofs(output_sem_2_non.c_str());
	// ofstream non_2_can_ofs(output_non_2_can.c_str());
	// ofstream non_2_sem_ofs(output_non_2_sem.c_str());
	// ofstream non_2_non_ofs(output_non_2_non.c_str());	
	// ofstream can_2_can_can_ofs(output_can_2_can_can.c_str());
	// ofstream can_2_can_sem_ofs(output_can_2_can_sem.c_str());
	// ofstream can_2_can_non_ofs(output_can_2_can_non.c_str());
	// ofstream can_2_sem_sem_ofs(output_can_2_sem_sem.c_str());
	// ofstream can_2_sem_non_ofs(output_can_2_sem_non.c_str());
	// ofstream can_2_non_non_ofs(output_can_2_non_non.c_str());
	// ofstream sem_2_can_can_ofs(output_sem_2_can_can.c_str());
	// ofstream sem_2_can_sem_ofs(output_sem_2_can_sem.c_str());
	// ofstream sem_2_can_non_ofs(output_sem_2_can_non.c_str());
	// ofstream sem_2_sem_sem_ofs(output_sem_2_sem_sem.c_str());
	// ofstream sem_2_sem_non_ofs(output_sem_2_sem_non.c_str());
	// ofstream sem_2_non_non_ofs(output_sem_2_non_non.c_str());
	// ofstream non_2_can_can_ofs(output_non_2_can_can.c_str());
	// ofstream non_2_can_sem_ofs(output_non_2_can_sem.c_str());
	// ofstream non_2_can_non_ofs(output_non_2_can_non.c_str());
	// ofstream non_2_sem_sem_ofs(output_non_2_sem_sem.c_str());
	// ofstream non_2_sem_non_ofs(output_non_2_sem_non.c_str());
	// ofstream non_2_non_non_ofs(output_non_2_non_non.c_str());		


	int SJ_all = 0, SJ_can = 0, SJ_sem = 0, SJ_non = 0;
	int SJ_flkStrKept = 0, SJ_flkStrChanged = 0, SJ_flkStrChanged_hom = 0, SJ_flkStrChanged_het = 0;
	int SJ_flkStrChanged_can2can = 0, SJ_flkStrChanged_can2sem = 0, SJ_flkStrChanged_can2non = 0,
		SJ_flkStrChanged_sem2can = 0, SJ_flkStrChanged_sem2sem = 0, SJ_flkStrChanged_sem2non = 0,
		SJ_flkStrChanged_non2can = 0, SJ_flkStrChanged_non2sem = 0, SJ_flkStrChanged_non2non = 0;
	int SJ_flkStrChanged_can2canORcan = 0, SJ_flkStrChanged_can2canORsem = 0, SJ_flkStrChanged_can2canORnon = 0,
		SJ_flkStrChanged_can2semORsem = 0, SJ_flkStrChanged_can2semORnon = 0, SJ_flkStrChanged_can2nonORnon = 0;
	int SJ_flkStrChanged_sem2canORcan = 0, SJ_flkStrChanged_sem2canORsem = 0, SJ_flkStrChanged_sem2canORnon = 0,
		SJ_flkStrChanged_sem2semORsem = 0, SJ_flkStrChanged_sem2semORnon = 0, SJ_flkStrChanged_sem2nonORnon = 0;
	int SJ_flkStrChanged_non2canORcan = 0, SJ_flkStrChanged_non2canORsem = 0, SJ_flkStrChanged_non2canORnon = 0,
		SJ_flkStrChanged_non2semORsem = 0, SJ_flkStrChanged_non2semORnon = 0, SJ_flkStrChanged_non2nonORnon = 0;

	int SJ_all_min3 = 0, SJ_can_min3 = 0, SJ_sem_min3 = 0, SJ_non_min3 = 0;
	int SJ_flkStrKept_min3 = 0, SJ_flkStrChanged_min3 = 0, SJ_flkStrChanged_hom_min3 = 0, SJ_flkStrChanged_het_min3 = 0;
	int SJ_flkStrChanged_can2can_min3 = 0, SJ_flkStrChanged_can2sem_min3 = 0, SJ_flkStrChanged_can2non_min3 = 0,
		SJ_flkStrChanged_sem2can_min3 = 0, SJ_flkStrChanged_sem2sem_min3 = 0, SJ_flkStrChanged_sem2non_min3 = 0,
		SJ_flkStrChanged_non2can_min3 = 0, SJ_flkStrChanged_non2sem_min3 = 0, SJ_flkStrChanged_non2non_min3 = 0;
	int SJ_flkStrChanged_can2canORcan_min3 = 0, SJ_flkStrChanged_can2canORsem_min3 = 0, SJ_flkStrChanged_can2canORnon_min3 = 0,
		SJ_flkStrChanged_can2semORsem_min3 = 0, SJ_flkStrChanged_can2semORnon_min3 = 0, SJ_flkStrChanged_can2nonORnon_min3 = 0;
	int SJ_flkStrChanged_sem2canORcan_min3 = 0, SJ_flkStrChanged_sem2canORsem_min3 = 0, SJ_flkStrChanged_sem2canORnon_min3 = 0,
		SJ_flkStrChanged_sem2semORsem_min3 = 0, SJ_flkStrChanged_sem2semORnon_min3 = 0, SJ_flkStrChanged_sem2nonORnon_min3 = 0;
	int SJ_flkStrChanged_non2canORcan_min3 = 0, SJ_flkStrChanged_non2canORsem_min3 = 0, SJ_flkStrChanged_non2canORnon_min3 = 0,
		SJ_flkStrChanged_non2semORsem_min3 = 0, SJ_flkStrChanged_non2semORnon_min3 = 0, SJ_flkStrChanged_non2nonORnon_min3 = 0;

	int SJ_all_min5 = 0, SJ_can_min5 = 0, SJ_sem_min5 = 0, SJ_non_min5 = 0;
	int SJ_flkStrKept_min5 = 0, SJ_flkStrChanged_min5 = 0, SJ_flkStrChanged_hom_min5 = 0, SJ_flkStrChanged_het_min5 = 0;
	int SJ_flkStrChanged_can2can_min5 = 0, SJ_flkStrChanged_can2sem_min5 = 0, SJ_flkStrChanged_can2non_min5 = 0,
		SJ_flkStrChanged_sem2can_min5 = 0, SJ_flkStrChanged_sem2sem_min5 = 0, SJ_flkStrChanged_sem2non_min5 = 0,
		SJ_flkStrChanged_non2can_min5 = 0, SJ_flkStrChanged_non2sem_min5 = 0, SJ_flkStrChanged_non2non_min5 = 0;
	int SJ_flkStrChanged_can2canORcan_min5 = 0, SJ_flkStrChanged_can2canORsem_min5 = 0, SJ_flkStrChanged_can2canORnon_min5 = 0,
		SJ_flkStrChanged_can2semORsem_min5 = 0, SJ_flkStrChanged_can2semORnon_min5 = 0, SJ_flkStrChanged_can2nonORnon_min5 = 0;
	int SJ_flkStrChanged_sem2canORcan_min5 = 0, SJ_flkStrChanged_sem2canORsem_min5 = 0, SJ_flkStrChanged_sem2canORnon_min5 = 0,
		SJ_flkStrChanged_sem2semORsem_min5 = 0, SJ_flkStrChanged_sem2semORnon_min5 = 0, SJ_flkStrChanged_sem2nonORnon_min5 = 0;
	int SJ_flkStrChanged_non2canORcan_min5 = 0, SJ_flkStrChanged_non2canORsem_min5 = 0, SJ_flkStrChanged_non2canORnon_min5 = 0,
		SJ_flkStrChanged_non2semORsem_min5 = 0, SJ_flkStrChanged_non2semORnon_min5 = 0, SJ_flkStrChanged_non2nonORnon_min5 = 0;		

	log_ofs << "start to scanning each junction and its flankstrings regarding ref, hap1 and hap2" << endl;
	ifstream SJ_ifs(inputSJfile.c_str());
	while(!SJ_ifs.eof())
	{
		string tmpSJstr;
		getline(SJ_ifs, tmpSJstr);
		if(tmpSJstr == "")
			break;
		string tmpSJ_chrName;
		int tmpSJ_startPos;
		int tmpSJ_endPos;
		int tmpSJ_supNum;
		string tmpOtherStr;
		extractSJchrNamePosSupNumFromStr(tmpSJstr, tmpSJ_chrName, tmpSJ_startPos, 
			tmpSJ_endPos, tmpSJ_supNum, tmpOtherStr);
		int tmpSJ_size = tmpSJ_endPos - tmpSJ_startPos - 1;
		if((tmpSJ_size < SJ_size_min)||(tmpSJ_size > SJ_size_max))
			continue;		
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrName);
		if(tmpSJ_chrNameInt < 0)
		{
			cout << "SJ chrNameInt < 0 !" << endl;
			exit(1);
		}
		string tmpSJ_flankString_raw = indexInfo->returnFlankString(
			tmpSJ_chrNameInt, tmpSJ_startPos, tmpSJ_endPos);
		string tmpSJ_flankString_hap1 = indexInfo_hap1->returnFlankString(
			tmpSJ_chrNameInt, tmpSJ_startPos, tmpSJ_endPos);
		string tmpSJ_flankString_hap2 = indexInfo_hap2->returnFlankString(
			tmpSJ_chrNameInt, tmpSJ_startPos, tmpSJ_endPos);
		string tmpSJinfoStr = tmpSJ_chrName + "\t" + int_to_str(tmpSJ_startPos) + "\t" + int_to_str(tmpSJ_endPos) 
			+ "\t" + int_to_str(tmpSJ_supNum) + "\t" + tmpSJ_flankString_raw + "\t" + tmpSJ_flankString_hap1 + "\t"
			+ tmpSJ_flankString_hap2 + "\t" + tmpOtherStr;
		all_ofs << tmpSJinfoStr << endl;
		SJ_all ++;
		if(tmpSJ_supNum >= 3)
			SJ_all_min3 ++;
		if(tmpSJ_supNum >= 5)
			SJ_all_min5 ++;		
		if(flkStr_can_or_not(tmpSJ_flankString_raw))
		{
			SJ_can ++;
			if(tmpSJ_supNum >= 3)
				SJ_can_min3 ++;
			if(tmpSJ_supNum >= 5)
				SJ_can_min5 ++;
		}
		else if(flkStr_sem_or_not(tmpSJ_flankString_raw))
		{
			SJ_sem ++;
			if(tmpSJ_supNum >= 3)
				SJ_sem_min3 ++;
			if(tmpSJ_supNum >= 5)
				SJ_sem_min5 ++;			
		}
		else
		{
			SJ_non ++;
			if(tmpSJ_supNum >= 3)
				SJ_non_min3 ++;
			if(tmpSJ_supNum >= 5)
				SJ_non_min5 ++;			
		}
		if((tmpSJ_flankString_raw == tmpSJ_flankString_hap1)&&(tmpSJ_flankString_hap1 == tmpSJ_flankString_hap2)) // no change
		{
			SJ_flkStrKept ++;
			if(tmpSJ_supNum >= 3)
				SJ_flkStrKept_min3 ++;
			if(tmpSJ_supNum >= 5)
				SJ_flkStrKept_min5 ++;
			kept_ofs << tmpSJinfoStr << endl;
		}			
		else // flank string changed
		{
			SJ_flkStrChanged ++;
			if(tmpSJ_supNum >= 3)
				SJ_flkStrChanged_min3 ++;
			if(tmpSJ_supNum >= 5)
				SJ_flkStrChanged_min5 ++;
			//changed_ofs << tmpSJinfoStr << endl;
			if(tmpSJ_flankString_hap1 == tmpSJ_flankString_hap2) // hom
			{
				SJ_flkStrChanged_hom ++;
				if(tmpSJ_supNum >= 3)
					SJ_flkStrChanged_hom_min3 ++;
				if(tmpSJ_supNum >= 5)
					SJ_flkStrChanged_hom_min5 ++;				
				//changed_hom_ofs << tmpSJinfoStr << endl;
				if(flkStr_can_or_not(tmpSJ_flankString_raw)) // from can
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)) // to can
					{
						SJ_flkStrChanged_can2can ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2can_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2can_min5 ++;
						//can_2_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Can2Can" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)) // to sem
					{
						SJ_flkStrChanged_can2sem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2sem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2sem_min5 ++;
						//can_2_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Can2Sem" << endl;
					}
					else // to non
					{
						SJ_flkStrChanged_can2non ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2non_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2non_min5 ++;
						//can_2_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Can2Non" << endl;
					}
				}
				else if(flkStr_sem_or_not(tmpSJ_flankString_raw)) // from sem
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)) // to can
					{
						SJ_flkStrChanged_sem2can ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2can_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2can_min5 ++;
						//sem_2_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Sem2Can" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)) // to sem
					{
						SJ_flkStrChanged_sem2sem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2sem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2sem_min5 ++;
						//sem_2_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Sem2Sem" << endl;
					}
					else // to non
					{
						SJ_flkStrChanged_sem2non ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2non_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2non_min5 ++;
						//sem_2_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Sem2Non" << endl;
					}
				}
				else // from non
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)) // to can
					{
						SJ_flkStrChanged_non2can ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2can_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2can_min5 ++;						
						//non_2_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Non2Can" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)) // to sem
					{
						SJ_flkStrChanged_non2sem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2sem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2sem_min5 ++;
						//non_2_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Non2Sem" << endl;
					}
					else // to non
					{
						SJ_flkStrChanged_non2non ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2non_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2non_min5 ++;						
						//non_2_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHom_Non2Non" << endl;
					}					
				}
			}
			else // het
			{
				SJ_flkStrChanged_het ++;
				if(tmpSJ_supNum >= 3)
					SJ_flkStrChanged_het_min3 ++;
				if(tmpSJ_supNum >= 5)
					SJ_flkStrChanged_het_min5 ++;
				//changed_het_ofs << tmpSJinfoStr << endl;
				if(flkStr_can_or_not(tmpSJ_flankString_raw)) // from can
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2)) // to can / can
					{
						SJ_flkStrChanged_can2canORcan ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2canORcan_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2canORcan_min5 ++;			
						//can_2_can_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Can/Can" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))
						||(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / sem
					{
						SJ_flkStrChanged_can2canORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2canORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2canORsem_min5 ++;						
						//can_2_can_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Can/Sem" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / non
					{
						SJ_flkStrChanged_can2canORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2canORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2canORnon_min5 ++;
						//can_2_can_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Can/Non" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2)) // to sem / sem
					{
						SJ_flkStrChanged_can2semORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2semORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2semORsem_min5 ++;					
						//can_2_sem_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Sem/Sem" << endl;
					}
					else if((flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))) // to sem / non
					{
						SJ_flkStrChanged_can2semORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2semORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2semORnon_min5 ++;
						//can_2_sem_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Sem/Non" << endl;
					}
					else // to non / non
					{
						SJ_flkStrChanged_can2nonORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_can2nonORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_can2nonORnon_min5 ++;
						//can_2_non_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Can2Non/Non" << endl;
					}
				}
				else if(flkStr_sem_or_not(tmpSJ_flankString_raw)) // from sem
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2)) // to can / can
					{
						SJ_flkStrChanged_sem2canORcan ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2canORcan_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2canORcan_min5 ++;					
						//sem_2_can_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Can/Can" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))
						||(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / sem
					{
						SJ_flkStrChanged_sem2canORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2canORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2canORsem_min5 ++;
						//sem_2_can_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Can/Sem" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / non
					{
						SJ_flkStrChanged_sem2canORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2canORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2canORnon_min5 ++;
						//sem_2_can_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Can/Non" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2)) // to sem / sem
					{
						SJ_flkStrChanged_sem2semORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2semORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2semORsem_min5 ++;
						//sem_2_sem_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Sem/Sem" << endl;
					}
					else if((flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))) // to sem / non
					{
						SJ_flkStrChanged_sem2semORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2semORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2semORnon_min5 ++;
						//sem_2_sem_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Sem/Non" << endl;
					}
					else // to non / non
					{
						SJ_flkStrChanged_sem2nonORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_sem2nonORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_sem2nonORnon_min5 ++;
						//sem_2_non_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Sem2Non/Non" << endl;
					}
				}
				else // from non
				{
					if(flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2)) // to can / can
					{
						SJ_flkStrChanged_non2canORcan ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2canORcan_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2canORcan_min5 ++;
						//non_2_can_can_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Can/Can" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))
						||(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / sem
					{
						SJ_flkStrChanged_non2canORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2canORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2canORsem_min5 ++;						
						//non_2_can_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Can/Sem" << endl;
					}
					else if((flkStr_can_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_can_or_not(tmpSJ_flankString_hap2))) // to can / non
					{
						SJ_flkStrChanged_non2canORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2canORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2canORnon_min5 ++;
						//non_2_can_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Can/Non" << endl;
					}
					else if(flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2)) // to sem / sem
					{
						SJ_flkStrChanged_non2semORsem ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2semORsem_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2semORsem_min5 ++;
						//non_2_sem_sem_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Sem/Sem" << endl;
					}
					else if((flkStr_sem_or_not(tmpSJ_flankString_hap1)&&flkStr_non_or_not(tmpSJ_flankString_hap2))
						||(flkStr_non_or_not(tmpSJ_flankString_hap1)&&flkStr_sem_or_not(tmpSJ_flankString_hap2))) // to sem / non
					{
						SJ_flkStrChanged_non2semORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2semORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2semORnon_min5 ++;
						//non_2_sem_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Sem/Non" << endl;
					}
					else // to non / non
					{
						SJ_flkStrChanged_non2nonORnon ++;
						if(tmpSJ_supNum >= 3)
							SJ_flkStrChanged_non2nonORnon_min3 ++;
						if(tmpSJ_supNum >= 5)
							SJ_flkStrChanged_non2nonORnon_min5 ++;
						//non_2_non_non_ofs << tmpSJinfoStr << endl;
						changed_ofs << tmpSJinfoStr << "\tHet_Non2Non/Non" << endl;
					}
				}
			}
		}
	}
	SJ_ifs.close();
	log_ofs << "All jobs done!" << endl;

	stats_ofs << "Junction Total #:\t" << SJ_all << "\t" << SJ_all_min3 << "\t" << SJ_all_min5 << endl << endl;

	stats_ofs << "Junction Can #:\t" << SJ_can << "\t" << SJ_can_min3 << "\t" << SJ_can_min5 << endl;
	stats_ofs << "Junction Sem #:\t" << SJ_sem << "\t" << SJ_sem_min3 << "\t" << SJ_sem_min5 << endl;
	stats_ofs << "Junction Non #:\t" << SJ_non << "\t" << SJ_non_min3 << "\t" << SJ_non_min5 << endl << endl;

	stats_ofs << "Junction FlkStr Kept #:\t" << SJ_flkStrKept << "\t" << SJ_flkStrKept_min3 << "\t" << SJ_flkStrKept_min5 << endl << endl;
	stats_ofs << "Junction FlkStr Changed #:\t" << SJ_flkStrChanged << "\t" << SJ_flkStrChanged_min3 << "\t" << SJ_flkStrChanged_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Hom #:\t" << SJ_flkStrChanged_hom << "\t" << SJ_flkStrChanged_hom_min3 << "\t" << SJ_flkStrChanged_hom_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Het #:\t" << SJ_flkStrChanged_het << "\t" << SJ_flkStrChanged_het_min3 << "\t" << SJ_flkStrChanged_het_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Can #:\t" << SJ_flkStrChanged_can2can << "\t" << SJ_flkStrChanged_can2can_min3 << "\t" << SJ_flkStrChanged_can2can_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Sem #:\t" << SJ_flkStrChanged_can2sem << "\t" << SJ_flkStrChanged_can2sem_min3 << "\t" << SJ_flkStrChanged_can2sem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Non #:\t" << SJ_flkStrChanged_can2non << "\t" << SJ_flkStrChanged_can2non_min3 << "\t" << SJ_flkStrChanged_can2non_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Can #:\t" << SJ_flkStrChanged_sem2can << "\t" << SJ_flkStrChanged_sem2can_min3 << "\t" << SJ_flkStrChanged_sem2can_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Sem #:\t" << SJ_flkStrChanged_sem2sem << "\t" << SJ_flkStrChanged_sem2sem_min3 << "\t" << SJ_flkStrChanged_sem2sem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Non #:\t" << SJ_flkStrChanged_sem2non << "\t" << SJ_flkStrChanged_sem2non_min3 << "\t" << SJ_flkStrChanged_sem2non_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Can #:\t" << SJ_flkStrChanged_non2can << "\t" << SJ_flkStrChanged_non2can_min3 << "\t" << SJ_flkStrChanged_non2can_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Sem #:\t" << SJ_flkStrChanged_non2sem << "\t" << SJ_flkStrChanged_non2sem_min3 << "\t" << SJ_flkStrChanged_non2sem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Non #:\t" << SJ_flkStrChanged_non2non << "\t" << SJ_flkStrChanged_non2non_min3 << "\t" << SJ_flkStrChanged_non2non_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Can/Can #:\t" << SJ_flkStrChanged_can2canORcan << "\t" << SJ_flkStrChanged_can2canORcan_min3 << "\t" << SJ_flkStrChanged_can2canORcan_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Can/Sem #:\t" << SJ_flkStrChanged_can2canORsem << "\t" << SJ_flkStrChanged_can2canORsem_min3 << "\t" << SJ_flkStrChanged_can2canORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Can/Non #:\t" << SJ_flkStrChanged_can2canORnon << "\t" << SJ_flkStrChanged_can2canORnon_min3 << "\t" << SJ_flkStrChanged_can2canORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Sem/Sem #:\t" << SJ_flkStrChanged_can2semORsem << "\t" << SJ_flkStrChanged_can2semORsem_min3 << "\t" << SJ_flkStrChanged_can2semORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Sem/Non #:\t" << SJ_flkStrChanged_can2semORnon << "\t" << SJ_flkStrChanged_can2semORnon_min3 << "\t" << SJ_flkStrChanged_can2semORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Can to Non/Non #:\t" << SJ_flkStrChanged_can2nonORnon << "\t" << SJ_flkStrChanged_can2nonORnon_min3 << "\t" << SJ_flkStrChanged_can2nonORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Can/Can #:\t" << SJ_flkStrChanged_sem2canORcan << "\t" << SJ_flkStrChanged_sem2canORcan_min3 << "\t" << SJ_flkStrChanged_sem2canORcan_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Can/Sem #:\t" << SJ_flkStrChanged_sem2canORsem << "\t" << SJ_flkStrChanged_sem2canORsem_min3 << "\t" << SJ_flkStrChanged_sem2canORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Can/Non #:\t" << SJ_flkStrChanged_sem2canORnon << "\t" << SJ_flkStrChanged_sem2canORnon_min3 << "\t" << SJ_flkStrChanged_sem2canORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Sem/Sem #:\t" << SJ_flkStrChanged_sem2semORsem << "\t" << SJ_flkStrChanged_sem2semORsem_min3 << "\t" << SJ_flkStrChanged_sem2semORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Sem/Non #:\t" << SJ_flkStrChanged_sem2semORnon << "\t" << SJ_flkStrChanged_sem2semORnon_min3 << "\t" << SJ_flkStrChanged_sem2semORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Sem to Non/Non #:\t" << SJ_flkStrChanged_sem2nonORnon << "\t" << SJ_flkStrChanged_sem2nonORnon_min3 << "\t" << SJ_flkStrChanged_sem2nonORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Can/Can #:\t" << SJ_flkStrChanged_non2canORcan << "\t" << SJ_flkStrChanged_non2canORcan_min3 << "\t" << SJ_flkStrChanged_non2canORcan_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Can/Sem #:\t" << SJ_flkStrChanged_non2canORsem << "\t" << SJ_flkStrChanged_non2canORsem_min3 << "\t" << SJ_flkStrChanged_non2canORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Can/Non #:\t" << SJ_flkStrChanged_non2canORnon << "\t" << SJ_flkStrChanged_non2canORnon_min3 << "\t" << SJ_flkStrChanged_non2canORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Sem/Sem #:\t" << SJ_flkStrChanged_non2semORsem << "\t" << SJ_flkStrChanged_non2semORsem_min3 << "\t" << SJ_flkStrChanged_non2semORsem_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Sem/Non #:\t" << SJ_flkStrChanged_non2semORnon << "\t" << SJ_flkStrChanged_non2semORnon_min3 << "\t" << SJ_flkStrChanged_non2semORnon_min5 << endl;
	stats_ofs << "Junction FlkStr Changed Non to Non/Non #:\t" << SJ_flkStrChanged_non2nonORnon << "\t" << SJ_flkStrChanged_non2nonORnon_min3 << "\t" << SJ_flkStrChanged_non2nonORnon << endl;	
	
	stats_ofs << endl;
	stats_ofs << "Can #:\t" << SJ_can << "\t" << SJ_can_min3 << "\t" << SJ_can_min5 << endl;
	stats_ofs << "Can->Non #:\t" << SJ_flkStrChanged_can2non + SJ_flkStrChanged_can2nonORnon 
		<< "\t" << SJ_flkStrChanged_can2non_min3 + SJ_flkStrChanged_can2nonORnon_min3
		<< "\t" << SJ_flkStrChanged_can2non_min5 + SJ_flkStrChanged_can2nonORnon_min5 << endl;
	stats_ofs << "Can->Can/Non #:\t" << SJ_flkStrChanged_can2canORnon << "\t" 
		<< SJ_flkStrChanged_can2canORnon_min3 << "\t" << SJ_flkStrChanged_can2canORnon_min5 << endl;
	stats_ofs << "Non #:\t" << SJ_non << "\t" << SJ_non_min3 << "\t" << SJ_non_min5 << endl; 
	stats_ofs << "Non->Can #:\t" << SJ_flkStrChanged_non2can + SJ_flkStrChanged_non2canORcan 
		<< "\t" << SJ_flkStrChanged_non2can_min3 + SJ_flkStrChanged_non2canORcan_min3
		<< "\t" << SJ_flkStrChanged_non2can_min5 + SJ_flkStrChanged_non2canORcan_min5 << endl;
	stats_ofs << "Non->Can/Non #:\t" << SJ_flkStrChanged_non2canORnon << "\t"
		<< SJ_flkStrChanged_non2canORnon_min3 << "\t" << SJ_flkStrChanged_non2canORnon_min5 << endl;
	all_ofs.close();
	kept_ofs.close();
	changed_ofs.close();
	// changed_hom_ofs.close();
	// changed_het_ofs.close();
	// can_2_can_ofs.close();
	// can_2_sem_ofs.close();
	// can_2_non_ofs.close();
	// sem_2_can_ofs.close();
	// sem_2_sem_ofs.close();
	// sem_2_non_ofs.close();
	// non_2_can_ofs.close();
	// non_2_sem_ofs.close();
	// non_2_non_ofs.close();	
	// can_2_can_can_ofs.close();
	// can_2_can_sem_ofs.close();
	// can_2_can_non_ofs.close();
	// can_2_sem_sem_ofs.close();
	// can_2_sem_non_ofs.close();
	// can_2_non_non_ofs.close();
	// sem_2_can_can_ofs.close();
	// sem_2_can_sem_ofs.close();
	// sem_2_can_non_ofs.close();
	// sem_2_sem_sem_ofs.close();
	// sem_2_sem_non_ofs.close();
	// sem_2_non_non_ofs.close();
	// non_2_can_can_ofs.close();
	// non_2_can_sem_ofs.close();
	// non_2_can_non_ofs.close();
	// non_2_sem_sem_ofs.close();
	// non_2_sem_non_ofs.close();
	// non_2_non_non_ofs.close();

	delete indexInfo;
	delete indexInfo_hap1;
	delete indexInfo_hap2;
	log_ofs.close();
	return 0;
}