/*    
 *    alignAll.cpp
 *	  MapSplice3
 *
 *    Copyright (C) 2016 University of Kentucky and
 *                       Xinan Liu
 *
 *    Authors: Xinan Liu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"

using namespace std;

int getMaxSegLength(bool getSegInfo_bool_Nor1, Seg_Info* segInfo_Nor1, 
	bool getSegInfo_bool_Rcm1, Seg_Info* segInfo_Rcm1,
	bool getSegInfo_bool_Nor2, Seg_Info* segInfo_Nor2, 
	bool getSegInfo_bool_Rcm2, Seg_Info* segInfo_Rcm2)//, int& maxSegHitNum)
{
	int tmpMaxLength = 0;	
	if(getSegInfo_bool_Nor1)
	{	
		int tmpMaxLength_Nor1 = segInfo_Nor1->getMaxSegLength();
		if(tmpMaxLength_Nor1 > tmpMaxLength)
			tmpMaxLength = tmpMaxLength_Nor1;
	}
	if(getSegInfo_bool_Rcm1)
	{	
		int tmpMaxLength_Rcm1 = segInfo_Rcm1->getMaxSegLength();
		if(tmpMaxLength_Rcm1 > tmpMaxLength)
			tmpMaxLength = tmpMaxLength_Rcm1;
	}
	if(getSegInfo_bool_Nor2)
	{	
		int tmpMaxLength_Nor2 = segInfo_Nor2->getMaxSegLength();
		if(tmpMaxLength_Nor2 > tmpMaxLength)
			tmpMaxLength = tmpMaxLength_Nor2;
	}
	if(getSegInfo_bool_Rcm2)
	{	
		int tmpMaxLength_Rcm2 = segInfo_Rcm2->getMaxSegLength();
		if(tmpMaxLength_Rcm2 > tmpMaxLength)	
			tmpMaxLength = tmpMaxLength_Rcm2;
	}
	return tmpMaxLength;
}


int main(int argc, char**argv)
{
	if((argc != 7)&&(argc != 8))
	{
		cout << "#0 Executable" << endl;
		cout << "#1 indexFolder" << endl;
		cout << "#2 outputDir" << endl;
		cout << "#3 threads_num" << endl;
		cout << "#4 fa_or_fq" << endl;
		cout << "#5 SE_or_PE" << endl;
		cout << "#6 inputRead_1" << endl;
		cout << "(#7 inputRead_2)" << endl;
		exit(1);
	}
	int readLength_max = 100;
	//int hit_log2_max = 10;
	vector<int> maxSegLengthReadNumVec;
	for(int tmp = 0; tmp <= readLength_max; tmp++)
		maxSegLengthReadNumVec.push_back(0);
	// vector< vector<int> > maxSegLengthHitNumVecVec;
	// for(int tmp = 0; tmp <= readLength_max; tmp++)
	// {
	// 	vector<int> tmpHitNumVec;
	// 	for(int tmp2 = 0; tmp2 <= hit_log2_max; tmp2 ++)
	// 		tmpHitNumVec.push_back(0);
	// 	maxSegLengthHitNumVecVec.push_back(tmpHitNumVec);	
	// }
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	string indexStr = argv[1];
	string outputDirStr = argv[2];
	string threadsNumStr = argv[3];
	string fa_or_fq_str = argv[4];
	string SE_or_PE_str = argv[5];
	string InputReadFile = argv[6];
	string InputReadFile_PE;
	bool SE_or_PE_bool;
	if(SE_or_PE_str == "SE")
		SE_or_PE_bool = true;
	else
	{
		SE_or_PE_bool = false;	
		InputReadFile_PE = argv[7];
	}
	bool InputAsFastq;
	if((fa_or_fq_str == "fa")||(fa_or_fq_str == "Fa")||(fa_or_fq_str == "FA")
		||(fa_or_fq_str == "fasta")||(fa_or_fq_str == "Fasta")||(fa_or_fq_str == "FASTA"))
		InputAsFastq = false;
	else if((fa_or_fq_str == "fq")||(fa_or_fq_str == "Fq")||(fa_or_fq_str == "FQ")
		||(fa_or_fq_str == "fastq")||(fa_or_fq_str == "Fastq")||(fa_or_fq_str == "FASTQ"))
		InputAsFastq = true;
	bool fasta_or_fastq_bool = (!InputAsFastq);
	/////////////////////////////////////////////////
	int normalRecordNum = 500000; //1000000;//1500000;
	int readTotalNum = 0;
	int threads_num = atoi(threadsNumStr.c_str());
	omp_set_num_threads(threads_num);

	outputDirStr += "/";
	string cmd_mkdir = "mkdir -p " + outputDirStr;
	system(cmd_mkdir.c_str());
	string log_file = outputDirStr + "log.txt";
	string settings_file = outputDirStr + "settings.txt";
	ofstream log_ofs(log_file.c_str());
	ofstream settings_ofs(settings_file.c_str());
	settings_ofs << "#0 Executable" << endl << "#1 indexFolder: " << indexStr << endl 
		<< "#2 outputDir: " << outputDirStr << endl << "#3 threads_num: " << threadsNumStr << endl
		<< "#4 fa_or_fq: " << fa_or_fq_str << endl << "#5 SE_or_PE: " << SE_or_PE_str << endl
		<< "#6 InputReadFile: " << InputReadFile << endl;  
	if(!SE_or_PE_bool)
		settings_ofs << "#7 InputReadFile_PE: " << InputReadFile_PE << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////          LOAD INDEX         /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	log_ofs << "start to load index" << endl;
    string preIndexArrayPreStr = indexStr;
    preIndexArrayPreStr.append("/");
    indexStr.append("/");

	string preIndexMapLengthArrayStr = preIndexArrayPreStr; preIndexMapLengthArrayStr.append("_MapLength"); ifstream preIndexMapLengthArray_ifs(preIndexMapLengthArrayStr.c_str(), ios::binary);
	string preIndexIntervalStartArrayStr = preIndexArrayPreStr; preIndexIntervalStartArrayStr.append("_IntervalStart"); ifstream preIndexIntervalStartArray_ifs(preIndexIntervalStartArrayStr.c_str(), ios::binary);
	string preIndexIntervalEndArrayStr = preIndexArrayPreStr; preIndexIntervalEndArrayStr.append("_IntervalEnd"); ifstream preIndexIntervalEndArray_ifs(preIndexIntervalEndArrayStr.c_str(), ios::binary);
	int* preIndexMapLengthArray; preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int)); preIndexMapLengthArray_ifs.read((char*)preIndexMapLengthArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalStartArray; preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalStartArray_ifs.read((char*)preIndexIntervalStartArray, PreIndexSize * sizeof(int));
	unsigned int *preIndexIntervalEndArray; preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int)); preIndexIntervalEndArray_ifs.read((char*)preIndexIntervalEndArray, PreIndexSize * sizeof(int));

	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);
	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.close();
	parameter_file_ifs.close();

	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
 	
	string SA_file = indexStr; SA_file.append("_SA"); 
	string lcpCompress_file = indexStr; lcpCompress_file.append("_lcpCompress"); 
	string childTab_file = indexStr; childTab_file.append("_childTab"); 
	string verifyChild_file = indexStr; verifyChild_file.append("_detChild"); 	
    unsigned int *sa; sa = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexInfo->returnIndexSize()) * sizeof(unsigned int)); 
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexInfo->returnIndexSize()) * sizeof(BYTE)); 	
	ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
	SA_file_ifs.read((char*)sa, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	lcpCompress_file_ifs.read((char*)lcpCompress, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	childTab_file_ifs.read((char*)childTab, (indexInfo->returnIndexSize()) * sizeof(unsigned int));
	verifyChild_file_ifs.read((char*)verifyChild, (indexInfo->returnIndexSize()) * sizeof(BYTE));
	SA_file_ifs.close();
	lcpCompress_file_ifs.close();
	childTab_file_ifs.close();
	verifyChild_file_ifs.close();
	log_ofs << "end of loading index" << endl;
	//////////////////////////////////////////////////       finish LOADing INDEX           ///////////////////////
	//////////////////////////////////////////////////       finish LOADing INDEX           ///////////////////////
	settings_ofs << "start to get maxSegLength for each read" << endl;
	string maxSegLength_file = outputDirStr + "/maxSegLength.txt";
	string stats_file = outputDirStr + "/stats.txt";
	ifstream inputRead_ifs(InputReadFile.c_str());
	ifstream inputRead_PE_ifs(InputReadFile_PE.c_str());
	ofstream maxSegLength_ofs(maxSegLength_file.c_str());
	ofstream stats_ofs(stats_file.c_str());
	if(SE_or_PE_bool)
		inputRead_PE_ifs.close();
	else
	{}

    string line1, line2, line3, line4, line1_PE, line2_PE, line3_PE, line4_PE;
    string line2_afterProcess, line2_PE_afterProcess;
	bool EndOfRecord = false;
	int tmpTurn = 0;
	int realRecordNum;
	int readPairNum = 0;

	// *************** needed for both SE and PE *************//
	vector<string> readName1Vec(normalRecordNum);
	vector<string> readSeq1Vec(normalRecordNum);
	vector<string> readQualSeq1Vec(normalRecordNum);
	// *******************************************************//
	// ***************  only needed for PE  ******************//
	vector<string> readName2Vec(normalRecordNum);
	vector<string> readSeq2Vec(normalRecordNum);
	vector<string> readQualSeq2Vec(normalRecordNum);
	//vector<int> segNumVec(normalRecordNum);
	vector<int> maxSegLengthVec(normalRecordNum);
	//vector<int> maxSegHitNumVec(normalRecordNum);

	InputReadPreProcess* readPreProcessInfo = new InputReadPreProcess();
	for(tmpTurn = 0;
		//tmpTurn <= 300     //used to control # of rounds to process
		; tmpTurn++)
	{
		if(EndOfRecord)
			break;
		int recordNum = normalRecordNum;
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << "[" << asctime(local) << "start to read Fasta file, turn: " << tmpTurn + 1 << endl;
		realRecordNum = normalRecordNum;

		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		getline(inputRead_ifs, line1); // readName_1
    		if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
    		{
				realRecordNum = recordNumTmp;
				EndOfRecord = true;
				break;    			
    		}
    		readName1Vec[recordNumTmp] = line1.substr(1);
    		getline(inputRead_ifs, line2); // readSeq_1
    		line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
    		readSeq1Vec[recordNumTmp] = line2_afterProcess;
    		if(InputAsFastq)
    		{
    			getline(inputRead_ifs, line3);
    			getline(inputRead_ifs, line4);
    			readQualSeq1Vec[recordNumTmp] = line4;   
    		}
    			
    		if(!SE_or_PE_bool)
    		{	
	    		getline(inputRead_PE_ifs, line1_PE); // readName_2
	    		readName2Vec[recordNumTmp] = line1_PE.substr(1);
	    		getline(inputRead_PE_ifs, line2_PE); // readSeq_2
	    		line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
	    		readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;
	    		if(InputAsFastq)
	    		{
	    			getline(inputRead_PE_ifs, line3_PE);
	    			getline(inputRead_PE_ifs, line4_PE);
	    			readQualSeq2Vec[recordNumTmp] = line4_PE;
	    		}
    		}
		}
		readTotalNum += realRecordNum;

		log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
		nowtime = time(NULL);
		local = localtime(&nowtime);		
		log_ofs << endl << "[" << asctime(local) << "finish reading Fasta/Fastq file, turn: " << tmpTurn + 1 << endl;
		log_ofs << "start to fix reads in phase1, turn: " << tmpTurn + 1 << endl;

		omp_set_num_threads(threads_num);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
		{
			int threadNO = omp_get_thread_num();

			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(readName1Vec[tmpOpenMP], readName2Vec[tmpOpenMP],
				readSeq1Vec[tmpOpenMP], readSeq2Vec[tmpOpenMP],
				readQualSeq1Vec[tmpOpenMP], readQualSeq2Vec[tmpOpenMP], fasta_or_fastq_bool, SE_or_PE_bool);
	
			//cout << endl << "read_name_1: " << readName1Vec[tmpOpenMP] << endl;
			//cout << "read_name_2: " << readName2Vec[tmpOpenMP] << endl;
			//cout << "readSeq_1: " << readSeq1Vec[tmpOpenMP] << endl;
			//cout << "readSeq_2: " << readSeq2Vec[tmpOpenMP] << endl;
			Seg_Info* segInfo_Nor1 = new Seg_Info(); 
			Seg_Info* segInfo_Rcm1 = new Seg_Info(); 
			Seg_Info* segInfo_Nor2 = new Seg_Info(); 
			Seg_Info* segInfo_Rcm2 = new Seg_Info();
			char* read_Nor1 = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());
			bool getSegInfo_bool_Nor1 = segInfo_Nor1->mapMain_preIndex_getSegInfo(read_Nor1, sa, lcpCompress, childTab, chrom,
				verifyChild, readInfo.returnReadSeqLength_1(), indexInfo, preIndexMapLengthArray, 
				preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo.returnReadSeq_1());
			char* read_Rcm1 = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
			bool getSegInfo_bool_Rcm1 = segInfo_Rcm1->mapMain_preIndex_getSegInfo(read_Rcm1, sa, lcpCompress, childTab, chrom,
				verifyChild, readInfo.returnReadSeqLength_1(), indexInfo, preIndexMapLengthArray, 
				preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo.returnRcmReadSeq_1());
			char* read_Nor2 = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
			bool getSegInfo_bool_Nor2 = segInfo_Nor2->mapMain_preIndex_getSegInfo(read_Nor2, sa, lcpCompress, childTab, chrom,
				verifyChild, readInfo.returnReadSeqLength_2(), indexInfo, preIndexMapLengthArray, 
				preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo.returnReadSeq_2());
			char* read_Rcm2 = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
			bool getSegInfo_bool_Rcm2 = segInfo_Rcm2->mapMain_preIndex_getSegInfo(read_Rcm2, sa, lcpCompress, childTab, chrom,
				verifyChild, readInfo.returnReadSeqLength_2(), indexInfo, preIndexMapLengthArray, 
				preIndexIntervalStartArray, preIndexIntervalEndArray, readInfo.returnRcmReadSeq_2());
			//int tmpMaxSegHitNum = 0;
			int tmpMaxSegLength = getMaxSegLength(getSegInfo_bool_Nor1, segInfo_Nor1, getSegInfo_bool_Rcm1,
				 segInfo_Rcm1, getSegInfo_bool_Nor2, segInfo_Nor2, getSegInfo_bool_Rcm2, segInfo_Rcm2);//, tmpMaxSegHitNum);
			maxSegLengthVec[tmpOpenMP] = tmpMaxSegLength;
			//maxSegHitNumVec[tmpOpenMP] = tmpMaxSegHitNum;
			delete segInfo_Nor1;
			delete segInfo_Rcm1;
			delete segInfo_Nor2;
			delete segInfo_Rcm2;			
		} // read file end
		
		nowtime = time(NULL);
		local = localtime(&nowtime);
		log_ofs << endl << "[" << asctime(local) << "finish fixing reads in phase1, turn: " << tmpTurn + 1 << endl;
		log_ofs << "start to output ... turn: " << tmpTurn+1 << endl;
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			maxSegLength_ofs << readName1Vec[tmpOpenMP] 
				<< "\t" << maxSegLengthVec[tmpOpenMP] << endl;// << "\t" << maxSegHitNumVec[tmpOpenMP] << endl;
		for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
		{
			int tmpMaxSegLength = maxSegLengthVec[tmpOpenMP];
			// int tmpMaxSegHitNum = maxSegHitNumVec[tmpOpenMP];
			// int tmpMaxSegHitNum_log2 = log2(tmpMaxSegHitNum);
			maxSegLengthReadNumVec[tmpMaxSegLength] ++;
			// if(tmpMaxSegHitNum_log2 >= hit_log2_max)
			// 	tmpMaxSegHitNum_log2 = hit_log2_max;
			// (maxSegLengthHitNumVecVec[tmpMaxSegLength])[tmpMaxSegHitNum_log2] ++;
		}
	}
	//stats_ofs << "maxSegLength\t(total)%\t(1)%\t(2~3)%\t(4~7)%\t(8~15)%\t(16~31)%\t(32~63)%\t(64~127)%\t(128~255)%\t(256~511)%\t(512~1023)%\t(>=1024)%" << endl;
	for(int tmp = 0; tmp <= readLength_max; tmp++)
	{
		int tmpReadNum = maxSegLengthReadNumVec[tmp];
		double tmpReadPerc = ((double)tmpReadNum/(double)readTotalNum) * 100;
		stats_ofs << tmp << "\t" << tmpReadPerc << endl;
		// if(tmpReadNum == 0)
		// {
		// 	for(int tmp2 = 0; tmp2 <= hit_log2_max; tmp2++)
		// 		stats_ofs << "\tNULL";
		// }
		// else
		// {
		// 	for(int tmp2 = 0; tmp2 <= hit_log2_max; tmp2++)
		// 	{	
		// 		double tmpHitNum = (maxSegLengthHitNumVecVec[tmp])[tmp2];
		// 		double tmpHitPerc = ((double)tmpHitNum/(double)tmpReadNum) * 100;
		// 		stats_ofs << "\t" << tmpHitPerc;
		// 	}
		// }
		// stats_ofs << endl;
	}
	settings_ofs << "end of getting maxSegLength for all reads." << endl;
	delete readPreProcessInfo;
	settings_ofs << endl << "readTotalNum: " << readTotalNum << endl;
	inputRead_ifs.close();
	inputRead_PE_ifs.close();
	maxSegLength_ofs.close();
	stats_ofs.close();
	free(preIndexMapLengthArray); 
	free(preIndexIntervalStartArray); 
	free(preIndexIntervalEndArray);
	free(sa);free(lcpCompress);
	free(childTab);
	free(verifyChild);
	free(chrom);
	delete indexInfo;
	nowtime = time(NULL);
	local = localtime(&nowtime);
	//cout << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... all jobs done ......" << endl << endl;
	log_ofs.close();
	settings_ofs.close();
    return 0;
} //end main