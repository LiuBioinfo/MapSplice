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
#include "../../../../general/read_block_test.h"
#include "../../../../general/otherFunc.h"
#include "../../../../general/index_info.h"
time_t nowtime;
struct tm *local;

using namespace std;

void get_chrNameIndexAndPos_from_wholeGenomePos(unsigned int tmpWholeGenomePos, int& tmpChrNameIndex,
	unsigned int& tmpChrPos, vector<unsigned int>& wholeGenomePosVec)
{
	if(tmpWholeGenomePos <= wholeGenomePosVec[0])
	{
		tmpChrNameIndex = 0;
		tmpChrPos = tmpWholeGenomePos;
		return;
	}

	for(int tmp = 1; tmp < wholeGenomePosVec.size(); tmp++)
	{
		unsigned int tmpWholeGenomePos_start = wholeGenomePosVec[tmp-1] + 1;
		unsigned int tmpWholeGenomePos_end = wholeGenomePosVec[tmp];
		if((tmpWholeGenomePos >= tmpWholeGenomePos_start)&&(tmpWholeGenomePos <= tmpWholeGenomePos_end))
		{
			tmpChrNameIndex = tmp;
			tmpChrPos = tmpWholeGenomePos - tmpWholeGenomePos_start + 1;
			return;
		}
	}

	cout << "error! invalid wholeGenomePos: " << tmpWholeGenomePos << endl;
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputChrSeqFaFileListFile inputClassificationGroupNum outputFolder" << endl;
		exit(1);
	}
	string inputChrSeqFaFileListFile = argv[1];
	vector<string> chrSeqFaPathVec;
	ifstream chrSeqFaFileList_ifs(inputChrSeqFaFileListFile.c_str());
	while(!chrSeqFaFileList_ifs.eof())
	{
		string tmpStr;
		getline(chrSeqFaFileList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		chrSeqFaPathVec.push_back(tmpStr);
	}
	chrSeqFaFileList_ifs.close();

	string inputClassificationGroupNumStr = argv[2];
	int inputClassificationGroupNum = atoi(inputClassificationGroupNumStr.c_str());
	int genomeSubRegionNum = inputClassificationGroupNum - 2; // shared kmer group & alien kmer group
	int overlapAreaSize = 500000;
	string outputDirStr = argv[3];
	outputDirStr += "/";
  	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());
   	string log_file = outputDirStr + "log.txt";
   	ofstream log_ofs(log_file.c_str());
   	string stats_file = outputDirStr + "stats.txt";
   	ofstream stats_ofs(stats_file.c_str());

   	cout << "start to read genome fa" << endl;
   	vector<string> chrNameVec;
   	vector<string> chrSeqVec;
   	vector<unsigned int> chrSeqLengthVec;
   	vector<unsigned int> wholeGenomePosVec;
   	unsigned int wholeGenomeLength;
   	for(int tmp = 0; tmp < chrSeqFaPathVec.size(); tmp++)
   	{
   		string tmpChrSeqFaPath = chrSeqFaPathVec[tmp];
   		ifstream tmpChrSeqFa_ifs(tmpChrSeqFaPath.c_str());
   		string tmpChrName;
   		getline(tmpChrSeqFa_ifs, tmpChrName);
   		int tmpDotLoc = tmpChrName.find(".");
   		cout << "tmpChrName: " << tmpChrName.substr(1, tmpDotLoc - 1) << endl;
   		chrNameVec.push_back(tmpChrName.substr(1, tmpDotLoc - 1));
   		string tmpChrSeq = "";
   		while(!tmpChrSeqFa_ifs.eof())
   		{
   			string tmpStr;
   			getline(tmpChrSeqFa_ifs, tmpStr);
   			if(tmpStr == "")
   				break;
   			tmpChrSeq += tmpStr;
   		}
   		chrSeqVec.push_back(tmpChrSeq);
   		unsigned int tmpChrSeqLength = tmpChrSeq.length();
   		chrSeqLengthVec.push_back(tmpChrSeqLength);
   		wholeGenomeLength += tmpChrSeqLength;
   		wholeGenomePosVec.push_back(wholeGenomeLength);
   		tmpChrSeqFa_ifs.close();
   		log_ofs << "chrNameIndex:\t" << tmp << endl;
   		log_ofs << "chrName:\t" << chrNameVec[tmp] << endl;
   		log_ofs << "chrSeqLength:\t" << chrSeqLengthVec[tmp] << endl;
   		log_ofs << "endPosInWholeGenome:\t" << wholeGenomePosVec[tmp] << endl; 
   	}
   	unsigned int genomeSubRegionSize = wholeGenomeLength/genomeSubRegionNum;
   	//vector< pair<int, int> > chrNameIndexPairVec_start_end;
   	//vector< pair<unsigned int, unsigned int> > lengthPairVec_start_end;
	cout << "start to generate genomeSubRegion fa" << endl;
	log_ofs << "start to generate genomeSubRegion fa" << endl;
	for(int tmpRegionIndex = 0; tmpRegionIndex < genomeSubRegionNum; tmpRegionIndex++)
	{
		stats_ofs << "Region:\t" << tmpRegionIndex + 1 << endl;
		unsigned int tmpWholeGenomePos_start, tmpWholeGenomePos_end;
		if(tmpRegionIndex == 0)
		{
			tmpWholeGenomePos_start = 1;
			tmpWholeGenomePos_end = genomeSubRegionSize + overlapAreaSize;
		}
		else if(tmpRegionIndex == genomeSubRegionNum - 1)
		{
			tmpWholeGenomePos_start = (genomeSubRegionNum - 1) * genomeSubRegionSize + 1 - overlapAreaSize;
			tmpWholeGenomePos_end = wholeGenomeLength;
		}
		else
		{
			tmpWholeGenomePos_start = tmpRegionIndex * genomeSubRegionSize + 1 - overlapAreaSize; 
			tmpWholeGenomePos_end = (tmpRegionIndex + 1) * genomeSubRegionSize + overlapAreaSize;
		}
		int tmpChrNameIndex_start, tmpChrNameIndex_end; 
		unsigned int tmpChrPos_start, tmpChrPos_end;
		get_chrNameIndexAndPos_from_wholeGenomePos(tmpWholeGenomePos_start, tmpChrNameIndex_start, tmpChrPos_start, wholeGenomePosVec);
		get_chrNameIndexAndPos_from_wholeGenomePos(tmpWholeGenomePos_end, tmpChrNameIndex_end, tmpChrPos_end, wholeGenomePosVec);
		stats_ofs << "chrNameIndex_start:\t" << tmpChrNameIndex_start << endl;
		stats_ofs << "chrName_start:\t" << chrNameVec[tmpChrNameIndex_start] << endl;
		stats_ofs << "chrPos_start:\t" << tmpChrPos_start << endl;
		stats_ofs << "chrNameIndex_end:\t" << tmpChrNameIndex_end << endl;
		stats_ofs << "chrName_end:\t" << chrNameVec[tmpChrNameIndex_end] << endl;
		stats_ofs << "chrPos_end:\t" << tmpChrPos_end << endl;	
		string outputDirStr_subRegion = outputDirStr + "Region." + int_to_str(tmpRegionIndex + 1) + "/";
		string cmd_mkdir_subRegion = "mkdir " + outputDirStr_subRegion;
		system(cmd_mkdir_subRegion.c_str());
		if(tmpChrNameIndex_start == tmpChrNameIndex_end)
		{
			string tmpChrSeqFaFile = outputDirStr_subRegion + chrNameVec[tmpChrNameIndex_start] + ".fa";
			ofstream tmpChrSeqFa_ofs(tmpChrSeqFaFile.c_str());
			tmpChrSeqFa_ofs << ">" << chrNameVec[tmpChrNameIndex_start] << endl;
			string tmpRegionSeqFa = chrSeqVec[tmpChrNameIndex_start].substr(tmpChrPos_start - 1, genomeSubRegionSize);
			vector<string> tmpSeqVec;
			int fiftyBaseLineNum = genomeSubRegionSize/50;
			for(int tmpLindeIndex = 0; tmpLindeIndex < fiftyBaseLineNum; tmpLindeIndex++)
				tmpSeqVec.push_back(tmpRegionSeqFa.substr(tmpLindeIndex * 50, 50));
			if(fiftyBaseLineNum * 50 < genomeSubRegionSize)
			{
				int lastLineLength = genomeSubRegionSize - fiftyBaseLineNum * 50;
				tmpSeqVec.push_back(tmpRegionSeqFa.substr(fiftyBaseLineNum * 50, lastLineLength));
			}
			for(int tmpSeqIndex = 0; tmpSeqIndex < tmpSeqVec.size(); tmpSeqIndex++)
				tmpChrSeqFa_ofs << tmpSeqVec[tmpSeqIndex] << endl;
			tmpChrSeqFa_ofs.close();
		}
		else
		{
			// output startPos chrSeqFa
			string tmpChrSeqFaFile_start = outputDirStr_subRegion + chrNameVec[tmpChrNameIndex_start] + ".fa";
			ofstream tmpChrSeqFa_ofs_start(tmpChrSeqFaFile_start.c_str());
			tmpChrSeqFa_ofs_start << ">" << chrNameVec[tmpChrNameIndex_start] << endl;
			string tmpRegionSeqFa_start = chrSeqVec[tmpChrNameIndex_start].substr(tmpChrPos_start - 1);
			int tmpRegionSeqFaLength_start = tmpRegionSeqFa_start.length();
			vector<string> tmpSeqVec_start;
			int fiftyBaseLineNum_start = tmpRegionSeqFaLength_start/50;
			for(int tmpLindeIndex = 0; tmpLindeIndex < fiftyBaseLineNum_start; tmpLindeIndex++)
				tmpSeqVec_start.push_back(tmpRegionSeqFa_start.substr(tmpLindeIndex * 50, 50));
			if(fiftyBaseLineNum_start * 50 < tmpRegionSeqFaLength_start)
			{
				int lastLineLength = tmpRegionSeqFaLength_start - fiftyBaseLineNum_start * 50;
				tmpSeqVec_start.push_back(tmpRegionSeqFa_start.substr(fiftyBaseLineNum_start * 50, lastLineLength));
			}
			for(int tmpSeqIndex = 0; tmpSeqIndex < tmpSeqVec_start.size(); tmpSeqIndex++)
				tmpChrSeqFa_ofs_start << tmpSeqVec_start[tmpSeqIndex] << endl;
			tmpChrSeqFa_ofs_start.close();
			// output inter chrSeqFa
			if((tmpChrNameIndex_start + 1) <= (tmpChrNameIndex_end - 1))
			{
				for(int tmpChrIndex = tmpChrNameIndex_start + 1; tmpChrIndex <= tmpChrNameIndex_end - 1; tmpChrIndex++)
				{
					string tmpChrSeqFaFile = outputDirStr_subRegion + chrNameVec[tmpChrIndex] + ".fa";
					ofstream tmpChrSeqFa_ofs(tmpChrSeqFaFile.c_str());
					tmpChrSeqFa_ofs << ">" << chrNameVec[tmpChrIndex] << endl;
					string tmpRegionSeqFa = chrSeqVec[tmpChrIndex];
					int tmpChrSeqLength = chrSeqLengthVec[tmpChrIndex];
					vector<string> tmpSeqVec;
					int fiftyBaseLineNum  = tmpChrSeqLength/50;
					for(int tmpLindeIndex = 0; tmpLindeIndex < fiftyBaseLineNum; tmpLindeIndex++)
						tmpSeqVec.push_back(tmpRegionSeqFa.substr(tmpLindeIndex * 50, 50));
					if(fiftyBaseLineNum * 50 < tmpChrSeqLength)
					{
						int lastLineLength = tmpChrSeqLength - fiftyBaseLineNum * 50;
						tmpSeqVec.push_back(tmpRegionSeqFa.substr(fiftyBaseLineNum * 50, lastLineLength));
					}
					for(int tmpSeqIndex = 0; tmpSeqIndex < tmpSeqVec.size(); tmpSeqIndex++)
						tmpChrSeqFa_ofs << tmpSeqVec[tmpSeqIndex] << endl;
					tmpChrSeqFa_ofs.close();
				}
			}
			// output endPos chrSeqFa
			string tmpChrSeqFaFile_end = outputDirStr_subRegion + chrNameVec[tmpChrNameIndex_end] + ".fa";
			ofstream tmpChrSeqFa_ofs_end(tmpChrSeqFaFile_end.c_str());
			tmpChrSeqFa_ofs_end << ">" << chrNameVec[tmpChrNameIndex_end] << endl;
			string tmpRegionSeqFa_end = chrSeqVec[tmpChrNameIndex_end].substr(0, tmpChrPos_end);
			int tmpRegionSeqFaLength_end = tmpRegionSeqFa_end.length();
			vector<string> tmpSeqVec_end;
			int fiftyBaseLineNum_end = tmpRegionSeqFaLength_end/50;
			for(int tmpLindeIndex = 0; tmpLindeIndex < fiftyBaseLineNum_end; tmpLindeIndex++)
				tmpSeqVec_end.push_back(tmpRegionSeqFa_end.substr(tmpLindeIndex * 50, 50));
			if(fiftyBaseLineNum_end * 50 < tmpRegionSeqFaLength_end)
			{
				int lastLineLength = tmpRegionSeqFaLength_end - fiftyBaseLineNum_end * 50;
				tmpSeqVec_end.push_back(tmpRegionSeqFa_end.substr(fiftyBaseLineNum_end * 50, lastLineLength));
			}
			for(int tmpSeqIndex = 0; tmpSeqIndex < tmpSeqVec_end.size(); tmpSeqIndex++)
				tmpChrSeqFa_ofs_end << tmpSeqVec_end[tmpSeqIndex] << endl;
			tmpChrSeqFa_ofs_end.close();
		}
	}
	stats_ofs.close();
	//parameter_ifs.close();
	//chrom_bit_file_ifs.close();
	log_ofs.close();
	//delete indexInfo;
	//free(chrom);
	return 0;
}