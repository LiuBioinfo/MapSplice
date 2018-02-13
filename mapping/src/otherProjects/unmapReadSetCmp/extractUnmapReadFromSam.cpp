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

#include "../../general/otherFunc.h"

using namespace std;

void unmappedAlignemnt2ReadFile(const string& SAMfile,
	string& outputUnmappedReadFq_1, string& outputUnmappedReadFq_2, ofstream& output_log_ofs)
{
	string outputFileStr = SAMfile + ".unmapped.read";
	string outputFileStr_read_end1 = outputUnmappedReadFq_1;
	string outputFileStr_read_end2 = outputUnmappedReadFq_2; 	

	FILE *fp;
	fp = fopen(SAMfile.c_str(), "r");
	ofstream output_read_end1_ofs(outputFileStr_read_end1.c_str());
	ofstream output_read_end2_ofs(outputFileStr_read_end2.c_str());

	int oriAlignmentNum = 0;
	int unmappedAlignmentNum = 0;

	string tmpOriReadNameStr;
	char tmpSamStrChar[2000];
	char read_name[100], flag[3], rname[50],pos[20],mapq[3],cigar[100],rnext[50];
	char pnext[20], tlen[10], seq[700], qual[700], others[200];

	bool end1_or_end2_bool = false;
	while(!feof(fp))
	{
		fgets(tmpSamStrChar, sizeof(tmpSamStrChar),fp);
		if(feof(fp))
			break;	
		if(tmpSamStrChar[0] == '@')
			continue;
		oriAlignmentNum ++;

		string tmpOriSamLineStr = tmpSamStrChar;
		int firstTabPos = tmpOriSamLineStr.find("\t");// << endl;
		string tmpClippedSamLine = tmpOriSamLineStr.substr(firstTabPos);
		tmpClippedSamLine = "seq" + tmpClippedSamLine;
		tmpOriReadNameStr = tmpOriSamLineStr.substr(0, firstTabPos);

		sscanf(tmpClippedSamLine.c_str(),"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,others);  
		end1_or_end2_bool = (!end1_or_end2_bool);
		int tmpFlag = atoi(flag);
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(mappedOrNot_bool)
		{}
		else
		{
			unmappedAlignmentNum ++;
			string tmpReadName = tmpOriReadNameStr;
			string tmpReadSeq = seq;
			string tmpQualSeq = qual;
			if(tmpQualSeq.length() != tmpReadSeq.length())
			{
				tmpQualSeq = "I";
				for(int tmp = 0;  tmp < tmpReadSeq.length() - 1; tmp++)
				{
					tmpQualSeq += "I";
				}
			}
			if(end1_or_end2_bool)
			{
				output_read_end1_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
			else
			{
				output_read_end2_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
		}
	}

	output_log_ofs << "oriAlignmentNum: " << oriAlignmentNum << endl;
	output_log_ofs << "unmappedAlignmentNum: " << unmappedAlignmentNum 
		<< "\tperc: " << ((double)unmappedAlignmentNum/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "generated read pair number: " << unmappedAlignmentNum/2 << endl;
	fclose(fp);
	//output_log_ofs.close();
	output_read_end1_ofs.close();
	output_read_end2_ofs.close();
	return;	
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputSam outputFolder" << endl;
		exit(1);
	}
	string inputSam = argv[1];
	cout << "creating folder ......" << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

	log_ofs << "start to initiate outputfiles ..." << endl;
	string outputFilePrefix = outputFolderStr + "/unmappedRead";
	string outputUnmappedReadFq_1 = outputFilePrefix + ".1.fq";
	string outputUnmappedReadFq_2 = outputFilePrefix + ".2.fq";
	unmappedAlignemnt2ReadFile(inputSam,
		outputUnmappedReadFq_1, outputUnmappedReadFq_2, log_ofs);

	log_ofs << "job ends ...." << endl;
	log_ofs.close();
	return 0;
}