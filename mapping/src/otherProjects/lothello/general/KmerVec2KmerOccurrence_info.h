// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef KMERVEC2KMEROCCURRENCE_INFO_H
#define KMERVEC2KMEROCCURRENCE_INFO_H
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

using namespace std;

class KmerVec2KmerOccurrence_Info
{
private:
	int total_Kmer_length; 
	int grouped_Kmer_prefix_length;

public:
	KmerVec2KmerOccurrence_Info()
	{}

	void initiate(int tmp_total_Kmer_length, int tmp_grouped_Kmer_prefix_length)
	{
		total_Kmer_length = tmp_total_Kmer_length;
		grouped_Kmer_prefix_length = tmp_grouped_Kmer_prefix_length;
	}

	void printoutKmerOccurrenceNum(char* toKeepKmerCount_charArray, 
		string& tmpOutputKmerOccurrenceNumFile, int tmpBin, bool outputTotalKmerOrSolidKmer_bool)
	{
		ofstream KmerOccurrenceNum_ofs(tmpOutputKmerOccurrenceNumFile.c_str());
		int toKeepCount_Kmer_length = total_Kmer_length - grouped_Kmer_prefix_length;
		unsigned long long int toKeepCount_bin_Kmer_num = pow(4, toKeepCount_Kmer_length);
		string KmerPrefix = this->int2Kmer((unsigned long long)tmpBin, grouped_Kmer_prefix_length);
		if(outputTotalKmerOrSolidKmer_bool)
		{
			for(unsigned long long int tmp = 0; tmp < toKeepCount_bin_Kmer_num; tmp++)
			{
				string tmpToKeepCountKmer = this->int2Kmer(tmp, toKeepCount_Kmer_length);
				string tmpRawKmer = KmerPrefix + tmpToKeepCountKmer;
				int tmpKmerCount = toKeepKmerCount_charArray[tmp];
				KmerOccurrenceNum_ofs << tmpRawKmer << "\t" << tmpKmerCount << endl;
			}
		}
		else
		{
			unsigned long long int tmpSolidKmerNum = 0;
			for(unsigned long long int tmp = 0; tmp < toKeepCount_bin_Kmer_num; tmp++)
			{
				int tmpKmerCount = (int)toKeepKmerCount_charArray[tmp];
				if(tmpKmerCount > 0)
				{
					//cout << "tmp: " << tmp << endl;
					//cout << "tmpKmerCount: " << tmpKmerCount << endl;
					tmpSolidKmerNum ++;
					string tmpToKeepCountKmer = this->int2Kmer(tmp, toKeepCount_Kmer_length);
					string tmpRawKmer = KmerPrefix + tmpToKeepCountKmer;
					KmerOccurrenceNum_ofs << tmpRawKmer << "\t" << tmpKmerCount << endl;
				}
			}
			cout << "tmpSolidKmerNum: " << tmpSolidKmerNum << endl;
		}
		KmerOccurrenceNum_ofs.close();
	}

	void printoutSetSpecificKmer(char* toKeepKmerCount_charArray, string& toAddId,
		string& tmpInputRawKmerFile, string& tmpOutputSetSpecificKmerFile)
	{
		ifstream rawKmer_ifs(tmpInputRawKmerFile.c_str());
		ofstream setSpecificKmer_ofs(tmpOutputSetSpecificKmerFile.c_str());
		while(!rawKmer_ifs.eof())
		{
			string tmpStr;
			getline(rawKmer_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpKmer = tmpStr.substr(0, tabLoc);
			unsigned long long int tmpKmerInt = this->toKeepCount_Kmer2int(tmpKmer);
			if(toKeepKmerCount_charArray[tmpKmerInt] == 1)
				setSpecificKmer_ofs << tmpKmer << "\t" << toAddId << endl;
		}
		setSpecificKmer_ofs.close();
		rawKmer_ifs.close();
	}

	string int2Kmer(unsigned long long tmpInt, int tmp_Kmer_length)
	{
		string tmpStr = "";
		for(int tmpBase = tmp_Kmer_length - 1; tmpBase >= 0; tmpBase --)
		{
			unsigned long long toAndBaseNum = 3;
			toAndBaseNum <<=(2*tmpBase);
			unsigned long long toAndResults = tmpInt & toAndBaseNum;
			toAndResults >>=(2*tmpBase);
			if(toAndResults == 0)
				tmpStr += "A";
			else if(toAndResults == 1)
				tmpStr += "C";
			else if(toAndResults == 2)
				tmpStr += "G";
			else if(toAndResults == 3)
				tmpStr += "T";
			else 
			{
				cout << "error! toAndResults: " << toAndResults << endl;
				exit(1);
			}					
		}
		return tmpStr;
	}	

	void updateKmerOccurrenceNum(char* toKeepKmerCount_charArray, string& tmpKmerFile)
	{
		ifstream tmpKmer_ifs(tmpKmerFile.c_str());
		while(!tmpKmer_ifs.eof())
		{
			string tmpStr;
			getline(tmpKmer_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpKmerStr = tmpStr.substr(0, tabLoc);
			if(tmpKmerStr.length() != total_Kmer_length)
			{
				cout << "error! (tmpKmerStr.length() != total_Kmer_length)" << endl;
				exit(1);
			}
			unsigned long long tmpKmerInt = this->toKeepCount_Kmer2int(tmpKmerStr);
			if(toKeepKmerCount_charArray[tmpKmerInt] == 255)
			{}
			else
				toKeepKmerCount_charArray[tmpKmerInt] ++;
		}
		tmpKmer_ifs.eof();
	}

	unsigned long long toKeepCount_Kmer2int(string& tmpRawKmer)
	{
		unsigned long long tmpInt = 0;
		for(int tmp = grouped_Kmer_prefix_length; tmp < total_Kmer_length; tmp++)
		{
			tmpInt <<=2;
			switch(tmpRawKmer.at(tmp))
			{
				case 'A':
					//tmpInt ++;
					break;
				case 'C':
					tmpInt ++;
					break;
				case 'G':
					tmpInt +=2;
					break;
				case 'T':
					tmpInt +=3;
					break;
				default:
					cout << "invalid tmpRawKmer.at(tmp): " << tmpRawKmer.at(tmp) << endl;
					exit(1);											
			}
		}
		return tmpInt;
	}
};
#endif