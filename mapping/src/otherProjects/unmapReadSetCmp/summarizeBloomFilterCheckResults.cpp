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
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputResultsFilePrefix dataNum outputFolder" << endl;
		exit(1);
	}
	string inputResultsFilePrefix = argv[1];
	string dataNumStr = argv[2];
	int dataNum = atoi(dataNumStr.c_str());
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	string output_uniqueKmer_summarize_file = outputFolderStr + "uniqueKmer.summarize.txt";
	string output_frequencyKmer_summarize_file = outputFolderStr + "frequencyKmer.summarize.txt";
	ofstream uniqueKmer_summarize_ofs(output_uniqueKmer_summarize_file.c_str());
	ofstream frequencyKmer_summarize_ofs(output_frequencyKmer_summarize_file.c_str());
	log_ofs << "start to examine bloomFilterCheckResults" << endl;
	cout << "start to examine bloomFilterCheckResults" << endl;
	vector< vector< pair<int,int> > > uniqueKmerNumPairVecVec; //< hitNumUnique, missNumUnique>
	vector< vector< pair<int,int> > > freqKmerNumPairVecVec; //< hitNumFreq, missNumFreq>
	for(int tmp1 = 1; tmp1 <= dataNum; tmp1++)
	{
		vector< pair<int,int> > tmpUniqueKmerNumPairVec;
		vector< pair<int,int> > tmpFreqKmerNumPairVec;
		for(int tmp2 = 1; tmp2 <= dataNum; tmp2++)
		{
			tmpUniqueKmerNumPairVec.push_back(pair<int,int>(-1,-1));
			tmpFreqKmerNumPairVec.push_back(pair<int,int>(-1,-1));
		}
		uniqueKmerNumPairVecVec.push_back(tmpUniqueKmerNumPairVec);
		freqKmerNumPairVecVec.push_back(tmpFreqKmerNumPairVec);
	}			
	for(int tmp1 = 1; tmp1 <= dataNum; tmp1++)
	{
		for(int tmp2 = 1; tmp2 <= dataNum; tmp2++)
		{
			cout << "tmp1: " << tmp1 << endl;
			cout << "tmp2: " << tmp2 << endl;
			string tmpFolderPath = inputResultsFilePrefix + int_to_str(tmp1) + "/" + int_to_str(tmp2) + "/";
			string tmp_log = tmpFolderPath + "log.txt";
			//cout << "tmp_log: " << endl << tmp_log << endl;
			ifstream tmp_log_ifs(tmp_log.c_str());
			if(!tmp_log_ifs.good())
				continue;
			vector<string> tmpLogStrVec;
			while(!tmp_log_ifs.eof())
			{
				string tmpStr;
				getline(tmp_log_ifs, tmpStr);
				//if(tmpStr == "")
				//	break;
				tmpLogStrVec.push_back(tmpStr);
			}
			//cout << "tmpLogStrVec[11]: " << endl << tmpLogStrVec[11] << endl;
			if(tmpLogStrVec.size() < 12)
				continue;
			else if(tmpLogStrVec[11] != "end of query")
				continue;
			else
			{}
			tmp_log_ifs.close();
			string tmp_miss_query_file = tmpFolderPath + "miss_query.txt";
			string tmp_hit_query_file = tmpFolderPath + "hit_query.txt";
			string tmp_summarize_file = tmpFolderPath + "summarize.txt";
			ifstream tmp_miss_ifs(tmp_miss_query_file.c_str());
			ifstream tmp_hit_ifs(tmp_hit_query_file.c_str());
			ofstream tmp_summarize_ofs(tmp_summarize_file.c_str());
			int tmp_uniqueKmerCount_hit = 0, tmp_uniqueKmerCount_miss = 0,
				tmp_freqKmerCount_hit = 0, tmp_freqKmerCount_miss = 0;
			while(!tmp_hit_ifs.eof())
			{
				string tmpHitStr_name;
				getline(tmp_hit_ifs, tmpHitStr_name);
				if(tmpHitStr_name == "")
					break;
				string tmpHitStr_seq;
				getline(tmp_hit_ifs, tmpHitStr_seq);
				tmp_uniqueKmerCount_hit ++;
				string tmpHitStr_freq = tmpHitStr_name.substr(1);
				int tmpHitInt_freq = atoi(tmpHitStr_freq.c_str());
				tmp_freqKmerCount_hit += tmpHitInt_freq;
			}
			while(!tmp_miss_ifs.eof())
			{
				string tmpMissStr_name;
				getline(tmp_miss_ifs, tmpMissStr_name);
				if(tmpMissStr_name == "")
					break;
				string tmpMissStr_seq;
				getline(tmp_miss_ifs, tmpMissStr_seq);
				tmp_uniqueKmerCount_miss ++;
				string tmpMissStr_freq = tmpMissStr_name.substr(1);
				int tmpMissInt_freq = atoi(tmpMissStr_freq.c_str());
				tmp_freqKmerCount_miss += tmpMissInt_freq;
			}
			(uniqueKmerNumPairVecVec[tmp1-1])[tmp2-1].first = tmp_uniqueKmerCount_hit;
			(uniqueKmerNumPairVecVec[tmp1-1])[tmp2-1].second = tmp_uniqueKmerCount_miss;
			(freqKmerNumPairVecVec[tmp1-1])[tmp2-1].first = tmp_freqKmerCount_hit;
			(freqKmerNumPairVecVec[tmp1-1])[tmp2-1].second = tmp_freqKmerCount_miss;
			tmp_summarize_ofs << "unique_kmer_hit: " << endl << tmp_uniqueKmerCount_hit 
				<< endl << "unique_kmer_miss: " << endl << tmp_uniqueKmerCount_miss
				<< endl << "frequency_kmer_hit: " << endl << tmp_freqKmerCount_hit
				<< endl << "frequency_kmer_miss: " << endl << tmp_freqKmerCount_miss << endl;
			tmp_summarize_ofs.close();
			tmp_hit_ifs.close();
			tmp_miss_ifs.close();
		}
	}
	for(int tmp1 = 1; tmp1 <= dataNum; tmp1++)
	{
		for(int tmp2 = 1; tmp2 <= dataNum; tmp2++)
		{
			int tmp_uniqueKmerCount_hit = ((uniqueKmerNumPairVecVec[tmp1-1])[tmp2-1]).first;
			int tmp_uniqueKmerCount_miss = ((uniqueKmerNumPairVecVec[tmp1-1])[tmp2-1]).second;
			int tmp_freqKmerCount_hit = ((freqKmerNumPairVecVec[tmp1-1])[tmp2-1]).first;
			int tmp_freqKmerCount_miss = ((freqKmerNumPairVecVec[tmp1-1])[tmp2-1]).second;
			if(tmp_uniqueKmerCount_hit < 0)
			{
				if(tmp2 == 1)
				{
					uniqueKmer_summarize_ofs << "N";
					frequencyKmer_summarize_ofs << "N";
				}
				else if(tmp2 == dataNum)
				{
					uniqueKmer_summarize_ofs << "\tN" << endl;
					frequencyKmer_summarize_ofs << "\tN" << endl;				
				}
				else
				{
					uniqueKmer_summarize_ofs << "\tN";
					frequencyKmer_summarize_ofs << "\tN";			
				}				
			}
			else
			{
				double tmp_uniqueKmer_hitPerc = (double)tmp_uniqueKmerCount_hit/(double)(tmp_uniqueKmerCount_hit + tmp_uniqueKmerCount_miss);
				double tmp_freqKmer_hitPerc = (double)tmp_freqKmerCount_hit/(double)(tmp_freqKmerCount_hit + tmp_freqKmerCount_miss);
				if(tmp2 == 1)
				{
					uniqueKmer_summarize_ofs << tmp_uniqueKmer_hitPerc;
					frequencyKmer_summarize_ofs << tmp_freqKmer_hitPerc;
				}
				else if(tmp2 == dataNum)
				{
					uniqueKmer_summarize_ofs << "\t" << tmp_uniqueKmer_hitPerc << endl;
					frequencyKmer_summarize_ofs << "\t" << tmp_freqKmer_hitPerc << endl;				
				}
				else
				{
					uniqueKmer_summarize_ofs << "\t" << tmp_uniqueKmer_hitPerc;
					frequencyKmer_summarize_ofs << "\t" << tmp_freqKmer_hitPerc;				
				}
			}
		}
	}		

	log_ofs.close();
	uniqueKmer_summarize_ofs.close();
	frequencyKmer_summarize_ofs.close();
	return 0;
}