// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef KMER2SORTEDKMER_INFO_H
#define KMER2SORTEDKMER_INFO_H
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

class Kmer2sortedKmer_Info
{
private:
	string inputRawKmerFile;
	string outputSortedKmerFile;
	string temporaryDir;
public:
	Kmer2sortedKmer_Info()
	{}

	void initiate_mkTmpDir(string& tmpInputRawKmerFile, string& tmpOutputSortedKmerFile, string& tmpTemporaryDir)
	{
		inputRawKmerFile = tmpInputRawKmerFile;
		outputSortedKmerFile = tmpOutputSortedKmerFile;
		temporaryDir = tmpTemporaryDir;
		string cmd_mkdir_temporaryDir = "mkdir " + temporaryDir;
		system(cmd_mkdir_temporaryDir.c_str());
	}

	void initiate_noMkTmpDir(string& tmpInputRawKmerFile, string& tmpOutputSortedKmerFile, string& tmpTemporaryDir)
	{
		inputRawKmerFile = tmpInputRawKmerFile;
		outputSortedKmerFile = tmpOutputSortedKmerFile;
		temporaryDir = tmpTemporaryDir;
		//string cmd_mkdir_temporaryDir = "mkdir " + temporaryDir;
		//system(cmd_mkdir_temporaryDir.c_str());
	}	

	void rm_tmpDir()
	{
		string cmd_rm_tmpDir = "rm -r " + temporaryDir;
		system(cmd_rm_tmpDir.c_str());
	}

	void sort_Kmer_keepInterDir()
	{
		this->sort_Kmer_keepInterDir(1);
	}

	void sort_Kmer_keepInterDir(int threads_num)
	{
		int bin_num = 64;
		string output_bin_dir = temporaryDir + "/bin/";
		string cmd_mkdir_output_bin_dir = "mkdir " + output_bin_dir;
		system(cmd_mkdir_output_bin_dir.c_str());

		string output_tmpSortDir = temporaryDir + "/tmpSortDir/";
		string cmd_mkdir_output_tmpSortDir = "mkdir " + output_tmpSortDir;
		system(cmd_mkdir_output_tmpSortDir.c_str());
		vector<string> output_tmpSortDir_threadVec;
		for(int tmp = 0; tmp < threads_num; tmp++)
		{
			string output_tmpSortDir_tmpThread = temporaryDir + "/tmpSortDir/thread_" + int_to_str(tmp) + "/";
			output_tmpSortDir_threadVec.push_back(output_tmpSortDir_tmpThread);
			string cmd_mkdir_output_tmpSortDir_tmpThread = "mkdir " + output_tmpSortDir_tmpThread;
			system(cmd_mkdir_output_tmpSortDir_tmpThread.c_str());
		}
		
		string output_sortedKmerBin = temporaryDir + "/sortedKmerBin/";
		string cmd_mkdir_output_sortedKmerBin = "mkdir " + output_sortedKmerBin;
		system(cmd_mkdir_output_sortedKmerBin.c_str());
		// start to initiate bin files
		vector<ofstream*> KmerBinOfsVec;
		for(int tmp = 0; tmp < bin_num; tmp ++)
		{
			string tmpBinFile = output_bin_dir + "bin." + int_to_str(tmp) + ".Kmer";
			ofstream *tmpBin_ofs = new ofstream(tmpBinFile.c_str());
			KmerBinOfsVec.push_back(tmpBin_ofs);			
		}
		// start to divide Kmers into bins
		ifstream Kmer_ifs(inputRawKmerFile.c_str());
		while(!Kmer_ifs.eof())
		{
			string tmpStr;
			getline(Kmer_ifs, tmpStr);
			if(tmpStr == "")
				break;
			string prefix3charStr = tmpStr.substr(0, 3);
			int tmpStrPrefixIndex = this->tripleChar2int(prefix3charStr);
			(*KmerBinOfsVec[tmpStrPrefixIndex]) << tmpStr << endl;
		}
		// release and close files
		for(int tmp = 0; tmp < bin_num; tmp ++)
		{
			(*KmerBinOfsVec[tmp]).close();
			delete KmerBinOfsVec[tmp];			
		}
		// start to sort each binKmer file
		omp_set_num_threads(threads_num);
		for(int tmp = 0; tmp < bin_num; tmp ++)
		{
			int tmpThreadNO = omp_get_thread_num();
			string tmpSourceKmerFile = output_bin_dir + "bin." + int_to_str(tmp) + ".Kmer";
			string tmpDestSortedKmerFile = output_sortedKmerBin + "bin." + int_to_str(tmp) + ".sortedKmer";
			string tmp_cmd_sortKmer = "sort -k1 -T " + output_tmpSortDir_threadVec[tmpThreadNO] + "/ " + tmpSourceKmerFile
				+ " > " + tmpDestSortedKmerFile;
			system(tmp_cmd_sortKmer.c_str());
		}
		// cat all sorted binKmer files
		string cmd_cat_sortedKmer = "cat";// + output_sortedKmerBin + "/*.sortedKmer > " + outputSortedKmerFile;
		for(int tmp = 0; tmp < bin_num; tmp++)
		{
			cmd_cat_sortedKmer += " ";
			string tmpSortedKmerFile = output_sortedKmerBin + "/bin." + int_to_str(tmp) + ".sortedKmer";
			cmd_cat_sortedKmer += tmpSortedKmerFile;
		}
		cmd_cat_sortedKmer += " > ";
		cmd_cat_sortedKmer += outputSortedKmerFile;
		system(cmd_cat_sortedKmer.c_str());
		// remove intermediate files
		//string cmd_rm_output_bin_dir = "rm -r " + output_bin_dir;
		//system(cmd_rm_output_bin_dir.c_str());
		//string cmd_rm_tmpSortDir = "rm -r " + output_tmpSortDir;
		//system(cmd_rm_tmpSortDir.c_str());
		//string cmd_rm_output_sortedKmerBin = "rm -r " + output_sortedKmerBin;
		//system(cmd_rm_output_sortedKmerBin.c_str());
	}

	void sort_Kmer_rmInterDir()
	{
		this->sort_Kmer_rmInterDir(1);
	}

	void sort_Kmer_rmInterDir(int threads_num)
	{
		this->sort_Kmer_keepInterDir(threads_num);
		// remove intermediate files
		string output_bin_dir = temporaryDir + "/bin/";
		string cmd_rm_output_bin_dir = "rm -r " + output_bin_dir;
		system(cmd_rm_output_bin_dir.c_str());
		string output_tmpSortDir = temporaryDir + "/tmpSortDir/";
		string cmd_rm_tmpSortDir = "rm -r " + output_tmpSortDir;
		system(cmd_rm_tmpSortDir.c_str());
		string output_sortedKmerBin = temporaryDir + "/sortedKmerBin/";
		string cmd_rm_output_sortedKmerBin = "rm -r " + output_sortedKmerBin;
		system(cmd_rm_output_sortedKmerBin.c_str());
	}

	int tripleChar2int(string& tripleCharStr)
	{
		if(tripleCharStr.length() != 3)
		{
			cout << "tripleCharStr.length() != 3 in tripleChar2int" << endl;
			exit(1); 
		}
		if(tripleCharStr == "AAA")
			return 0;
		else if(tripleCharStr == "AAC")
			return 1;
		else if(tripleCharStr == "AAG")
			return 2;
		else if(tripleCharStr == "AAT")
			return 3;	
		else if(tripleCharStr == "ACA")
			return 4;
		else if(tripleCharStr == "ACC")
			return 5;
		else if(tripleCharStr == "ACG")
			return 6;					
		else if(tripleCharStr == "ACT")
			return 7;		
		else if(tripleCharStr == "AGA")
			return 8;
		else if(tripleCharStr == "AGC")
			return 9;
		else if(tripleCharStr == "AGG")
			return 10;					
		else if(tripleCharStr == "AGT")
			return 11;	
		else if(tripleCharStr == "ATA")
			return 12;
		else if(tripleCharStr == "ATC")
			return 13;
		else if(tripleCharStr == "ATG")
			return 14;	
		else if(tripleCharStr == "ATT")
			return 15;
		else if(tripleCharStr == "CAA")
			return 16;
		else if(tripleCharStr == "CAC")
			return 17;
		else if(tripleCharStr == "CAG")
			return 18;
		else if(tripleCharStr == "CAT")
			return 19;	
		else if(tripleCharStr == "CCA")
			return 20;
		else if(tripleCharStr == "CCC")
			return 21;
		else if(tripleCharStr == "CCG")
			return 22;					
		else if(tripleCharStr == "CCT")
			return 23;		
		else if(tripleCharStr == "CGA")
			return 24;
		else if(tripleCharStr == "CGC")
			return 25;
		else if(tripleCharStr == "CGG")
			return 26;					
		else if(tripleCharStr == "CGT")
			return 27;	
		else if(tripleCharStr == "CTA")
			return 28;
		else if(tripleCharStr == "CTC")
			return 29;
		else if(tripleCharStr == "CTG")
			return 30;	
		else if(tripleCharStr == "CTT")
			return 31;
		else if(tripleCharStr == "GAA")
			return 32;
		else if(tripleCharStr == "GAC")
			return 33;
		else if(tripleCharStr == "GAG")
			return 34;
		else if(tripleCharStr == "GAT")
			return 35;	
		else if(tripleCharStr == "GCA")
			return 36;
		else if(tripleCharStr == "GCC")
			return 37;
		else if(tripleCharStr == "GCG")
			return 38;					
		else if(tripleCharStr == "GCT")
			return 39;		
		else if(tripleCharStr == "GGA")
			return 40;
		else if(tripleCharStr == "GGC")
			return 41;
		else if(tripleCharStr == "GGG")
			return 42;					
		else if(tripleCharStr == "GGT")
			return 43;	
		else if(tripleCharStr == "GTA")
			return 44;
		else if(tripleCharStr == "GTC")
			return 45;
		else if(tripleCharStr == "GTG")
			return 46;	
		else if(tripleCharStr == "GTT")
			return 47;		
		else if(tripleCharStr == "TAA")
			return 48;
		else if(tripleCharStr == "TAC")
			return 49;
		else if(tripleCharStr == "TAG")
			return 50;
		else if(tripleCharStr == "TAT")
			return 51;	
		else if(tripleCharStr == "TCA")
			return 52;
		else if(tripleCharStr == "TCC")
			return 53;
		else if(tripleCharStr == "TCG")
			return 54;					
		else if(tripleCharStr == "TCT")
			return 55;		
		else if(tripleCharStr == "TGA")
			return 56;
		else if(tripleCharStr == "TGC")
			return 57;
		else if(tripleCharStr == "TGG")
			return 58;					
		else if(tripleCharStr == "TGT")
			return 59;	
		else if(tripleCharStr == "TTA")
			return 60;
		else if(tripleCharStr == "TTC")
			return 61;
		else if(tripleCharStr == "TTG")
			return 62;	
		else if(tripleCharStr == "TTT")
			return 63;
		else
		{
			cout << "invalid char in tripleChar2int: " << tripleCharStr << endl;
			exit(1);
		}
	}

	string return_outputSortedKmerFile()
	{
		return outputSortedKmerFile;
	}
};
#endif