#ifndef KMER2GROUPEDKMER_INFO_H
#define KMER2GROUPEDKMER_INFO_H
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

class Kmer2groupedKmer_Info
{
private:
	string inputRawKmerFile;
	string outputGroupedKmerDir;
	int Kmer_prefix_length;

	vector<string> outputGroupedKmerFileVec;
	int group_num;
public:
	Kmer2groupedKmer_Info()
	{}

	string returnTmpBinFile(int tmpBin)
	{
		return outputGroupedKmerFileVec[tmpBin];
	}

	void initiate_mkdir(string& tmpInputRawKmerFile, 
		string& tmpOutputGroupedKmerDir, int tmp_Kmer_prefix_length)
	{
		inputRawKmerFile = tmpInputRawKmerFile;
		outputGroupedKmerDir = tmpOutputGroupedKmerDir + "/";
		Kmer_prefix_length = tmp_Kmer_prefix_length;

		string cmd_mkdir = "mkdir " + outputGroupedKmerDir;
		system(cmd_mkdir.c_str());

		group_num = pow(4, Kmer_prefix_length);

		for(int tmp = 0; tmp < group_num; tmp++)
		{
			string tmpBinKmerFile = outputGroupedKmerDir + "bin." + int_to_str(tmp) + ".Kmer";
			outputGroupedKmerFileVec.push_back(tmpBinKmerFile);
		}
	}

	void groupKmer()
	{
		// start to initiate bin files
		vector<ofstream*> KmerBinOfsVec;
		for(int tmp = 0; tmp < group_num; tmp ++)
		{
			string tmpBinFile = outputGroupedKmerFileVec[tmp];
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
			string prefix3charStr = tmpStr.substr(0, Kmer_prefix_length);
			int tmpStrPrefixIndex = this->KmerPrefix2int(prefix3charStr);
			(*KmerBinOfsVec[tmpStrPrefixIndex]) << tmpStr << endl;
		}
		// release and close files
		for(int tmp = 0; tmp < group_num; tmp ++)
		{
			(*KmerBinOfsVec[tmp]).close();
			delete KmerBinOfsVec[tmp];			
		}
	}

	int KmerPrefix2int(string& KmerPrefixStr)
	{
		int tmpLength = KmerPrefixStr.length();
		if(tmpLength != 3)
		{
			cout << "error in KmerPrefix2int, tmpLength != 3" << endl;
			exit(1);
		}
		else
			return (this->tripleChar2int(KmerPrefixStr));
	}

	int tripleChar2int(string& tripleCharStr)
	{
		if(tripleCharStr.at(0) == 'A')
		{
			if(tripleCharStr.at(1) == 'A')
			{
				if(tripleCharStr.at(2) == 'A')
					return 0;
				else if(tripleCharStr.at(2) == 'C')
					return 1;
				else if(tripleCharStr.at(2) == 'G')
					return 2;
				else if(tripleCharStr.at(2) == 'T')
					return 3;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'C')
			{
				if(tripleCharStr.at(2) == 'A')
					return 4;
				else if(tripleCharStr.at(2) == 'C')
					return 5;
				else if(tripleCharStr.at(2) == 'G')
					return 6;
				else if(tripleCharStr.at(2) == 'T')
					return 7;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'G')
			{
				if(tripleCharStr.at(2) == 'A')
					return 8;
				else if(tripleCharStr.at(2) == 'C')
					return 9;
				else if(tripleCharStr.at(2) == 'G')
					return 10;
				else if(tripleCharStr.at(2) == 'T')
					return 11;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'T')
			{
				if(tripleCharStr.at(2) == 'A')
					return 12;
				else if(tripleCharStr.at(2) == 'C')
					return 13;
				else if(tripleCharStr.at(2) == 'G')
					return 14;
				else if(tripleCharStr.at(2) == 'T')
					return 15;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else
			{
				cout << "error in tripleChar2int, tripleCharStr.at(1): " << tripleCharStr.at(1) << endl;
				exit(1);
			}			
		}
		else if(tripleCharStr.at(0) == 'C')
		{
			if(tripleCharStr.at(1) == 'A')
			{
				if(tripleCharStr.at(2) == 'A')
					return 16;
				else if(tripleCharStr.at(2) == 'C')
					return 17;
				else if(tripleCharStr.at(2) == 'G')
					return 18;
				else if(tripleCharStr.at(2) == 'T')
					return 19;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'C')
			{
				if(tripleCharStr.at(2) == 'A')
					return 20;
				else if(tripleCharStr.at(2) == 'C')
					return 21;
				else if(tripleCharStr.at(2) == 'G')
					return 22;
				else if(tripleCharStr.at(2) == 'T')
					return 23;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'G')
			{
				if(tripleCharStr.at(2) == 'A')
					return 24;
				else if(tripleCharStr.at(2) == 'C')
					return 25;
				else if(tripleCharStr.at(2) == 'G')
					return 26;
				else if(tripleCharStr.at(2) == 'T')
					return 27;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'T')
			{
				if(tripleCharStr.at(2) == 'A')
					return 28;
				else if(tripleCharStr.at(2) == 'C')
					return 29;
				else if(tripleCharStr.at(2) == 'G')
					return 30;
				else if(tripleCharStr.at(2) == 'T')
					return 31;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else
			{
				cout << "error in tripleChar2int, tripleCharStr.at(1): " << tripleCharStr.at(1) << endl;
				exit(1);
			}			
		}
		else if(tripleCharStr.at(0) == 'G')
		{
			if(tripleCharStr.at(1) == 'A')
			{
				if(tripleCharStr.at(2) == 'A')
					return 32;
				else if(tripleCharStr.at(2) == 'C')
					return 33;
				else if(tripleCharStr.at(2) == 'G')
					return 34;
				else if(tripleCharStr.at(2) == 'T')
					return 35;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'C')
			{
				if(tripleCharStr.at(2) == 'A')
					return 36;
				else if(tripleCharStr.at(2) == 'C')
					return 37;
				else if(tripleCharStr.at(2) == 'G')
					return 38;
				else if(tripleCharStr.at(2) == 'T')
					return 39;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'G')
			{
				if(tripleCharStr.at(2) == 'A')
					return 40;
				else if(tripleCharStr.at(2) == 'C')
					return 41;
				else if(tripleCharStr.at(2) == 'G')
					return 42;
				else if(tripleCharStr.at(2) == 'T')
					return 43;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'T')
			{
				if(tripleCharStr.at(2) == 'A')
					return 44;
				else if(tripleCharStr.at(2) == 'C')
					return 45;
				else if(tripleCharStr.at(2) == 'G')
					return 46;
				else if(tripleCharStr.at(2) == 'T')
					return 47;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else
			{
				cout << "error in tripleChar2int, tripleCharStr.at(1): " << tripleCharStr.at(1) << endl;
				exit(1);
			}			
		}
		else if(tripleCharStr.at(0) == 'T')
		{
			if(tripleCharStr.at(1) == 'A')
			{
				if(tripleCharStr.at(2) == 'A')
					return 48;
				else if(tripleCharStr.at(2) == 'C')
					return 49;
				else if(tripleCharStr.at(2) == 'G')
					return 50;
				else if(tripleCharStr.at(2) == 'T')
					return 51;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'C')
			{
				if(tripleCharStr.at(2) == 'A')
					return 52;
				else if(tripleCharStr.at(2) == 'C')
					return 53;
				else if(tripleCharStr.at(2) == 'G')
					return 54;
				else if(tripleCharStr.at(2) == 'T')
					return 55;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'G')
			{
				if(tripleCharStr.at(2) == 'A')
					return 56;
				else if(tripleCharStr.at(2) == 'C')
					return 57;
				else if(tripleCharStr.at(2) == 'G')
					return 58;
				else if(tripleCharStr.at(2) == 'T')
					return 59;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else if(tripleCharStr.at(1) == 'T')
			{
				if(tripleCharStr.at(2) == 'A')
					return 60;
				else if(tripleCharStr.at(2) == 'C')
					return 61;
				else if(tripleCharStr.at(2) == 'G')
					return 62;
				else if(tripleCharStr.at(2) == 'T')
					return 63;
				else
				{
					cout << "error in tripleChar2int, tripleCharStr.at(2): " << tripleCharStr.at(2) << endl;
					exit(1);
				}	
			}
			else
			{
				cout << "error in tripleChar2int, tripleCharStr.at(1): " << tripleCharStr.at(1) << endl;
				exit(1);
			}			
		}
		else
		{
			cout << "error in tripleChar2int, tripleCharStr.at(0): " << tripleCharStr.at(0) << endl;
			exit(1);
		}
	}

	int tripleChar2int_old(string& tripleCharStr)
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
};
#endif