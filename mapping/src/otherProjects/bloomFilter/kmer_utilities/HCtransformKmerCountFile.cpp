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

using namespace std;

string returnTransformedKmerSeq(string& tmpRawKmerSeq)
{
	string tmpTransformedKmerSeq = "";
	int tmpRawKmerSeqLen = tmpRawKmerSeq.length();
	string lastBaseCharStr = "X";
	for(int tmp = 0; tmp < tmpRawKmerSeqLen; tmp++)
	{
		string tmpCharStr = tmpRawKmerSeq.substr(tmp, 1);
		if(tmpCharStr != lastBaseCharStr)
		{
			tmpTransformedKmerSeq += tmpCharStr;
			lastBaseCharStr = tmpCharStr;
		}
	}
	return tmpTransformedKmerSeq;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputRawKmerCountFile outputTransformedKmerCountFile" << endl;
		exit(1);
	}
	bool remove_intermediate_dir_file_bool = true;

	string inputRawKmerCountFile = argv[1];
	string outputTransformedKmerCountFile = argv[2];
	string outputRawKmerCountFile_transformed = outputTransformedKmerCountFile + ".tmp.transformed";
	string outputRawKmerCountFile_transformed_sorted = outputRawKmerCountFile_transformed + ".sorted";
	string tmp_sort_dir = outputRawKmerCountFile_transformed + ".sort.tmpDir";
	
	// start to transform each Kmer in raw kmer count file
	ifstream rawKmerCount_ifs(inputRawKmerCountFile.c_str());
	ofstream rawKmerCount_transformed_ofs(outputRawKmerCountFile_transformed.c_str());
	while(!rawKmerCount_ifs.eof())
	{
		string tmpStr;
		getline(rawKmerCount_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpRawKmerSeq;
		string tmpRawKmerCountStr;
		string tmpTransformedKmerSeq;
		if(tabLoc == string::npos) // no other field, only k-mer exists
		{
			cout << "error ! count field does not exist !" << endl;
			exit(1);
		}
		else
		{
			tmpRawKmerSeq = tmpStr.substr(0, tabLoc);
			tmpTransformedKmerSeq = returnTransformedKmerSeq(tmpRawKmerSeq);
			tmpRawKmerCountStr = tmpStr.substr(tabLoc + 1);
			rawKmerCount_transformed_ofs << tmpTransformedKmerSeq << "\t" << tmpRawKmerCountStr << endl;
		}
	}
	rawKmerCount_transformed_ofs.close();
	rawKmerCount_ifs.close();

	// start to sort transformed kmer count file
	string mkdir_cmd = "mkdir " + tmp_sort_dir;
	system(mkdir_cmd.c_str());
	string sort_cmd = "sort -k1 -T " + tmp_sort_dir + " " + outputRawKmerCountFile_transformed
		+ " > " + outputRawKmerCountFile_transformed_sorted;
	system(sort_cmd.c_str());

	// merge duplicate transformed K-mers
	ifstream rawKmerCountFile_transformed_sorted_ifs(outputRawKmerCountFile_transformed_sorted.c_str());
	ofstream transformedKmerCount_ofs(outputTransformedKmerCountFile.c_str());
	string lastTransformedKmer_uniq_infoStr;
	string lastTransformedKmer_uniq;
	int lastTransformedKmer_uniq_count;
	getline(rawKmerCountFile_transformed_sorted_ifs, lastTransformedKmer_uniq_infoStr);
	if(lastTransformedKmer_uniq_infoStr == "")
	{
		cout << "error ! invalid file, the 1st line is NULL" << endl;
		exit(1);
	}
	else
	{	
		int tabLoc = lastTransformedKmer_uniq_infoStr.find("\t");
		lastTransformedKmer_uniq = lastTransformedKmer_uniq_infoStr.substr(0, tabLoc);
		string lastTransformedKmer_uniq_countStr = lastTransformedKmer_uniq_infoStr.substr(tabLoc + 1);
		lastTransformedKmer_uniq_count = atoi(lastTransformedKmer_uniq_countStr.c_str());
	}
	while(!rawKmerCountFile_transformed_sorted_ifs.eof())
	{
		string tmpStr;
		getline(rawKmerCountFile_transformed_sorted_ifs, tmpStr);
		if(tmpStr == "")
		{
			transformedKmerCount_ofs << lastTransformedKmer_uniq 
				<< "\t" << lastTransformedKmer_uniq_count << endl;			
			break;
		}
		int tabLoc = tmpStr.find("\t");
		string tmpTransformedKmerSeq;
		string tmpTransformedKmerCountStr;
		if(tabLoc == string::npos) // no other field, only k-mer exists
		{
			cout << "error ! count field does not exist !" << endl;
			exit(1);
		}
		else
		{
			tmpTransformedKmerSeq = tmpStr.substr(0, tabLoc);
			tmpTransformedKmerCountStr = tmpStr.substr(tabLoc + 1);
			int tmpTransformedKmerCount = atoi(tmpTransformedKmerCountStr.c_str());
			if(tmpTransformedKmerSeq == lastTransformedKmer_uniq)
			{
				lastTransformedKmer_uniq_count += tmpTransformedKmerCount;
			}
			else // (tmpTransformedKmerSeq != lastTransformedKmer_uniq) 
			{
				transformedKmerCount_ofs << lastTransformedKmer_uniq 
					<< "\t" << lastTransformedKmer_uniq_count << endl;
				lastTransformedKmer_uniq = tmpTransformedKmerSeq;
				lastTransformedKmer_uniq_count = tmpTransformedKmerCount;
			}
		}
	}
	transformedKmerCount_ofs.close();
	rawKmerCountFile_transformed_sorted_ifs.close();

	// remove intermediate dir file
	if(remove_intermediate_dir_file_bool)
	{
		string rm_transformed_file = "rm " + outputRawKmerCountFile_transformed;
		string rm_transformed_sorted_file = "rm " + outputRawKmerCountFile_transformed_sorted;
		string rm_tmp_sort_dir = "rm -r " + tmp_sort_dir;
		//string rm_transformed_sorted_uniq_file = "rm " + outputRawKmerCountFile_transformed_sorted_uniq;
		system(rm_transformed_file.c_str());
		system(rm_transformed_sorted_file.c_str());
		system(rm_tmp_sort_dir.c_str());
		//system(rm_transformed_sorted_uniq_file.c_str());
	}
	return 0;
}