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

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputFaOrFqFile outputSeqOnlyFile" << endl;
		exit(1);
	}
	bool fa_or_fq_bool;
	string inputFaOrFqFile = argv[1];
	string outputSeqOnlyFile = argv[2];
	int inputFaOrFqFile_length = inputFaOrFqFile.length();	
	if(inputFaOrFqFile_length > 6)
	{
		string tmpSuffix_len3 = inputFaOrFqFile.substr(inputFaOrFqFile_length-3);
		string tmpSuffix_len6 = inputFaOrFqFile.substr(inputFaOrFqFile_length-6);		
		if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa")
			||(tmpSuffix_len6 == ".fasta")||(tmpSuffix_len6 == ".FASTA")||(tmpSuffix_len6 == ".Fasta"))
			fa_or_fq_bool = true;
		else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq")
				||(tmpSuffix_len6 == ".fastq")||(tmpSuffix_len6 == ".FASTQ")||(tmpSuffix_len6 == ".Fastq"))
			fa_or_fq_bool = false;
		else
		{
			cout << inputFaOrFqFile << " is not a regular fasta or fastq format file name !" << endl;
			exit(1);
		}		
	}
	else if(inputFaOrFqFile_length > 3)
	{
		string tmpSuffix_len3 = inputFaOrFqFile.substr(inputFaOrFqFile_length-3);
		if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa"))
			fa_or_fq_bool = true;
		else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq"))
			fa_or_fq_bool = false;
		else
		{
			cout << inputFaOrFqFile << " is not a regular fasta or fastq format file name !" << endl;
			exit(1);
		}
	}
	else
	{
		cout << inputFaOrFqFile << " is not a regular fasta or fastq format file name !" << endl;
		exit(1);		
	}

	ifstream fa_fq_ifs(inputFaOrFqFile.c_str());	
	ofstream seq_ofs(outputSeqOnlyFile.c_str());
	while(!fa_fq_ifs.eof())
	{	
		string tmpName;
		getline(fa_fq_ifs, tmpName);
		if(tmpName == "")
			break;
		string tmpSeq;
		getline(fa_fq_ifs, tmpSeq);
		seq_ofs << tmpSeq << endl;
		if(!fa_or_fq_bool)
		{
			string tmpComment, tmpQual;
			getline(fa_fq_ifs, tmpComment);
			getline(fa_fq_ifs, tmpQual);
		}
	}
	seq_ofs.close();
	fa_fq_ifs.close();
	return 0;
}