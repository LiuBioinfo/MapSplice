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

time_t nowtime;
struct tm *local;

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 15)
	{
		cout << "Executable inputGenomeTester4binFolder outputFolder outputScriptFile inputFa_1 inputFa_2 inputFa_3 inputFa_4 inputFa_5 inputFa_6 inputFa_7 inputFa_8 kmerLength freq_cutoff num_threads" << endl;
		exit(1);
	}
	string inputGenomeTester4binFolder = argv[1];
	inputGenomeTester4binFolder += "/";

	string outputFolderPath = argv[2];
	outputFolderPath += "/";
	
	string kmersetManipulateFile = argv[3];//outputFolderStr + "kmerset_manipulate.txt";
	ofstream kmersetManipulate_ofs(kmersetManipulateFile.c_str());

	vector<string> faFileVec;
	for(int tmp = 0; tmp < 8; tmp++)
		faFileVec.push_back(argv[4+tmp]);

	string kmerLengthStr = argv[12];
	string freq_cutoff_str = argv[13];
	string num_threads_str = argv[14];

	kmersetManipulate_ofs << "cd " << inputGenomeTester4binFolder << endl << endl;
	
	// generate kmer set for each fasta file
	for(int tmp = 0; tmp < 8; tmp ++)
		kmersetManipulate_ofs << "./glistmaker " << faFileVec[tmp] << " --outputname " << outputFolderPath << "/kmer." << tmp + 1
			<< " --wordlength " << kmerLengthStr << " --num_threads " << num_threads_str << " --cutoff " << freq_cutoff_str << endl;
	kmersetManipulate_ofs << endl;
	// merge kmer set file
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.1" << "_" << kmerLengthStr << ".list "
		<< outputFolderPath << "/kmer.2" << "_" << kmerLengthStr << ".list --union -o " << outputFolderPath << "/kmer.1_2" << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.3" << "_" << kmerLengthStr << ".list "
		<< outputFolderPath << "/kmer.4" << "_" << kmerLengthStr << ".list --union -o " << outputFolderPath << "/kmer.3_4" << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.5" << "_" << kmerLengthStr << ".list "
		<< outputFolderPath << "/kmer.6" << "_" << kmerLengthStr << ".list --union -o " << outputFolderPath << "/kmer.5_6" << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.7" << "_" << kmerLengthStr << ".list "
		<< outputFolderPath << "/kmer.8" << "_" << kmerLengthStr << ".list --union -o " << outputFolderPath << "/kmer.7_8" << endl;

	kmersetManipulate_ofs << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.1" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.1" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.1" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.2" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.2" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.2" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.3" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.3" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.3" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.4" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.4" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.4" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.5" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.5" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.5" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.6" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.6" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.6" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.7" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.7" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.7" << "_" << kmerLengthStr << ".list" << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.8" << "_" << kmerLengthStr << ".list > "
		<< outputFolderPath << "/kmer.8" << "_" << kmerLengthStr << ".list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.8" << "_" << kmerLengthStr << ".list" << endl;

	kmersetManipulate_ofs << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.1_2_" << kmerLengthStr << "_union.list "
		<< outputFolderPath << "/kmer.3_4_" << kmerLengthStr << "_union.list --union -o " << outputFolderPath << "/kmer.1_2_3_4" << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.5_6_" << kmerLengthStr << "_union.list "
		<< outputFolderPath << "/kmer.7_8_" << kmerLengthStr << "_union.list --union -o " << outputFolderPath << "/kmer.5_6_7_8" << endl;

	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.1_2_" << kmerLengthStr << "_union.list" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.3_4_" << kmerLengthStr << "_union.list" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.5_6_" << kmerLengthStr << "_union.list" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.7_8_" << kmerLengthStr << "_union.list" << endl;

	kmersetManipulate_ofs << endl;
	kmersetManipulate_ofs << "./glistcompare " << outputFolderPath << "/kmer.1_2_3_4_" << kmerLengthStr << "_union.list "
		<< outputFolderPath << "/kmer.5_6_7_8_" << kmerLengthStr << "_union.list --union -o " << outputFolderPath << "kmer.1_2_3_4_5_6_7_8" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.1_2_3_4_" << kmerLengthStr << "_union.list" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.5_6_7_8_" << kmerLengthStr << "_union.list" << endl;

	kmersetManipulate_ofs << endl;
	kmersetManipulate_ofs << "./glistquery " << outputFolderPath << "/kmer.1_2_3_4_5_6_7_8_" << kmerLengthStr << "_union.list > "
	 	<< outputFolderPath << "/kmer.1_2_3_4_5_6_7_8_" << kmerLengthStr << "_union.list.readable" << endl;
	kmersetManipulate_ofs << "rm " << outputFolderPath << "/kmer.1_2_3_4_5_6_7_8_" << kmerLengthStr << "_union.list" << endl;
	
	kmersetManipulate_ofs.close();
	return 0;
}