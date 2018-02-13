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

using namespace std;

int return_index_phylumNameVec(string& tmpPhylumStr)
{
	int tmpPhylumStrLength = tmpPhylumStr.length();
	switch(tmpPhylumStrLength)
	{
		case 6:
			if(tmpPhylumStr == "Random")
				return 17;
			else
				return -1;
			break;
		case 8:
			if(tmpPhylumStr == "Chlorobi")
				return 3;
			else
				return -1;
			break;
		case 9:
			if(tmpPhylumStr == "Simulated")
				return 16;
			else
				return -1;
			break;
		case 10:
			if(tmpPhylumStr == "Eukaryotes")
				return 8;
			else if(tmpPhylumStr == "Firmicutes")
				return 10;
			else
				return -1;
			break;
		case 11:
			if(tmpPhylumStr == "Chloroflexi")
				return 4;
			if(tmpPhylumStr == "Nitrospirae")
				return 12;
			else
				return -1;
			break;
		case 13:
			if(tmpPhylumStr == "Acidobacteria")
				return 0;
			else if(tmpPhylumStr == "Bacteroidetes")
				return 2;
			else if(tmpPhylumStr == "Crenarchaeota")
				return 5;
			else if(tmpPhylumStr == "Cyanobacteria")
				return 6;
			else if(tmpPhylumStr == "Elusimicrobia")
				return 7;
			else if(tmpPhylumStr == "Euryarchaeota")
				return 9;	
			else
				return -1;
			break;
		case 14:
			if(tmpPhylumStr == "Actinobacteria")
				return 1;
			else if(tmpPhylumStr == "Planctomycetes")
				return 13;
			else if(tmpPhylumStr == "Proteobacteria")
				return 14;							
			else
				return -1;
			break;
		case 15:
			if(tmpPhylumStr == "Verrucomicrobia")
				return 15;
			else
				return -1;
			break;
		case 16:
			if(tmpPhylumStr == "Gemmatimonadetes")
				return 11;
			else
				return -1;
			break;
		default:
			return -1;
	}
}

void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
{
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("_", startLoc);
		if(tabLoc == string::npos)
			break;
		string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
		tmpFieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpFieldVec.push_back(tmpStr.substr(startLoc));
}	

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputFq_1 inputFq_2 outputFolder" << endl;
		exit(1);
	}
	string inputFq_1 = argv[1];
	string inputFq_2 = argv[2];

	string outputFolder = argv[3];
	outputFolder += "/";
	string mkdir_outputFolder = "mkdir " + outputFolder;
	system(mkdir_outputFolder.c_str());

	vector<string> phylumNameVec;
	phylumNameVec.push_back("Acidobacteria"); // 0
	phylumNameVec.push_back("Actinobacteria"); // 1
	phylumNameVec.push_back("Bacteroidetes"); // 2
	phylumNameVec.push_back("Chlorobi"); // 3
	phylumNameVec.push_back("Chloroflexi"); 
	phylumNameVec.push_back("Crenarchaeota"); //5
	phylumNameVec.push_back("Cyanobacteria"); //6
	phylumNameVec.push_back("Elusimicrobia"); //7
	phylumNameVec.push_back("Eukaryotes"); 
	phylumNameVec.push_back("Euryarchaeota");
	phylumNameVec.push_back("Firmicutes"); //10
	phylumNameVec.push_back("Gemmatimonadetes"); // 11
	phylumNameVec.push_back("Nitrospirae"); // 12
	phylumNameVec.push_back("Planctomycetes"); //13
	phylumNameVec.push_back("Proteobacteria");// 14
	phylumNameVec.push_back("Verrucomicrobia"); // 15
	phylumNameVec.push_back("Simulated");
	phylumNameVec.push_back("Random"); //17

	vector<ofstream*> phylum_ofs_vec_1;
	vector<ofstream*> phylum_ofs_vec_2;
	for(int tmp = 0; tmp < phylumNameVec.size(); tmp++)
	{
		string tmpPhylum = phylumNameVec[tmp];
		string tmpPhylum_fq_1 = outputFolder + phylumNameVec[tmp] + "_1.fq";
		string tmpPhylum_fq_2 = outputFolder + phylumNameVec[tmp] + "_2.fq";
		ofstream* tmp_fq_1_ofs = new ofstream(tmpPhylum_fq_1.c_str());
		ofstream* tmp_fq_2_ofs = new ofstream(tmpPhylum_fq_2.c_str());
		phylum_ofs_vec_1.push_back(tmp_fq_1_ofs);
		phylum_ofs_vec_2.push_back(tmp_fq_2_ofs);
	}
	string other_fq_1 = outputFolder + "other_1.fq";
	string other_fq_2 = outputFolder + "other_2.fq";
	ofstream other_1_ofs(other_fq_1.c_str());
	ofstream other_2_ofs(other_fq_2.c_str());

	ifstream fq_ifs_1(inputFq_1.c_str());
	ifstream fq_ifs_2(inputFq_2.c_str());
	while(!fq_ifs_1.eof())
	{
		string tmpStr_1;
		getline(fq_ifs_1, tmpStr_1);
		if(tmpStr_1 == "")
			break;
		vector<string> tmpFieldVec;
		parseStr2fieldVec(tmpFieldVec, tmpStr_1);
		string tmpPhylumStr = tmpFieldVec[tmpFieldVec.size() - 2];
		int tmpOfsIndex = return_index_phylumNameVec(tmpPhylumStr);
		string tmpStr_2, tmpStr_3, tmpStr_4, tmpStr_5, tmpStr_6, tmpStr_7, tmpStr_8;
		getline(fq_ifs_1, tmpStr_2);
		getline(fq_ifs_1, tmpStr_3);
		getline(fq_ifs_1, tmpStr_4);
		getline(fq_ifs_2, tmpStr_5);
		getline(fq_ifs_2, tmpStr_6);
		getline(fq_ifs_2, tmpStr_7);
		getline(fq_ifs_2, tmpStr_8);
		if(tmpOfsIndex < 0)
		{
			other_1_ofs << tmpStr_1 << endl << tmpStr_2 << endl << tmpStr_3 << endl << tmpStr_4 << endl;
			other_2_ofs << tmpStr_5 << endl << tmpStr_6 << endl << tmpStr_7 << endl << tmpStr_8 << endl;
		}
		else
		{	
			(*phylum_ofs_vec_1[tmpOfsIndex]) << tmpStr_1 << endl << tmpStr_2 << endl << tmpStr_3 << endl << tmpStr_4 << endl;
			(*phylum_ofs_vec_2[tmpOfsIndex]) << tmpStr_5 << endl << tmpStr_6 << endl << tmpStr_7 << endl << tmpStr_8 << endl;
		}
	}
	other_1_ofs.close();
	other_2_ofs.close();
	fq_ifs_1.close();
	fq_ifs_2.close();
	for(int tmp = 0; tmp < phylumNameVec.size(); tmp++)
	{	
		(*phylum_ofs_vec_1[tmp]).close();
		delete phylum_ofs_vec_1[tmp];
		(*phylum_ofs_vec_2[tmp]).close();
		delete phylum_ofs_vec_2[tmp];		
	}
	return 0;
}