// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef BUILD_SNPMERINDEX_H
#define BUILD_SNPMERINDEX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sys/types.h>    
#include <dirent.h>    
#include <stdio.h>    
#include <errno.h>
#include <set>
#include "../buildIndex/buildGlobalAndLocalIndex/buildIndexParameter.h"

using namespace std;

class Build_snpMerIndex
{
private:
	int Range;
	//int LONGLCP;
	unsigned int NULL_NUM;
	//#define PREINDEX_STRINGLENGTH 14
	//#define CANDALILOC 100
	//#define SEGMENTNUM 20
	//#define minValSegLength 20
	//#define PreIndexSize 268435456

public:
	Build_snpMerIndex()
	{
		Range = 6;
		//LONGLCP = 255;
		NULL_NUM = 4294967290;
	}

	bool build_snpMerIndex(string& inputMergedFastaFile, string& outputIndexFolder)
	{
		BuildIndexParameter_info* buildIndexInfo = new BuildIndexParameter_info();
		string OutputIndexFolder = outputIndexFolder;
		buildIndexInfo->OutputIndexFolderStr = OutputIndexFolder + "/";		
		//cout << "mkdir " << endl;
		string mkdir_cmd_outputIndex = "mkdir -p " + OutputIndexFolder;
		system(mkdir_cmd_outputIndex.c_str());
		string parameter_file = buildIndexInfo->OutputIndexFolderStr + "/_parameter";
		ofstream parameter_file_ofs(parameter_file.c_str());
		string chrom_file = buildIndexInfo->OutputIndexFolderStr + "/_chromMerge";
		ofstream chrom_file_ofs(chrom_file.c_str(),ios::binary);
		string log_file = buildIndexInfo->OutputIndexFolderStr + "/_progress.log";
		ofstream log_ofs(log_file.c_str());

		log_ofs << "start to read merged fasta file" << endl;
		string mergedFaFile = inputMergedFastaFile;
		ifstream mergedFa_ifs(mergedFaFile.c_str());
		while(!mergedFa_ifs.eof())
		{
			string tmpChrNameStr, tmpChrSeqStr;
			getline(mergedFa_ifs, tmpChrNameStr);
			if(tmpChrNameStr == "")
				break;
			tmpChrNameStr += ".fa";
			getline(mergedFa_ifs, tmpChrSeqStr);
			(buildIndexInfo->chrNameStrVec).push_back(tmpChrNameStr.substr(1));
			int tmpChrSeqLength = tmpChrSeqStr.length();
			(buildIndexInfo->chromLengthVec).push_back(tmpChrSeqLength);
			(buildIndexInfo->chromSeqVec).push_back(tmpChrSeqStr);
		}
		mergedFa_ifs.close();
		log_ofs << endl << " ****************   getChrEndPosVec() starts ... ***********************" << endl << endl;
		buildIndexInfo->getChrEndPosVec(log_ofs);
		log_ofs << endl << " ****************   outputChromSeq() starts ...   **************************" << endl << endl; 
		buildIndexInfo->outputChromSeq_mergedFa(chrom_file_ofs, log_ofs);
		log_ofs << endl << " ****************   outputIndexParameter() starts ... **********************" << endl << endl;	
		buildIndexInfo->outputIndexParameter(parameter_file_ofs);
	
		/////////////////////////////////////////////////
		//////generate original size index
		/////////////////////////////////////////////////
		log_ofs << endl << " *****************   start to generate original size index   **********************" << endl << endl;
		unsigned int MAX = (buildIndexInfo->chrEndPosInGenomeVec)[(buildIndexInfo->chrNameStrVec).size()-1] + 1 + 1;
		log_ofs << "MAX: " << MAX << endl;
		string chrom_file_str = chrom_file;
		string SA_file = buildIndexInfo->OutputIndexFolderStr; SA_file.append("/_SA"); ofstream SA_file_ofs(SA_file.c_str(),ios::binary); 
		string lcp_file = buildIndexInfo->OutputIndexFolderStr; lcp_file.append("/_lcp"); ofstream lcp_file_ofs(lcp_file.c_str(),ios::binary);
		string up_file = buildIndexInfo->OutputIndexFolderStr; up_file.append("/_up"); ofstream up_file_ofs(up_file.c_str(),ios::binary);
		string down_file = buildIndexInfo->OutputIndexFolderStr; down_file.append("/_down"); ofstream down_file_ofs(down_file.c_str(),ios::binary);
		string next_file = buildIndexInfo->OutputIndexFolderStr; next_file.append("/_next"); ofstream next_file_ofs(next_file.c_str(),ios::binary);
		string chrom_bit_file = buildIndexInfo->OutputIndexFolderStr; chrom_bit_file.append("/_chrom"); ofstream chrom_bit_file_ofs(chrom_bit_file.c_str(),ios::binary);
	  	FILE *fp = fopen(chrom_file_str.c_str(), "r");
	    unsigned int *r = (unsigned int*)malloc(MAX * sizeof(unsigned int));
		char ch;
		char head[100];
		char base[6] = {'X','A','C','G','T','N'};//char base[Range] = {'X','A','C','G','T','N'};
		char *chrom = (char*)malloc(MAX * sizeof(char));	

		unsigned int chrom_base_num = 0;
		while((ch = fgetc(fp)) != EOF)
		{   
			//printf("ch = %c\n",ch); 
			if((ch == 'A')||(ch == 'a')) {chrom[chrom_base_num] = 'A'; r[chrom_base_num] = 1; chrom_base_num++;}
			else if((ch == 'C')||(ch == 'c')) {chrom[chrom_base_num] = 'C'; r[chrom_base_num] = 2; chrom_base_num++;}
			else if((ch == 'G')||(ch == 'g')) {chrom[chrom_base_num] = 'G'; r[chrom_base_num] = 3; chrom_base_num++;}
			else if((ch == 'T')||(ch == 't')) {chrom[chrom_base_num] = 'T'; r[chrom_base_num] = 4; chrom_base_num++;}
			else if((ch == 'N')||(ch == 'n')) 
			{
				chrom[chrom_base_num] = 'N'; r[chrom_base_num] = 5; chrom_base_num++;
			}
			else if((ch == 'X')) {chrom[chrom_base_num] = 'X'; r[chrom_base_num] = 6; chrom_base_num++;}
			else if((ch == '\t')||(ch == '\n')) {continue;}
			else 
			{
				printf("\n found illegal input is '%c', For now, replace it with N",ch); 
				log_ofs << endl << "illegal input in provided reference Fasta file: " << ch 
					<< endl << "For now, replace it with N" << endl;
				//break;
				chrom[chrom_base_num] = 'N'; r[chrom_base_num] = 5; chrom_base_num++;
			}
		}
		fclose(fp);
		log_ofs << "the number of bases in Chromo is "<< chrom_base_num << endl;
		log_ofs << "chrom is ready" << endl;
		r[MAX-1] = Range;
		chrom[MAX-1] = 'X';	
		log_ofs << "chrom[MAX-3]: " << chrom[MAX-3] << endl;	
		log_ofs << "chrom[MAX-2]: " << chrom[MAX-2] << endl;
		log_ofs << "chrom[MAX-1]: " << chrom[MAX-1] << endl;
		chrom_bit_file_ofs.write((const char*) chrom, MAX * sizeof(char));
	    unsigned int *sa = (unsigned int*)malloc(MAX * sizeof(unsigned int));
		log_ofs << "start to build SA array" << endl;
		//Xinan's SA generation algorithm
		//da(r,sa,MAX,Range+1); 
		//Kyle's SA generation algorithm
		//cout << "MAX: " << MAX << endl;
		//cout << "Range + 1: " << Range + 1 << endl;
		bool tmpSuffixArrayBool = this->suffixArray_bool(r, sa, MAX, Range+1);
		if(!tmpSuffixArrayBool)
			return false;
		log_ofs << "SA is ready" << endl;	
	    //////////////////////////////////////////////////
	    //unsigned int rank[MAX]={0}, lcp[MAX]={0}, up[MAX] = {0}, down[MAX] = {0}, next[MAX] = {0};  
	    unsigned int *rank = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	    unsigned int *lcp = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	    unsigned int *up = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	    unsigned int *down = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	    unsigned int *next = (unsigned int*)malloc(MAX * sizeof(unsigned int));

	    log_ofs << "start to build LCP array"<< endl;
		this->build_lcp(r, sa, lcp, rank, MAX); //build lcp array
		free(r); free(rank);
		log_ofs << "lcp is ready" << endl;	
		log_ofs << "start to build up & down array" << endl;
		this->build_up_down(lcp, up, down, MAX); // build up and down array
		log_ofs << "up & down are ready" << endl;
		log_ofs << "start to build next array" << endl;		
		this->build_next(lcp, next, MAX); //build next array
		log_ofs << "next is ready" << endl;	

		log_ofs << "start to output Index to file "<< endl;
		SA_file_ofs.write((const char*) sa, MAX * sizeof(unsigned int));
		lcp_file_ofs.write((const char*) lcp, MAX * sizeof(unsigned int));
		up_file_ofs.write((const char*) up, MAX * sizeof(unsigned int));
		down_file_ofs.write((const char*) down, MAX * sizeof(unsigned int));
		next_file_ofs.write((const char*) next, MAX * sizeof(unsigned int));
		//free(sa);
		log_ofs << "finish writing index to files " << endl;
		/////////////////////////////////////////////////////////////////////
		/////  compress index for whole genome original size index
		/////////////////////////////////////////////////////////////////////
		unsigned int indexSize = MAX;
		string childTab_file = buildIndexInfo->OutputIndexFolderStr; 
		childTab_file.append("_childTab"); ofstream childTab_file_ofs(childTab_file.c_str(),ios::binary);
		string detChild_file = buildIndexInfo->OutputIndexFolderStr; 
		detChild_file.append("_detChild"); ofstream detChild_file_ofs(detChild_file.c_str(), ios::binary);

	    unsigned int *childTab;
	    childTab = (unsigned int*)malloc(indexSize * sizeof(unsigned int));
	   	BYTE *verifyChild;
		verifyChild = (BYTE*)malloc(indexSize * sizeof(BYTE));
		log_ofs << "start to compress original size index " << endl;
		BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();
		bool tmpUpDownNextCompress_success_bool = 
			tmpBuildIndexInfo->compressUpDownNext2ChildtabVerifyChild_returnSuccessOrNotBool(up, down, next, childTab, verifyChild, indexSize);
		if(!tmpUpDownNextCompress_success_bool)
		{
			delete tmpBuildIndexInfo;
			free(lcp); free(up); free(down); free(next); free(childTab); free(verifyChild);
			lcp = NULL;
			up = NULL;
			down = NULL;
			next = NULL;
			childTab = NULL;
			verifyChild = NULL;
			return false;
		}
		//log_ofs << "finish compressing index" << endl;
		childTab_file_ofs.write((const char*) childTab, indexSize * sizeof(unsigned int));
		detChild_file_ofs.write((const char*) verifyChild, indexSize * sizeof(BYTE));	
		childTab_file_ofs.close();
		detChild_file_ofs.close();
		log_ofs << "finish compressing index" << endl; 
		free(up); free(down); free(next); //free(childTab); free(verifyChild);
		string lcpCompress_file = buildIndexInfo->OutputIndexFolderStr; 
		lcpCompress_file.append("_lcpCompress"); ofstream lcpCompress_file_ofs(lcpCompress_file.c_str(),ios::binary);
		//BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();
		log_ofs << "start to compress Lcp array" << endl;
		BYTE *lcpCompress = (BYTE*)malloc(indexSize * sizeof(BYTE));
		tmpBuildIndexInfo->compressLcp2Lcpcompress(lcp, lcpCompress, indexSize);
		log_ofs << "finish compressing Lcp array " << endl;
		lcpCompress_file_ofs.write((const char*) lcpCompress, indexSize * sizeof(BYTE));	
		free(lcp);
		//free(lcpCompress);
		lcpCompress_file_ofs.close();
		free(childTab);free(chrom);free(verifyChild);	

		string cmd_delete_lcp = "rm " + lcp_file;
		string cmd_delete_up = "rm " + up_file;
		string cmd_delete_down = "rm " + down_file;
		string cmd_delete_next = "rm " + next_file;
		string cmd_delete_chrom = "rm " + chrom_file;
		system(cmd_delete_lcp.c_str());
		system(cmd_delete_up.c_str());
		system(cmd_delete_down.c_str());
		system(cmd_delete_next.c_str());
		system(cmd_delete_chrom.c_str());
	 	log_ofs << endl << "all index building jobs done" << endl;

	 	SA_file_ofs.close();
	 	lcp_file_ofs.close();
		up_file_ofs.close();
		down_file_ofs.close();
		next_file_ofs.close();
		chrom_bit_file_ofs.close();
		delete tmpBuildIndexInfo;
		tmpBuildIndexInfo = NULL;

		delete buildIndexInfo;
		buildIndexInfo = NULL;
		return true;
	}

	// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
	bool radixPass_bool(unsigned int* a, unsigned int* b, unsigned int* r, unsigned int n, int K)
	{
		//cout << "radixPass starts ...." << endl;
		// count occurrences
		unsigned int* c = new unsigned int[K + 1]; // counter array
		//cout << "process 1 starts " << endl;
		for (int i = 0; i <= K; i++)
			c[i] = 0; // reset counters
		//cout << "process 2 starts " << endl;
		//cout << "K:" << K << endl;
		//cout << "n: " << n << endl;
		for (unsigned int i = 0; i < n; i++)
		{
			int tmpAi = a[i];
			if(tmpAi < 0)
				return false;
			int tmpRAi = r[a[i]];
			if(tmpRAi < 0)
				return false;
			c[r[a[i]]]++; // count occurrences
		}
		//cout << "process 3 starts " << endl;
		for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
		{
			int t = c[i];
			c[i] = sum;
			sum += t;
		}
		//cout << "process 4 starts " << endl;
		for (unsigned int i = 0; i < n; i++)
			b[c[r[a[i]]]++] = a[i]; // sort
		//cout << "radix pass ends " << endl;
		delete[] c;
		return true;
	}

	// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
	// static void radixPass(unsigned int* a, unsigned int* b, unsigned int* r, unsigned int n, int K)
	// {
	// 	cout << "radixPass starts ...." << endl;
	// 	// count occurrences
	// 	unsigned int* c = new unsigned int[K + 1]; // counter array
	// 	cout << "process 1 starts " << endl;
	// 	for (int i = 0; i <= K; i++)
	// 		c[i] = 0; // reset counters
	// 	cout << "process 2 starts " << endl;
	// 	cout << "K:" << K << endl;
	// 	cout << "n: " << n << endl;
	// 	for (unsigned int i = 0; i < n; i++)
	// 	{
	// 		int tmpAi = a[i];
	// 		int tmpRAi = r[a[i]];
	// 		cout << "tmpAi: " << tmpAi << endl;
	// 		cout << "tmpRAi: " << tmpRAi << endl;
	// 		c[r[a[i]]]++; // count occurrences
	// 	}
	// 	cout << "process 3 starts " << endl;
	// 	for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
	// 	{
	// 		int t = c[i];
	// 		c[i] = sum;
	// 		sum += t;
	// 	}
	// 	//cout << "process 4 starts " << endl;
	// 	for (unsigned int i = 0; i < n; i++)
	// 		b[c[r[a[i]]]++] = a[i]; // sort
	// 	//cout << "radix pass ends " << endl;
	// 	delete[] c;
	// }

	bool suffixArray_bool(unsigned int* T, unsigned int* SA, unsigned int n, int K) 
	{
		//cout << "suffixArray starts ...." << endl;
		unsigned int n0 = (n + 2) / 3,
			n1 = (n + 1) / 3,
			n2 = n / 3,
			n02 = n0 + n2;
		//cout << "n0: " << n0 << endl;
		//cout << "n1: " << n1 << endl;
		//cout << "n2: " << n2 << endl;
		//cout << "n02: " << n02 << endl;
		unsigned int* R = new unsigned int[n02 + 3];
		R[n02] = R[n02 + 1] = R[n02 + 2] = 0;
		unsigned int* SA12 = new unsigned int[n02 + 3];
		SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
		unsigned int* R0 = new unsigned int[n0];
		unsigned int* SA0 = new unsigned int[n0];
		//cout << "step0 starts ......" << endl;
		//******* Step 0: Construct sample ********
		// generate positions of mod 1 and mod 2 suffixes
		// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
		for (unsigned int i = 0, j = 0; i < n + (n0 - n1); i++)
			if (i % 3 != 0)
				R[j++] = i;
		//cout << "step 0 ends " << endl;
		//******* Step 1: Sort sample suffixes ********
		// lsb radix sort the mod 1 and mod 2 triples
		bool tmpRadixPass_bool_1 = this->radixPass_bool(R, SA12, T + 2, n02, K);
		//cout << "1st radixPass ends " << endl;
		if(!tmpRadixPass_bool_1)
			return false;
		bool tmpRadixPass_bool_2 = this->radixPass_bool(SA12, R, T + 1, n02, K);
		if(!tmpRadixPass_bool_2)
			return false;
		//cout << "2nd radixPass ends " << endl;
		bool tmpRadixPass_bool_3 = this->radixPass_bool(R, SA12, T, n02, K);
		if(!tmpRadixPass_bool_3)
			return false;
		//cout << "3rd radixPass ends " << endl;
		
		// find lexicographic names of triples and
		// write them to correct places in R
		unsigned int name = 0, c0 = 0-1, c1 = 0-1, c2 = 0-1;
		for (unsigned int i = 0; i < n02; i++)
		{
			if (T[SA12[i]] != c0 || T[SA12[i] + 1] != c1 || T[SA12[i] + 2] != c2)
			{
				name++;
				c0 = T[SA12[i]];
				c1 = T[SA12[i] + 1];
				c2 = T[SA12[i] + 2];
			}

			if (SA12[i] % 3 == 1)
			{
				// write to R1
				R[SA12[i] / 3] = name;
			}
			else
			{
				// write to R2
				R[SA12[i] / 3 + n0] = name;
			}
		}
		//cout << "find lexicographic names of triples and write them to correct places in R ends " << endl;
		// recurse if names are not yet unique
		if (name < n02)
		{
			bool tmpSuffixArray_bool = this->suffixArray_bool(R, SA12, n02, name);
			if(!tmpSuffixArray_bool)
				return false;
			// store unique names in R using the suffix array
			for (unsigned int i = 0; i < n02; i++)
				R[SA12[i]] = i + 1;
		}
		else // generate the suffix array of R directly
			for (unsigned int i = 0; i < n02; i++)
				SA12[R[i] - 1] = i;
		//cout << "recurse if names are not yet unique ends" << endl;
		//******* Step 2: Sort nonsample suffixes ********
		// stably sort the mod 0 suffixes from SA12 by their first character
		for (unsigned int i = 0, j = 0; i < n02; i++)
			if (SA12[i] < n0)
				R0[j++] = 3 * SA12[i];
		bool tmpRadixPass_bool_4 = radixPass_bool(R0, SA0, T, n0, K);
		if(!tmpRadixPass_bool_4)
			return false;
		//cout << "step2 ends ...." << endl;
		//******* Step 3: Merge ********
		// merge sorted SA0 suffixes and sorted SA12 suffixes
		for (unsigned int p = 0, t = n0 - n1, k = 0; k < n; k++)
		{
			#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
			unsigned int i = GetI(); // pos of current offset 12 suffix
			unsigned int j = SA0[p]; // pos of current offset 0 suffix
			if (SA12[t] < n0 // different compares for mod 1 and mod 2 suffixes
					? leq(T[i], R[SA12[t] + n0], T[j], R[j / 3])
					: leq(T[i], T[i + 1], R[SA12[t] - n0 + 1], T[j], T[j + 1], R[j / 3 + n0]))
			{ 	// suffix from SA12 is smaller

				SA[k] = i;
				t++;
				if (t == n02) // done --- only SA0 suffixes left
					for (k++; p < n0; p++, k++)
						SA[k] = SA0[p];
			}
			else // suffix from SA0 is smaller
			{
				SA[k] = j;
				p++;
				if (p == n0) // done --- only SA12 suffixes left
					for (k++; t < n02; t++, k++)
						SA[k] = GetI();
			}
		}
		//cout << "step3 ends ...." << endl;
		delete[] R;
		delete[] SA12;
		delete[] SA0;
		delete[] R0;
		return true;
	}

	inline bool leq(int a1, int a2, int b1, int b2) // lexicographic order
	{
		return (a1 < b1 || (a1 == b1 && a2 <= b2));
	} // for pairs

	inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3)
	{
		return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
	} // and triples


	void build_up_down(unsigned int *lcptab, unsigned int *up, unsigned int *down, unsigned int n)
	{
	    unsigned int lastIndex = NULL_NUM;	
		stack<unsigned int> up_down;
		up_down.push(0);
		unsigned int i;
		for (i = 0; i < n; i++)
		{
			while (lcptab[i] < lcptab[up_down.top()])
			{	lastIndex = up_down.top();
				up_down.pop();
				if((lcptab[i] <= lcptab[up_down.top()]) && (lcptab[up_down.top()] != lcptab[lastIndex]))
					down[up_down.top()] = lastIndex;
			}
			// now lcptab[i] >= lcptab[up_down.top()] holds
			if(lastIndex != NULL_NUM)
			{
				up[i] = lastIndex;
				lastIndex = NULL_NUM;
			}
			up_down.push(i);	
		}
		return;
	}

	void build_next(unsigned int *lcptab, unsigned int *next, unsigned int n)
	{
		stack<unsigned int> nextIndex;
		unsigned int lastIndex;
		unsigned int j;
		nextIndex.push(0);
		for(j = 0; j < n; j++)
		{
			while(lcptab[j] < lcptab[nextIndex.top()])
				nextIndex.pop();
			if(lcptab[j] == lcptab[nextIndex.top()])
			{
				lastIndex = nextIndex.top();
				nextIndex.pop();
				next[lastIndex] = j;
			}
			nextIndex.push(j);
		}	
		return;
	}

	void build_lcp(unsigned int *r, unsigned int *sa, unsigned int *lcp, unsigned int *rank, unsigned int n)
	{

		unsigned int i, j; 
		int k=0;
		for (i = 0; i < n; i++) rank[sa[i]] = i;
			//cout << "build lcp - step 1" << endl;
		for (i = 0; i < n; lcp[rank[i++]] = k) 
		for (k?k--:0, (rank[i] == 0)?(j=0):(j=sa[rank[i]-1]); r[i+k] == r[j+k]; k++);
		lcp[0] = 0;
			//cout << "build lcp - step 2" << endl; 
		return;
	}
};
#endif