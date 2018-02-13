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

void getFileIncludeOrExcludeCodeVec(int val, int binBoolVecSize, vector<bool>& binBoolVec)  
{
    for(int i = binBoolVecSize - 1; i >= 0; i--)  
    {  
        if(val & (1 << i))  
            binBoolVec.push_back(1);  
        else  
            binBoolVec.push_back(0);  
    }  
}

int main(int argc, char** argv)
{
	if(argc < 8)
	{
		cout << "Executable inputGenomeTester4binFolder inputSetNum outputFolder outputScriptFile inputKmerSet_1 inputKmerSet_2 ..." << endl;
		exit(1);
	}
	string inputGenomeTester4binFolder = argv[1];
	inputGenomeTester4binFolder += "/";

	string inputSetNumStr = argv[2];
	int setNum = atoi(inputSetNumStr.c_str());
	int kmerSetFileNum = argc - 5;
	if(setNum != kmerSetFileNum)
	{
		cout << "setNum != kmerSetFileNum" << endl;
		cout << "setNum: " << setNum << endl;
		cout << "kmerSetFileNum: " << kmerSetFileNum << endl;
		exit(1);
	}

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());
	log_ofs << "Command Line:" << endl;
	for(int tmp = 0; tmp < argc; tmp++)
		log_ofs << " " << argv[tmp];
	log_ofs << endl;

	string outputCmdFile = outputFolderStr + "cmd.txt";
	ofstream cmd_ofs(outputCmdFile.c_str());

	string kmersetManipulateFile = argv[4];//outputFolderStr + "kmerset_manipulate.txt";
	ofstream kmersetManipulate_ofs(kmersetManipulateFile.c_str());

	log_ofs << endl << "start to load Kmer set files " << endl;
	vector<string> kmerFilePathVec;
	for(int tmp = 0; tmp < setNum; tmp++)
	{
		string tmpKmerSetFile = argv[5 + tmp];
		kmerFilePathVec.push_back(tmpKmerSetFile);
	}

	log_ofs << endl << "start to generate K-mer set manipulating function code " << endl;
 	int totalGeneratedClassNum = pow(2, setNum) - 1;
  	cout << "totalGeneratedClassNum: " << totalGeneratedClassNum << endl;
  	log_ofs << "totalGeneratedClassNum: " << totalGeneratedClassNum << endl;
  	vector< vector<bool> > fileIncludeOrExcludeCodeVecVec;
	for(int tmpManipulateDecCode = 1; tmpManipulateDecCode <= totalGeneratedClassNum; tmpManipulateDecCode++)
	{
	    //cout << "tmpManipulateDecCode: " << tmpManipulateDecCode << endl;
	    log_ofs << "tmpManipulateDecCode: " << tmpManipulateDecCode << endl;
        vector<bool> tmpFileIncludeOrExcludeCodeVec;
        getFileIncludeOrExcludeCodeVec(tmpManipulateDecCode, setNum, tmpFileIncludeOrExcludeCodeVec);
        fileIncludeOrExcludeCodeVecVec.push_back(tmpFileIncludeOrExcludeCodeVec);
        //for(int tmpKmersetIndex = 0; tmpKmersetIndex < kmersetNum; tmpKmersetIndex ++)
        //    cout << tmpFileIncludeOrExcludeCodeVec[tmpKmersetIndex];
        //cout << endl;
        for(int tmpKmersetIndex = 0; tmpKmersetIndex < setNum; tmpKmersetIndex ++)
            log_ofs << tmpFileIncludeOrExcludeCodeVec[tmpKmersetIndex];
        log_ofs << endl;        
	}

	log_ofs << endl << "start to get all union & intersection of kmer sets for all combinations " << endl;
	kmersetManipulate_ofs << "cd " << inputGenomeTester4binFolder << endl << endl;
	for(int tmpManipulateDecCode = 1; tmpManipulateDecCode <= totalGeneratedClassNum; tmpManipulateDecCode++)
	{		
		int tmpFileIncludeOrExcludeCodeVecVec_index = tmpManipulateDecCode - 1;
		kmersetManipulate_ofs << "#Class: " << tmpManipulateDecCode << endl << "#";
		//kmersetManipulate_ofs << "\t";
		for(int tmpKmersetIndex = 0; tmpKmersetIndex < setNum; tmpKmersetIndex ++)
			kmersetManipulate_ofs << (fileIncludeOrExcludeCodeVecVec[tmpFileIncludeOrExcludeCodeVecVec_index])[tmpKmersetIndex];
		kmersetManipulate_ofs << endl;
		vector<int> tmp_intersection_member_kmerSetIndex_vec_included;
		vector<int> tmp_union_member_kmerSetIndex_vec_excluded;
		for(int tmpKmersetIndex = 0; tmpKmersetIndex < setNum; tmpKmersetIndex ++)
		{
			if((fileIncludeOrExcludeCodeVecVec[tmpFileIncludeOrExcludeCodeVecVec_index])[tmpKmersetIndex])
				tmp_intersection_member_kmerSetIndex_vec_included.push_back(tmpKmersetIndex);
			else
				tmp_union_member_kmerSetIndex_vec_excluded.push_back(tmpKmersetIndex);
		}
		kmersetManipulate_ofs << "#To included intersection kmer set members: " << endl << "#";
		int tmp_intersection_member_kmerSetIndex_vec_included_num 
			= tmp_intersection_member_kmerSetIndex_vec_included.size();
		for(int tmp = 0; tmp < tmp_intersection_member_kmerSetIndex_vec_included_num; tmp++)
			kmersetManipulate_ofs << tmp_intersection_member_kmerSetIndex_vec_included[tmp] << ",";
		kmersetManipulate_ofs << endl << "#To excluded union kmer set members: " << endl << "#";
		int tmp_union_member_kmerSetIndex_vec_excluded_num
			= tmp_union_member_kmerSetIndex_vec_excluded.size();
		for(int tmp = 0; tmp < tmp_union_member_kmerSetIndex_vec_excluded_num; tmp++)
			kmersetManipulate_ofs << tmp_union_member_kmerSetIndex_vec_excluded[tmp] << ",";
		kmersetManipulate_ofs << endl;

		// get tmp_final_intersection_included_kmer
		kmersetManipulate_ofs << "#Command line to get tmp_final_intersection_included_kmer:" << endl;
		if(tmp_intersection_member_kmerSetIndex_vec_included_num <= 0)
		{
			cout << "error ! tmp_intersection_member_kmerSetIndex_vec_included !" << endl;
			exit(1);
		}
		else if(tmp_intersection_member_kmerSetIndex_vec_included_num == 1)
			kmersetManipulate_ofs << "cp " << kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[0]]
				<< " " << outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << "_20_intrsec.list" << endl;
		else if(tmp_intersection_member_kmerSetIndex_vec_included_num == 2)
			kmersetManipulate_ofs << "./glistcompare  " 
				<< kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[0]] 
				<< " " << kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[1]] 
				<< " --intersection --outputname " << outputFolderStr 
				<< "final.intersection_included." << tmpManipulateDecCode << endl;
		else
		{
			int tmpTempFileNum = tmp_intersection_member_kmerSetIndex_vec_included_num - 2;
			//intersection
			// get 1st temporary intersection kmer set file
			kmersetManipulate_ofs << "./glistcompare  "
				<< kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[0]] 
				<< " " << kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[1]]
				<< " --intersection --outputname " << outputFolderStr
				<< "final.intersection_included." << tmpManipulateDecCode << ".tmp.0" << endl;
			// get other temporary intersection kmer set file
			if(tmpTempFileNum > 1)
			{
				for(int tmpTempFileIndex = 1; tmpTempFileIndex < tmpTempFileNum; tmpTempFileIndex ++)
					kmersetManipulate_ofs << "./glistcompare "
						<< outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << ".tmp." << tmpTempFileIndex - 1
						<< "_20_intrsec.list " << kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[1+tmpTempFileIndex]]
						<< " --intersection --outputname " << outputFolderStr << "final.intersection_included." 
						<< tmpManipulateDecCode << ".tmp." << tmpTempFileIndex << endl;
			}
			// get final intersection kmer set file from the last temp kmer set intersection file
			kmersetManipulate_ofs << "./glistcompare " << outputFolderStr 
				<< "final.intersection_included." << tmpManipulateDecCode << ".tmp." << tmpTempFileNum - 1 << "_20_intrsec.list "
				<< kmerFilePathVec[tmp_intersection_member_kmerSetIndex_vec_included[tmp_intersection_member_kmerSetIndex_vec_included.size()-1]]
				<< " --intersection --outputname " << outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << endl;
			//remove temp files			
			for(int tmpTempFileIndex = 0; tmpTempFileIndex < tmpTempFileNum; tmpTempFileIndex ++)
				kmersetManipulate_ofs << "rm " << outputFolderStr << "final.intersection_included." 
					<< tmpManipulateDecCode << ".tmp." << tmpTempFileIndex << "_20_intrsec.list" << endl;
		}

		// get tmp_final_union_excluded_kmer
		kmersetManipulate_ofs << "#Command line to get tmp_final_union_excluded_kmer:" << endl;
		if(tmp_union_member_kmerSetIndex_vec_excluded_num < 0)
		{
			cout << "error ! tmp_union_member_kmerSetIndex_vec_excluded_num !" << endl;
			exit(1);
		}
		else if(tmp_union_member_kmerSetIndex_vec_excluded_num == 0)
		{}
		else if(tmp_union_member_kmerSetIndex_vec_excluded_num == 1)
			kmersetManipulate_ofs << "cp " << kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[0]]
				<< " " << outputFolderStr << "final.union_excluded." << tmpManipulateDecCode << "_20_union.list" << endl;
		else if(tmp_union_member_kmerSetIndex_vec_excluded_num == 2)
			kmersetManipulate_ofs << "./glistcompare "
				<< kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[0]] 
				<< " " << kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[1]] 
				<< " --union --outputname " << outputFolderStr
				<< "final.union_excluded." << tmpManipulateDecCode << endl;
		else
		{
			int tmpTempFileNum = tmp_union_member_kmerSetIndex_vec_excluded_num - 2;
			//union
			// get 1st temporary union kmer set file
			kmersetManipulate_ofs << "./glistcompare "
				<< kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[0]] 
				<< " " << kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[1]]
				<< " --union --outputname " << outputFolderStr
				<< "final.union_excluded." << tmpManipulateDecCode << ".tmp.0" << endl;
			// get other temporary union kmer set file
			if(tmpTempFileNum > 1)
			{
				for(int tmpTempFileIndex = 1; tmpTempFileIndex < tmpTempFileNum; tmpTempFileIndex ++)
					kmersetManipulate_ofs << "./glistcompare "
						<< outputFolderStr << "final.union_excluded." << tmpManipulateDecCode << ".tmp." << tmpTempFileIndex - 1
						<< "_20_union.list " << kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[1+tmpTempFileIndex]]
						<< " --union --outputname " << outputFolderStr << "final.union_excluded."
						<< tmpManipulateDecCode << ".tmp." << tmpTempFileIndex << endl;
			}
			// get final union kmer set file from the last temp kmer set union file
			kmersetManipulate_ofs << "./glistcompare " << outputFolderStr
				<< "final.union_excluded." << tmpManipulateDecCode << ".tmp." << tmpTempFileNum - 1 << "_20_union.list "
				<< kmerFilePathVec[tmp_union_member_kmerSetIndex_vec_excluded[tmp_union_member_kmerSetIndex_vec_excluded.size()-1]]
				<< " --union --outputname " << outputFolderStr << "final.union_excluded." << tmpManipulateDecCode << endl;
			//remove temp files
			for(int tmpTempFileIndex = 0; tmpTempFileIndex < tmpTempFileNum; tmpTempFileIndex ++)
				kmersetManipulate_ofs << "rm " << outputFolderStr << "final.union_excluded." 
					<< tmpManipulateDecCode << ".tmp." << tmpTempFileIndex << "_20_union.list" << endl;
		}

		// substract tmp_final_union_excluded_kmer from tmp_final_intersection_included_kmer
		kmersetManipulate_ofs << "#Command line to substract tmp_final_union_excluded_kmer from tmp_final_intersection_included_kmer:" << endl;
		
		if(tmp_union_member_kmerSetIndex_vec_excluded_num == 0)
		{
			kmersetManipulate_ofs << "./cp " 
				<< outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << "_20_intrsec.list "
				<< outputFolderStr << "final.class." << tmpManipulateDecCode << endl;
			kmersetManipulate_ofs << "rm " << outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << "_20_intrsec.list" << endl;
		}
		else
		{
			kmersetManipulate_ofs << "./glistcompare "
				<< outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << "_20_intrsec.list "
				<< outputFolderStr << "final.union_excluded." << tmpManipulateDecCode << "_20_union.list --difference --outputname "
				<< outputFolderStr << "final.class." << tmpManipulateDecCode << endl;
			kmersetManipulate_ofs << "rm " << outputFolderStr << "final.intersection_included." << tmpManipulateDecCode << "_20_intrsec.list" << endl;
			kmersetManipulate_ofs << "rm " << outputFolderStr << "final.union_excluded." << tmpManipulateDecCode << "_20_union.list" << endl;
		}
		kmersetManipulate_ofs << endl;
	}

	kmersetManipulate_ofs << endl << "#Command line to get all union & intersection of kmer sets for all combinations in readable format" << endl;
	for(int tmpManipulateDecCode = 1; tmpManipulateDecCode <= totalGeneratedClassNum; tmpManipulateDecCode++)
	{
		kmersetManipulate_ofs << "./glistquery " << outputFolderStr << "final.class." << tmpManipulateDecCode
			<< "_20_0_diff1.list > " << outputFolderStr << "final.class." << tmpManipulateDecCode << ".readable" << endl;
		kmersetManipulate_ofs << "rm " << outputFolderStr << "final.class." << tmpManipulateDecCode << "_20_0_diff1.list" << endl;
	}
	// string mergedReadableKmerFile = outputFolderStr + "final.class.merged";
	// ofstream mergedReadableKmer_ofs(mergedReadableKmerFile.c_str());
	// for(int tmpManipulateDecCode = 1; tmpManipulateDecCode <= totalGeneratedClassNum; tmpManipulateDecCode++)
	// {	
	// 	string tmp_kmer_readable_file = outputFolderStr + "final.class." + tmpManipulateDecCode + ".readable";
	// 	ifstream tmp_kmer_readable_ifs(tmp_kmer_readable_file.c_str());
	// 	while(!tmp_kmer_readable_ifs.eof())
	// 	{
	// 		string tmpStr;
	// 		getline(tmp_kmer_readable_ifs, tmpStr);
	// 		if(tmpStr == "")
	// 			break;
	// 		int tabLoc = tmpStr.find("\t");
	// 		string tmpKmerStr = tmpStr.substr(0, tabLoc);
	// 		mergedReadableKmer_ofs << tmpKmerStr << tmpManipulateDecCode << endl;
	// 	}
	// 	tmp_kmer_readable_ifs.close();
	// }
	// mergedReadableKmer_ofs.close();
	log_ofs.close();
	cmd_ofs.close();
	kmersetManipulate_ofs.close();
	return 0;
}