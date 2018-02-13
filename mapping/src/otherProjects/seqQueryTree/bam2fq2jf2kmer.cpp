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
#include <sstream>

using namespace std;

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

string getMinCovStr(string& filePath)
{
	long long fileSize = filesize(filePath.c_str());
	long long fileSize_M = fileSize/1048576;
	string tmp_mincov_str;
	if(fileSize_M <= 300)
		tmp_mincov_str = "1";
	else if(fileSize_M <= 500)
		tmp_mincov_str = "3";
	else if(fileSize_M <= 1024)
		tmp_mincov_str = "10";
	else if(fileSize_M <= 3072)
		tmp_mincov_str = "20";
	else
		tmp_mincov_str = "50";
	return tmp_mincov_str;
}

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 tophat2path (bam2fa)" << endl;
		cout << "#2 jellyfishPath (fa2jf2kmer)" << endl;
		cout << "#3 idList" << endl;
		cout << "#4 bamFilePathList" << endl;
		cout << "#5 outputDir" << endl;
		cout << "#6 threads_num" << endl;
		cout << "#7 Kmer_length" << endl;
		cout << "#8 bf_size" << endl;
		exit(1);
	}
	string bf_size_str = argv[8];
	string Kmer_length_str = argv[7];
	string threads_num_str = argv[6];


	string tophat2path = argv[1]; tophat2path += "/";
	string jellyfishPath = argv[2]; jellyfishPath += "/";
	string jellyFishBinPath = jellyfishPath;
	string idList = argv[3];
	string bamFilePathList = argv[4];
	string outputDir = argv[5]; outputDir += "/";

	vector<string> idVec;
	ifstream id_ifs(idList.c_str());
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		idVec.push_back(tmpStr);
	}
	id_ifs.close();

	vector<string> bamFilePathVec;
	ifstream bamFilePath_ifs(bamFilePathList.c_str());
	while(!bamFilePath_ifs.eof())
	{
		string tmpStr;
		getline(bamFilePath_ifs, tmpStr);
		if(tmpStr == "")
			break;
		bamFilePathVec.push_back(tmpStr);
	}
	bamFilePath_ifs.close();	

	if(idVec.size() != bamFilePathVec.size())
	{
		cout << "idVec.size() != bamFilePathVec.size()" << endl;
		cout << "idVec.size(): " << idVec.size() << endl;
		cout << "bamFilePathVec.size(): " << bamFilePathVec.size() << endl;
		exit(1);
	}

	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmp_id = idVec[tmp];
		string tmp_bam_path = bamFilePathVec[tmp];
		string tmp_fa_path = outputDir + tmp_id + ".fa";
		string tmp_jf_path = outputDir + tmp_id + ".jf";
		string tmp_kmer_path = outputDir + tmp_id + ".kmer";
		cout << "tmp_id: " << tmp_id << endl;
		// bam 2 fa
		string cmd_bam2fa = tophat2path + "bam2fastx --all --fasta -o " + tmp_fa_path + " " + tmp_bam_path; 
		cout << "start to do bam2fa: " << cmd_bam2fa << endl;
		system(cmd_bam2fa.c_str());
		// rm bam
		string cmd_rm_bam = "rm " + tmp_bam_path;
		cout << "start to do rm_bam: " << cmd_rm_bam << endl;
		system(cmd_rm_bam.c_str());
		// fa 2 jf
		cout << "file_size: " << filesize(tmp_fa_path.c_str()) << endl;
		string count_min_str = getMinCovStr(tmp_fa_path);
		cout << "count_min: " << count_min_str << endl;
		string cmd_fa2jf = jellyfishPath + "jellyfish count -o " + tmp_jf_path + " -m " 
			+ Kmer_length_str + " -t " + threads_num_str + " -s " + bf_size_str + " -C -L " 
			+ count_min_str + " " + tmp_fa_path;
		cout << "start to do fa2jf: " << cmd_fa2jf << endl;
		system(cmd_fa2jf.c_str());	
		// rm fa
		string cmd_rm_fa = "rm " + tmp_fa_path;
		cout << "start to do rm_fa: " << cmd_rm_fa << endl;
		system(cmd_rm_fa.c_str());
		// jf 2 kmer
		string cmd_jf2kmer = jellyFishBinPath + "jellyfish dump -t -c -o " + tmp_kmer_path + " " + tmp_jf_path;
		cout << "start to do jf2kmer: " << cmd_jf2kmer << endl;
		system(cmd_jf2kmer.c_str());	
		// rm jf
		string cmd_rm_jf = "rm " + tmp_jf_path;
		cout << "start to do rm_jf: " << cmd_rm_jf << endl << endl;
		system(cmd_rm_jf.c_str());
	
		//script_ofs << cmd_bam2fa << endl << cmd_rm_bam << endl 
		//	<< cmd_fa2jf << endl << cmd_rm_fa << endl
		//	<< cmd_jf2kmer << endl << cmd_rm_jf << endl << endl;
	}
	return 0;
}