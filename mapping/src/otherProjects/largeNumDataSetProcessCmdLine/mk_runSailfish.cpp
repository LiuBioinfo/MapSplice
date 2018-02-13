// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
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

string getSRRidFromFqPath(string& fq_path)
{
	// /scratch/lji226/projects/xal_seq/fq/SRR988510.fastq
	int fqSlashLoc = fq_path.find("fq/SRR");
	int dotFastqLoc = fq_path.find(".fastq");
	string SRRid = fq_path.substr(fqSlashLoc + 3, dotFastqLoc - 1 - fqSlashLoc - 3 + 1);
	return SRRid;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 fastq_list_file" << endl;
		cout << "#2 outputDir" << endl;
		cout << "#3 scriptFile" << endl;
		exit(1);
	}

	string fastq_list_file = argv[1];
	string outputDir = argv[2]; outputDir += "/";
	string scriptFile = argv[3];

	vector<string> fqVec;
	vector<string> srrVec;

	ifstream fq_ifs(fastq_list_file.c_str());
	while(!fq_ifs.eof())
	{
		string tmpStr;
		getline(fq_ifs, tmpStr);
		if(tmpStr == "")
			break;
		fqVec.push_back(tmpStr);
		string tmpSRRid = getSRRidFromFqPath(tmpStr);
		srrVec.push_back(tmpSRRid);
	}
	fq_ifs.close();

	ofstream script_ofs(scriptFile.c_str());
	script_ofs << "export DYLD_LIBRARY_PATH=/home/xli262/tools/SailfishBeta-0.9.0_DebianSqueeze/lib:$DYLD_LIBRARY_PATH" << endl;
	script_ofs << "export PATH=/home/xli262/tools/SailfishBeta-0.9.0_DebianSqueeze/bin:$PATH" << endl << endl;
	script_ofs << "cd /home/xli262/tools/SailfishBeta-0.9.0_DebianSqueeze/bin/" << endl << endl;
	for(int tmp = 0; tmp < fqVec.size(); tmp++)
	{
		string tmp_fq_path = fqVec[tmp];
		string tmp_srr_id = srrVec[tmp];
		string tmp_output_dir = outputDir + tmp_srr_id;
		script_ofs << "mkdir " << tmp_output_dir << endl;
		script_ofs << "./sailfish quant -i /home/xli262/chrom_Index/sailfish_index/human/ -l U -p 16 \\" << endl 
			<< "-r " << tmp_fq_path << " \\" << endl << "-o " << tmp_output_dir << endl << endl;    
	}
	script_ofs.close();
	return 0;
}