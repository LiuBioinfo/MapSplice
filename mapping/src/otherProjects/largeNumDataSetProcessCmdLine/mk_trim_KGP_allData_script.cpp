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

void getIdVec(string& dataListFile, vector<string>& idVec)
{
	ifstream id_ifs(dataListFile.c_str());
	while(!id_ifs.eof())
	{
		string tmpStr;
		getline(id_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		if(tabLoc == string::npos)
			idVec.push_back(tmpStr);
		else
			idVec.push_back(tmpStr.substr(0, tabLoc));
	}
	id_ifs.close();
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 dataList" << endl;
		cout << "#1 outputScript" << endl;
		exit(1);
	}
	string dataList = argv[1];
	string outputScript = argv[2];
	ofstream script_ofs(outputScript.c_str());
	script_ofs << "cd /home/xli262/yizhang/Trimmomatic-0.33/" << endl;
	vector<string> idVec;
	getIdVec(dataList, idVec);
	for(int tmp = 0; tmp < idVec.size(); tmp++)
	{
		string tmpId = idVec[tmp];
		script_ofs << endl << "mkdir -p /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << endl;  
		script_ofs << "java -jar /home/xli262/yizhang/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 16 \\" << endl;
		script_ofs << "/scratch/xli262/impsOn1000genome/CEU_462_rawData/data/" << tmpId << "/" << tmpId << "_1.fastq \\" << endl;
		script_ofs << "/scratch/xli262/impsOn1000genome/CEU_462_rawData/data/" << tmpId << "/" << tmpId << "_2.fastq \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/1_paired.fastq \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/1_unpaired.fastq \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/2_paired.fastq \\" << endl;
		script_ofs << "/scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpId << "/2_unpaired.fastq \\" << endl;		
		script_ofs << "ILLUMINACLIP:/home/xli262/yizhang/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa:2:30:10 MINLEN:50" << endl;
	}
	script_ofs.close();
	return 0;
}