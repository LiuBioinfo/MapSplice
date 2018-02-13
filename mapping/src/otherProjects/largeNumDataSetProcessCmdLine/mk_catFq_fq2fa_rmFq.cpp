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
	if(argc != 5)
	{
		cout << "Executable dataListFile rawDatFolderPath reformattedDataFolderPath outputScript" << endl;
		exit(1);
	}
	string dataListFile = argv[1];
	string rawDataFolderaPath = argv[2];
	string reformattedDataFolderPath = argv[3];
	string outputScript = argv[4];
	ofstream script_ofs(outputScript.c_str());
	ifstream dataList_ifs(dataListFile.c_str());
	while(!dataList_ifs.eof())
	{
		string tmpStr;
		getline(dataList_ifs, tmpStr);
		if(tmpStr == "")
			break;
		//dataIdVec.push_back(tmpStr);
		string tmpRawReadPath = rawDataFolderaPath + "/" + tmpStr + "/*";
		string tmp_cat_cmd = "cat " + tmpRawReadPath + " > " + reformattedDataFolderPath + "/" + tmpStr + ".fq"; 
		string tmp_fq2fa_cmd = "./fq2fa " + reformattedDataFolderPath + "/" + tmpStr + ".fq " 
			+ reformattedDataFolderPath + "/" + tmpStr + ".fa";
		string tmp_rmFq_cmd = "rm " + reformattedDataFolderPath + "/" + tmpStr + ".fq";
		script_ofs << tmp_cat_cmd << endl << tmp_fq2fa_cmd << endl << tmp_rmFq_cmd << endl;
	}
	dataList_ifs.close();
	script_ofs.close();
	return 0;
}