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
	if(argc != 7)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 cutoff" << endl;
		cout << "#2 threads" << endl;
		cout << "#3 hashfile" << endl;
		cout << "#4 bitmap_size" << endl;
		cout << "#5 fq_bfbv_list_file" << endl;
		cout << "#6 output_script" << endl;
		exit(1);
	}
	string cutoff_str = argv[1];
	string threads_str = argv[2];
	string hashfile = argv[3];
	string bitmap_size_str = argv[4];
	string fq_bfbv_list_file = argv[5];
	string output_script = argv[6];

	vector<string> fq_vec;
	vector<string> bfbv_vec;
	ifstream fq_bfbv_list_ifs(fq_bfbv_list_file.c_str());
	while(!fq_bfbv_list_ifs.eof())
	{
		string tmpStr;
		getline(fq_bfbv_list_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc = tmpStr.find("\t");
		string tmpFqFile = tmpStr.substr(0, tabLoc);
		string tmpBfbvFile = tmpStr.substr(tabLoc + 1);
		fq_vec.push_back(tmpFqFile);
		bfbv_vec.push_back(tmpBfbvFile);
	}
	fq_bfbv_list_ifs.close();

	ofstream script_ofs(output_script.c_str());
	for(int tmp = 0; tmp < fq_vec.size(); tmp++)
		script_ofs << "./bt count --cutoff " << cutoff_str << " --threads " 
			<< threads_str << " " << hashfile << " " << bitmap_size_str 
			<< " " << fq_vec[tmp] << " " << bfbv_vec[tmp] << endl;
	script_ofs.close();
	return 0;
}