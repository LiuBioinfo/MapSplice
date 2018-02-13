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
	if(argc != 4)
	{
		cout << "Executable localIndexFolderPath localIndexNum outputScript" << endl;
		exit(1);
	}
	string localIndexFolderPath = argv[1];
	string localIndexNumStr = argv[2];
	int localIndexNum = atoi(localIndexNumStr.c_str());
	string outputScript = argv[3];
	ofstream script_ofs(outputScript.c_str());
	for(int tmp = 1; tmp <= localIndexNum; tmp++)
	{
		script_ofs << "rm " << localIndexFolderPath << "/" << tmp << "/lcp" << endl;
		script_ofs << "rm " << localIndexFolderPath << "/" << tmp << "/up" << endl;
		script_ofs << "rm " << localIndexFolderPath << "/" << tmp << "/down" << endl;
		script_ofs << "rm " << localIndexFolderPath << "/" << tmp << "/next" << endl;
	}
	script_ofs.close();
	return 0;
}