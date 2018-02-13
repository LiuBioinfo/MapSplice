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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "Executable folderPathToCreate" << endl;
		exit(1);
	}
	string folderPathToCreate = argv[1];
	string mkdir_cmd = "mkdir -p " + folderPathToCreate;
	//system(mkdir_cmd);
	system(mkdir_cmd.c_str());
	return 0;
}