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

int main(int argc, char** argv)
{
	if(argc != 10)
	{
		cout << "Executable graftIndex hostIndex graftAnn hostAnn graftSam_sortedByReadName ";
		cout << "hostSam_sortedByReadName outputFolder totalReadOrPairNum PE_or_SE" << endl;
		exit(1);
	}
	string graftIndexPath = argv[1];
	string hostIndexPath = argv[2];
	string graftAnnPath = argv[3];
	string hostAnnPath = argv[4];
	string graftSamPath = argv[5];
	string hostSamPath = argv[6];
	string outputFolderPath = argv[7];
	string totalReadOrPairNumStr = argv[8];
	string PE_or_SE_str = argv[9];

	int totalReadOrPairNum = atoi(totalReadOrPairNumStr.c_str());
	bool PE_or_SE_bool;
	if((PE_or_SE_str == "PE")||(PE_or_SE_str == "pe")||(PE_or_SE_str == "Pe"))
		PE_or_SE_bool = true;
	else if((PE_or_SE_str == "SE")||(PE_or_SE_str == "se")||(PE_or_SE_str == "Se"))
		PE_or_SE_bool = false;
	else
	{
		cout << "invalid PE_or_SE_str: " << PE_or_SE_str << endl;
		exit(1);
	}

	int readOrPairNumInEachTurn = 1000000;
	int tmpTurn = 0;




	return 0;
}