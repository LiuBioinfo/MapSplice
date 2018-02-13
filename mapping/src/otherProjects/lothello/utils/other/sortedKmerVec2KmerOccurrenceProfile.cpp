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
#include "../../../general/read_block_test.h"
#include "../../../general/otherFunc.h"
#include "../../../general/index_info.h"
#include "../general/sortedKmerVec2KmerOccurrenceProfile_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 7)
	{
		cout << "Executable#0 inputSortedKmerFile_union#1 outputKmerOccurrenceProfile#2" << endl;
		cout << "temporaryDir#3 withToAddIdOrNot_bool#4 sortedKmerFile_num#5" << endl;
		cout << "sortedKmerFile_1#6 toAddId_1#7 (sortedKmerFile_2#8 toAddId_2#9 ...)" << endl;
		cout << " or " << endl;
		cout << "sortedKmerFile_1#6 (sortedKmerFile_2#7 ...)" << endl;
		exit(1);
	}
	string inputSortedKmerFile_union = argv[1];
	string outputKmerOccurrenceProfile = argv[2];
	string temporaryDir = argv[3]; temporaryDir += "/";
	
	bool withToAddIdOrNot_bool;
	string withToAddIdOrNot_bool_str = argv[4];
	if(withToAddIdOrNot_bool_str == "Y")
		withToAddIdOrNot_bool = true;
	else if(withToAddIdOrNot_bool_str == "N")
		withToAddIdOrNot_bool = false;
	else
	{
		cout << "error! withToAddIdOrNot_bool_str should be Y or N" << endl;
		exit(1); 
	}

	string sortedKmerFile_num_str = argv[5];
	int sortedKmerFile_num = atoi(sortedKmerFile_num_str.c_str());
	vector<string> sortedKmerFileVec;
	vector<string> toAddIdVec;
	if(withToAddIdOrNot_bool)
	{
		for(int tmp = 0; tmp < sortedKmerFile_num; tmp++)
		{
			sortedKmerFileVec.push_back(6 + tmp * 2);
			toAddIdVec.push_back(7 + tmp * 2);
		}
	}
	else
	{
		for(int tmp = 0; tmp < sortedKmerFile_num; tmp++)
		{
			sortedKmerFileVec.push_back(6 + tmp);
			toAddIdVec.push_back(int_to_str(tmp));
		}
	}

	SortedKmerVec2KmerOccurrenceProfile_Info sortedKmerVec2KmerOccurrenceProfileInfo;
	sortedKmerVec2KmerOccurrenceProfileInfo.initiate_mkTmpDir(inputSortedKmerFile_union,
		sortedKmerFileVec, toAddIdVec, temporaryDir, outputKmerOccurrenceProfile);
	sortedKmerVec2KmerOccurrenceProfileInfo.generateKmerOccurrenceProfile();
	return 0;
}