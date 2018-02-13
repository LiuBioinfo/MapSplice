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
#include "general/bloom_filter.hpp"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4){
		cout << "Executable inputFa2buildBF inputQueryFa outputFolder" << endl;
		exit(1);
	}
	string inputFa2buildBF = argv[1];
	string inputQueryFa = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogFile = outputFolderStr + "log.txt";
	ofstream log_ofs(outputLogFile.c_str());

    nowtime = time(NULL);
    local = localtime(&nowtime);
    cout << endl << "[" << asctime(local) << "start to initiate parameters" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to initiate parameters" << endl;
	bloom_parameters parameters;
	// How many elements roughly do we expect to insert?
	parameters.projected_element_count = 3000000000;
	// Maximum tolerable false positive probability? (0,1)
	parameters.false_positive_probability = 0.0001; // 1 in 10000
	// Simple randomizer (optional)
	parameters.random_seed = 0xA5A5A5A5;
	if(!parameters){
    	std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
    	return 1;
    }
	else
		parameters.compute_optimal_parameters();
    
    nowtime = time(NULL);
    local = localtime(&nowtime);
    cout << endl << "[" << asctime(local) << "start to initiate bloom filter" << endl;
	log_ofs << endl << "[" << asctime(local) << "start to initiate bloom filter" << endl;
    //Instantiate Bloom Filter
    bloom_filter filter(parameters);
    filter.buildFromFaFile(inputFa2buildBF, log_ofs);
    
    nowtime = time(NULL);
    local = localtime(&nowtime);
    cout << endl << "[" << asctime(local) << "start to do query" << endl;
    log_ofs << endl << "[" << asctime(local) << "start to do query" << endl;
    // query against filter
    string hit_query_file = outputFolderStr + "hit_query.txt";
    string miss_query_file = outputFolderStr + "miss_query.txt";
    ofstream hit_ofs(hit_query_file.c_str());
    ofstream miss_ofs(miss_query_file.c_str());
    ifstream query_ifs(inputQueryFa.c_str());
    int tmpLineNO = 0;
    while(!query_ifs.eof())
    {
    	string tmpName;
    	getline(query_ifs, tmpName);
        tmpLineNO ++;
        int tmpThousandIndex = tmpLineNO / 1000000;
        if(tmpLineNO == tmpThousandIndex * 1000000)
            cout << "Processed Line #: " << tmpLineNO << endl;
    	if(tmpName == "")
    		break;
        string tmpSeq;
        getline(query_ifs, tmpSeq);
    	if(filter.contains(tmpSeq))
    		hit_ofs << tmpName << endl << tmpSeq << endl;
    	else
    		miss_ofs << tmpName << endl << tmpSeq << endl;
    }
    nowtime = time(NULL);
    local = localtime(&nowtime);
    cout << endl << "[" << asctime(local) << "end of query" << endl;
    log_ofs << endl << "[" << asctime(local) << "end of query" << endl;
    query_ifs.close();
    miss_ofs.close();
    hit_ofs.close();
    log_ofs.close();
    return 0;
}