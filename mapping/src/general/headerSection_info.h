// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef HEADERSECTION_INFO_H
#define HEADERSECTION_INFO_H

#include <stdlib.h>
#include <string>
#include <string.h>
//temporary
#include <fstream>
#include <iostream>

using namespace std;

class HeaderSection_Info
{
private:
	string SQ_str;
	string headerSectionInfoStr;
public:
	string returnHeaderSectionInfoStr()
	{
		headerSectionInfoStr = SQ_str;
		return headerSectionInfoStr;
	}

	HeaderSection_Info()
	{
		SQ_str = "@SQ\tSN:Not_avaliable";
		//headerSectionInfoStr = 
	}

	HeaderSection_Info(bool dict_provided_bool, Index_Info* indexInfo, string dict_file_path//, Option_Info* optionInfo
		)
	{
        //No dict provided, create unordered header
        if (!dict_provided_bool) {
            int chromNum = indexInfo->returnChrNameStrSize();
		    SQ_str = "@SQ\tSN:" + indexInfo->returnChrNameStr(0)
			    + "\tLN:" + int_to_str(indexInfo->returnChromLength(0)) + "\n";
		    for(int tmp = 1; tmp < (
			    chromNum - 1
			    ); tmp++)
		    {
			    SQ_str = SQ_str
				    + "@SQ\tSN:" + indexInfo->returnChrNameStr(tmp)
				    + "\tLN:" + int_to_str(indexInfo->returnChromLength(tmp)) + "\n"; 
		    }
		    SQ_str = SQ_str
			    + "@SQ\tSN:" + indexInfo->returnChrNameStr(chromNum - 1)
			    + "\tLN:" + int_to_str(indexInfo->returnChromLength(chromNum - 1)-1);
        }
        //Dict provided, order header according to dict
		else {
            ifstream dict_ifs;
            dict_ifs.open(dict_file_path);
            if (dict_ifs.fail()) {
                cout << "Failed to open dict file, make sure file path is correct. Terminating..." << endl;
                exit(1);
            }
            string temp;
            SQ_str = "";
            while (!dict_ifs.eof()) {
                dict_ifs >> temp;
                if (temp == "@SQ") {
                    SQ_str += temp + '\t'; //SQ
                    dict_ifs >> temp;
                    SQ_str += temp + '\t'; //SN
                    dict_ifs >> temp;
                    SQ_str += temp + '\n'; //LN
                }
            }
            SQ_str.erase(SQ_str.length() - 1, 1); //remove extra \n at end
            dict_ifs.close();
        }
	}

};

#endif
