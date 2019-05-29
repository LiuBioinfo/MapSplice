// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
/*    
 *    buildIndex.cpp
 *	  MapSplice3
 *
 *    Copyright (C) 2016 University of Kentucky and
 *                       Xinan Liu
 *
 *    Authors: Xinan Liu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputChromSeqFolderPath outputIndexFolder\nor\nExecutable mergedFaFile outputIndexFolder" << endl;
		exit(1);
	}
	
	//Check if argv[1] refers to a merged fa file or directory with individual files
	bool merged = false;
	struct stat buf;
	stat(argv[1], &buf);
	if (S_ISREG(buf.st_mode)) {
	    merged = true;
	}
	
	//Build from directory
	if (merged == false) {
        string inputChromSeqFolderPath = argv[1]; inputChromSeqFolderPath += "/";
        string outputIndexFolder = argv[2]; outputIndexFolder += "/";
	    string outputIndexFolder_2ndLevel = outputIndexFolder + "2ndLevelIndex";
	    
        string buildWholeGenomeIndex_cmd = "./buildWholeGenome " + inputChromSeqFolderPath 
		    + " " + outputIndexFolder;
	    cout << "buildWholeGenomeIndex_cmd: " << buildWholeGenomeIndex_cmd << endl;
	    string build2ndLevelIndex_cmd = "./build2ndLevelIndex " + outputIndexFolder
		+ " " + outputIndexFolder_2ndLevel;
	    cout << "build2ndLevelIndex_cmd: " << build2ndLevelIndex_cmd << endl;
	    
	    system(buildWholeGenomeIndex_cmd.c_str());
	    system(build2ndLevelIndex_cmd.c_str());
    }
    //Build from merged fa
    else {
        string mergedFaFile = argv[1];
        string outputIndexFolder = argv[2]; outputIndexFolder += "/";
	    string outputIndexFolder_2ndLevel = outputIndexFolder + "2ndLevelIndex";
	    
        string buildWholeGenomeIndex_mergedFa_cmd = "./buildWholeGenome_mergedFa " + mergedFaFile
		    + " " + outputIndexFolder;
	    cout << "buildWholeGenomeIndex_mergedFa_cmd: " << buildWholeGenomeIndex_mergedFa_cmd << endl;
	    string build2ndLevelIndex_cmd = "./build2ndLevelIndex " + outputIndexFolder
		+ " " + outputIndexFolder_2ndLevel;
	    cout << "build2ndLevelIndex_cmd: " << build2ndLevelIndex_cmd << endl;
	    
	    system(buildWholeGenomeIndex_mergedFa_cmd.c_str());
	    system(build2ndLevelIndex_cmd.c_str());
    }
	
	cout << "All done !" << endl;
	return 0;
}
