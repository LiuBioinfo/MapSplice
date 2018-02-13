#ifndef BATCH_MANAGER_INFO_H
#define BATCH_MANAGER_INFO_H

#include <string>
#include <string.h>
#include <vector>

using namespace std;

class Batch_Manager_Info
{
private:
	vector<string> SNP_file_vec;
	vector<string> SNPmer_index_folder_vec;

	vector<int> totalReadNumVec;

	vector<bool> SNPlearnedSuccess_bool_vec;
	vector<string> read_file_1_vec;
	vector<string> read_file_2_vec;
	vector<string> result_folder_vec;
	int data_set_num;
public:
	vector<Stats_Info*> statsInfoVec;

	Batch_Manager_Info()
	{
		data_set_num = 0;
	}

	void initaite_readFileVec_withReadFileInCommandLine(string& InputReadFile, string& InputReadFile_PE)
	{
		read_file_1_vec.push_back(InputReadFile);
		read_file_2_vec.push_back(InputReadFile_PE);
	}

	void assign_SNPlearnedSuccessBool(bool tmpBool, int tmpData)
	{
		SNPlearnedSuccess_bool_vec[tmpData] = tmpBool;
	}

	void initiateStatsInfoVec_withDataSetNum(bool SE_or_PE_bool, int threads_num)
	{
		if(SE_or_PE_bool)
		{
			for(int tmp = 0; tmp < data_set_num; tmp++)
			{
				Stats_Info* tmpData_statsInfo = new Stats_Info();
				tmpData_statsInfo->initiate_stats_info_SE(threads_num);
				statsInfoVec.push_back(tmpData_statsInfo);
			}
		}
		else
		{
			for(int tmp = 0; tmp < data_set_num; tmp++)
			{
				Stats_Info* tmpData_statsInfo = new Stats_Info();
				tmpData_statsInfo->initiate_stats_info_PE(threads_num);
				statsInfoVec.push_back(tmpData_statsInfo);
			}
		}
	}

	bool returnSNPlearnedSuccessBool(int tmp)
	{
		return SNPlearnedSuccess_bool_vec[tmp];
	}

	int return_dataSet_num()
	{
		return data_set_num;
	}

	string returnSNPfile_withIndexInVec(int tmp)
	{
		return SNP_file_vec[tmp];
	}

	string returnSNPmerIndexFolderPath_withIndexInVec(int tmp)
	{
		return SNPmer_index_folder_vec[tmp];
	}

	void initiate_SNPfileVec_withSNPfileListFile(
		string& tmpSNPfileListFile, string& tmpSNPmerIndexFolderListFile)
	{
		ifstream SNP_ifs(tmpSNPfileListFile.c_str());
		while(!SNP_ifs.eof())
		{
			string tmpFileStr;
			getline(SNP_ifs, tmpFileStr);
			if(tmpFileStr == "")
				break;
			SNP_file_vec.push_back(tmpFileStr);				
		}
		SNP_ifs.close();
		ifstream SNPmer_index_ifs(tmpSNPmerIndexFolderListFile.c_str());
		while(!SNPmer_index_ifs.eof())
		{
			string tmpIndexPathStr;
			getline(SNPmer_index_ifs, tmpIndexPathStr);
			if(tmpIndexPathStr == "")
				break;
			SNPmer_index_folder_vec.push_back(tmpIndexPathStr);
		}
		SNPmer_index_ifs.close();
	}

	void initiate_SNPlearnedSuccessBoolVec_withDataSetNum()
	{
		for(int tmp = 0; tmp < data_set_num; tmp++)
		{
			bool tmp_SNPlearned_success_bool = false;
			SNPlearnedSuccess_bool_vec.push_back(tmp_SNPlearned_success_bool);
		}
	}

	void assignTotalReadNum(int tmp_total_read_num, int tmpData)
	{
		totalReadNumVec[tmpData] = tmp_total_read_num;
	}

	void initaite_totalReadNumVec_withDataSetNum()
	{
		for(int tmp = 0; tmp < data_set_num; tmp++)
		{
			totalReadNumVec.push_back(0);
		}
	}

	void initiate_readFileVec_withReadFileListFile(
		string& tmpReadFileListFile_1, string& tmpReadFileListFile_2)
	{
		ifstream fileList_1_ifs(tmpReadFileListFile_1.c_str());
		ifstream fileList_2_ifs(tmpReadFileListFile_2.c_str());
		while(!fileList_1_ifs.eof())
		{
			string tmpFileStr;
			getline(fileList_1_ifs, tmpFileStr);
			if(tmpFileStr == "")
				break;
			read_file_1_vec.push_back(tmpFileStr);
		}
		while(!fileList_2_ifs.eof())
		{
			string tmpFileStr;
			getline(fileList_2_ifs, tmpFileStr);
			if(tmpFileStr == "")
				break;
			read_file_2_vec.push_back(tmpFileStr);
		}		
		int read_file_1_vec_size = read_file_1_vec.size();
		int read_file_2_vec_size = read_file_2_vec.size();
		if(read_file_1_vec_size != read_file_2_vec_size)
		{
			cout << "read_file_1_vec_size != read_file_2_vec_size" << endl;
			cout << "read_file_1_vec_size: " << read_file_1_vec_size << endl;
			cout << "read_file_2_vec_size: " << read_file_2_vec_size << endl;
			exit(1); 
		}
		fileList_1_ifs.close();
		fileList_2_ifs.close();
	}

	void initiate_dataSetNum_accordingToReadFileVec()
	{
		data_set_num = read_file_1_vec.size();
	}

	void initiate_create_resultFolderVec_accordingToReadFileVec(string& tmpBatchResultTotalFolder)
	{
		string batchResultDir = tmpBatchResultTotalFolder;
 	  	string mkdirBatchResultDirCommand = "mkdir -p " + batchResultDir;
   		system(mkdirBatchResultDirCommand.c_str());
   		string readFile2folder_info_file_str = batchResultDir + "/read_file_2_result_folder_info.txt";
   		ofstream readFile2folder_info_ofs(readFile2folder_info_file_str.c_str());
   		for(int tmp = 0; tmp < data_set_num; tmp++)
   		{
   			string tmpReadFile_1 = read_file_1_vec[tmp];
   			string tmpReadFile_2 = read_file_2_vec[tmp];
   			string tmpFolderStr = batchResultDir + "/" + int_to_str(tmp + 1);
   			string mkdirTmpDataSetDirCommand = "mkdir -p " + tmpFolderStr;
   			system(mkdirTmpDataSetDirCommand.c_str());
   			readFile2folder_info_ofs << "tmp_data_set: " << tmp + 1 << endl; 
   			readFile2folder_info_ofs << "tmp_read_file:" << endl 
   				<< tmpReadFile_1 << endl << tmpReadFile_2 << endl;
   			readFile2folder_info_ofs << "tmp_folder: " << endl << tmpFolderStr << endl;
   			result_folder_vec.push_back(tmpFolderStr);
   		}
   		readFile2folder_info_ofs.close();
	}

	int returnTotalReadNum_withIndexInVec(int tmp)
	{
		return totalReadNumVec[tmp];
	}

	string returnReadFile_1_withIndexInVec(int tmp)
	{
		return read_file_1_vec[tmp];
	}

	string returnReadFile_2_withIndexInVec(int tmp)
	{
		return read_file_2_vec[tmp];
	}

	string returnResultFolder_withIndexInVec(int tmp)
	{
		return result_folder_vec[tmp];
	}
};
#endif