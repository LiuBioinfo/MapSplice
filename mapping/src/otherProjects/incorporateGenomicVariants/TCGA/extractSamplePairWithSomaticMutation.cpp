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
#include "../../../general/index_info.h"

using namespace std;

bool parseOutMafStr2Sample_BarCode_UUID(string& tmpMafStr, int& tmpChrNameInt, int& tmpChrPos, string& tmpAltBase_1, 
	string& tmpAltBase_2, string& tmpTumorSampleBarcode, string& tmpNormalSampleBarcode, 
	string& tmpTumorSampleUUID, string& tmpNormalSampleUUID, Index_Info* indexInfo)
{
	vector<string> mafFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 34; tmp++)
	{
		int tabLoc = tmpMafStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpMafStr find ......" << endl;
			exit(1);
		}
		string tmpMafField = tmpMafStr.substr(startLoc, tabLoc-startLoc);
		//cout << "tmpMafField: " << tmpMafField << "\t" << tmp+1 << endl;
		mafFieldVec.push_back(tmpMafField);
		startLoc = tabLoc + 1;
	}
	string tmpChrNameStr = "chr" + mafFieldVec[4];
	tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
	if(tmpChrNameInt < 0)
		return false;
	//cout << "tmpChrNameStr: " << tmpChrNameStr << endl;
	string tmpChrPosStr_start = mafFieldVec[5];
	string tmpChrPosStr_end = mafFieldVec[6];
	if(tmpChrPosStr_start != tmpChrPosStr_end) // not SNP
		return false;
	//cout << "tmpChrPosStr_start: " << tmpChrPosStr_start << endl;
	tmpChrPos = atoi(tmpChrPosStr_start.c_str());
	string tmpMutationType = mafFieldVec[9];
	if(tmpMutationType != "SNP")
		return false;
	//cout << "tmpMutationType: " << tmpMutationType << endl;
	string tmpRefBase = mafFieldVec[10];
	//cout << "tmpRefBase: " << tmpRefBase << endl;
	//cout << "tmpChrPos: " << tmpChrPos << endl;
	string tmpBaseInChrSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpChrPos, 1);
	//cout << "tmpBaseInChrSeq: " << tmpBaseInChrSeq << endl;
	if(tmpRefBase != tmpBaseInChrSeq)
	{
		cout << "tmpRefBase != tmpBaseInChrSeq" << endl;
		cout << "tmpMafStr: " << tmpMafStr << endl;
		exit(1); 
	}
	tmpAltBase_1 = mafFieldVec[11];
	tmpAltBase_2 = mafFieldVec[12];
	if((tmpRefBase == tmpAltBase_1)&&(tmpRefBase == tmpAltBase_2))
	{
		cout << "((tmpRefBase == tmpAltBase_1)&&(tmpRefBase == tmpAltBase_2))" << endl;
		exit(1);
	}
	if(tmpAltBase_1 == tmpAltBase_2)
	{
		cout << "tmpAltBase_1 == tmpAltBase_2" << endl;
		exit(1);
	}
	if((tmpAltBase_1.length() != 1)||(tmpAltBase_2.length() != 1))
	{
		cout << "(tmpAltBase_1.length() != 1)||(tmpAltBase_2.length() != 1)" << endl;
		cout << "tmpAltBase_1: " << tmpAltBase_1 << endl;
		cout << "tmpAltBase_2: " << tmpAltBase_2 << endl;
		exit(1);
	}
	tmpTumorSampleBarcode = mafFieldVec[15];
	tmpNormalSampleBarcode = mafFieldVec[16];
	tmpTumorSampleUUID = mafFieldVec[32];
	tmpNormalSampleUUID = mafFieldVec[33];
	return true;
}

void extractBarCodeUUIDvecFromMafFile(vector<string>& barcodeVec_T, vector<string>& barcodeVec_N,
	vector<string>& UUIDvec_T, vector<string>& UUIDvec_N, string& MAFfile, Index_Info* indexInfo)
{
	ifstream maf_ifs(MAFfile.c_str());
	string tmp1stLine;
	getline(maf_ifs, tmp1stLine);
	while(!maf_ifs.eof())
	{
		string tmpStr;
		getline(maf_ifs, tmpStr);
		//cout << "tmpStr: " << tmpStr << endl;
		if(tmpStr == "")
			break;
		int tmpChrNameInt, tmpChrPos;
		string tmpAltBase_1, tmpAltBase_2, tmpTumorSampleBarcode, tmpNormalSampleBarcode, 
			tmpTumorSampleUUID, tmpNormalSampleUUID;
		bool parseSuccess_bool = parseOutMafStr2Sample_BarCode_UUID(tmpStr, tmpChrNameInt, tmpChrPos, 
			tmpAltBase_1, tmpAltBase_2, tmpTumorSampleBarcode, tmpNormalSampleBarcode, 
			tmpTumorSampleUUID, tmpNormalSampleUUID, indexInfo);
		if(parseSuccess_bool)
		{	
			barcodeVec_T.push_back(tmpTumorSampleBarcode);
			barcodeVec_N.push_back(tmpNormalSampleBarcode);
			UUIDvec_T.push_back(tmpTumorSampleUUID);
			UUIDvec_N.push_back(tmpNormalSampleUUID);
		}
	}
	maf_ifs.close();
}

void generateUniqueVecFromRawVec(
	vector<string>& barcodeVec_T_raw, vector<string>& barcodeVec_N_raw, 
	vector<string>& UUIDvec_T_raw, vector<string>& UUIDvec_N_raw, 
	vector<string>& barcodeVec_uniq, vector<string>& UUIDvec_uniq)
{
	vector<string> barcodeVec_raw;
	vector<string> UUIDvec_raw;

	int rawVecSize_T = barcodeVec_T_raw.size();
	for(int tmp = 0; tmp < rawVecSize_T; tmp++)
	{
		string tmpBarcode_T_raw = barcodeVec_T_raw[tmp];
		string tmpBarcode_N_raw = barcodeVec_N_raw[tmp];
		string tmpUIUD_T_raw = UUIDvec_T_raw[tmp];
		string tmpUIUD_N_raw = UUIDvec_N_raw[tmp];
		barcodeVec_raw.push_back(tmpBarcode_T_raw);
		barcodeVec_raw.push_back(tmpBarcode_N_raw);
		UUIDvec_raw.push_back(tmpUIUD_T_raw);
		UUIDvec_raw.push_back(tmpUIUD_N_raw);
	}
	int rawVecSize = barcodeVec_raw.size();
	for(int tmp1 = 0; tmp1 < rawVecSize; tmp1++)
	{
		string tmpBarcode_raw = barcodeVec_raw[tmp1];		
		string tmpUIUD_raw = UUIDvec_raw[tmp1];		
		int currentUniqVecSize = barcodeVec_uniq.size();
		bool new_bool = true;
		for(int tmp2 = 0; tmp2 < currentUniqVecSize; tmp2 ++)
		{
			string tmpBarcode_uniq = barcodeVec_uniq[tmp2];		
			string tmpUIUD_uniq = UUIDvec_uniq[tmp2];
			if((tmpBarcode_raw == tmpBarcode_uniq)
				&&(tmpUIUD_raw == tmpUIUD_uniq))
			{
				new_bool = false;
				break;
			}
		}
		if(new_bool)
		{
			barcodeVec_uniq.push_back(tmpBarcode_raw);
			UUIDvec_uniq.push_back(tmpUIUD_raw);
		}
		else
		{}
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputTCGAmafFile outputFile" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	cout << "start to extractBarCodeUUIDvecFromMafFile ..." << endl;
	string TCGAmafFile = argv[2];
	vector<string> barcodeVec_T_raw;
	vector<string> barcodeVec_N_raw;
	vector<string> UUIDvec_T_raw;
	vector<string> UUIDvec_N_raw;	
	extractBarCodeUUIDvecFromMafFile(barcodeVec_T_raw, barcodeVec_N_raw,
		UUIDvec_T_raw, UUIDvec_N_raw, TCGAmafFile, indexInfo);

	cout << "start to generate unique BarCodeUUIDvec From raw vec ..." << endl;
	vector<string> barcodeVec_uniq;
	vector<string> UUIDvec_uniq;
	generateUniqueVecFromRawVec(barcodeVec_T_raw, barcodeVec_N_raw, 
		UUIDvec_T_raw, UUIDvec_N_raw, barcodeVec_uniq, UUIDvec_uniq);

	cout << "start to output unique BarCodeUUIDvec ..." << endl;
	int barcodeVec_uniq_size = barcodeVec_uniq.size();
	string pairedSampleIdFile = argv[3];
	ofstream pairedSampleId_ofs(pairedSampleIdFile.c_str());
	for(int tmp = 0; tmp < barcodeVec_uniq_size; tmp++)
		pairedSampleId_ofs << barcodeVec_uniq[tmp] << "\t" << UUIDvec_uniq[tmp] << endl;

	pairedSampleId_ofs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}