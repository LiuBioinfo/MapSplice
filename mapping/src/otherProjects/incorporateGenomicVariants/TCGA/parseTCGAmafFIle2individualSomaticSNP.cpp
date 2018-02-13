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
#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
using namespace std;

bool parseOutMafStr2SNP(string& tmpMafStr, int& tmpChrNameInt, int& tmpChrPos, string& tmpAltBase_1, 
	string& tmpAltBase_2, string& tmpTumorSampleId, string& tmpNormalSampleId, Index_Info* indexInfo)
{
	vector<string> mafFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 17; tmp++)
	{
		int tabLoc = tmpMafStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpMafStr find ......" << endl;
			exit(1);
		}
		string tmpMafField = tmpMafStr.substr(startLoc, tabLoc-startLoc);
		mafFieldVec.push_back(tmpMafField);
		startLoc = tabLoc + 1;
	}
	string tmpChrNameStr = "chr" + mafFieldVec[4];
	tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
	if(tmpChrNameInt < 0)
		return false;
	string tmpChrPosStr_start = mafFieldVec[5];
	string tmpChrPosStr_end = mafFieldVec[6];
	if(tmpChrPosStr_start != tmpChrPosStr_end) // not SNP
		return false;
	tmpChrPos = atoi(tmpChrPosStr_start.c_str());
	string tmpMutationType = mafFieldVec[9];
	if(tmpMutationType != "SNP")
		return false;
	string tmpRefBase = mafFieldVec[10];
	string tmpBaseInChrSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt, tmpChrPos, 1);
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
	tmpTumorSampleId = mafFieldVec[15];
	tmpNormalSampleId = mafFieldVec[16];
	return true;
}

void parseMafFile2SNPvec(string& TCGAmafFile, vector<int>& chrNameIntVec, vector<int>& chrPosVec, vector<string>& altBaseVec_1, 
	vector<string>& altBaseVec_2, vector<string>& tumorSampeIdVec, vector<string>& normalSampleIdVec, Index_Info* indexInfo,
	string& simplifiedSNPfile)
{
	ofstream simplifiedSNP_ofs(simplifiedSNPfile.c_str());
	ifstream TCGAmaf_ifs(TCGAmafFile.c_str());
	string tmp1stLine;
	getline(TCGAmaf_ifs, tmp1stLine);
	while(!TCGAmaf_ifs.eof())
	{
		string tmpStr;
		getline(TCGAmaf_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tmpChrNameInt, tmpChrPos; 
		string tmpAltBase_1, tmpAltBase_2, tmpTumorSampleId, tmpNormalSampleId;
		bool parseOutMafStr2SNP_bool = parseOutMafStr2SNP(tmpStr, tmpChrNameInt, tmpChrPos, 
			tmpAltBase_1, tmpAltBase_2, tmpTumorSampleId, tmpNormalSampleId, indexInfo);
		if(parseOutMafStr2SNP_bool)
		{
			chrNameIntVec.push_back(tmpChrNameInt);
			chrPosVec.push_back(tmpChrPos);
			altBaseVec_1.push_back(tmpAltBase_1);
			altBaseVec_2.push_back(tmpAltBase_2);
			tumorSampeIdVec.push_back(tmpTumorSampleId);
			normalSampleIdVec.push_back(tmpNormalSampleId);
			simplifiedSNP_ofs << indexInfo->returnChrNameStr(tmpChrNameInt) << "\t" << tmpChrPos << "\t" << tmpAltBase_1 
				<< "\t" << tmpAltBase_2 << "\t" << tmpTumorSampleId << "\t" << tmpNormalSampleId << endl;
		}
	}
	TCGAmaf_ifs.close();
	simplifiedSNP_ofs.close();
}	

int tumorNormalSampleId2matchedSamplePairId(string& tmpTumorSampleId, string& tmpNormalSampleId)
{
	vector<string> tumorIdFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpTumorSampleId.find("-", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpTumorSampleId find ......" << endl;
			exit(1);
		}
		string tumorIdField = tmpTumorSampleId.substr(startLoc, tabLoc-startLoc);
		tumorIdFieldVec.push_back(tumorIdField);
		startLoc = tabLoc + 1;
	}
	vector<string> normalIdFieldVec;
	startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpNormalSampleId.find("-", startLoc);
		if(tabLoc == string::npos)
		{
			cout << "invalid tmpTumorSampleId find ......" << endl;
			exit(1);
		}
		string normalIdField = tmpNormalSampleId.substr(startLoc, tabLoc-startLoc);
		normalIdFieldVec.push_back(normalIdField);
		startLoc = tabLoc + 1;
	}	
	string tmpTumorRedefinedIdStr = tumorIdFieldVec[2];
	string tmpNormalRedefinedIdStr = normalIdFieldVec[2];
	if(tmpTumorRedefinedIdStr != tmpNormalRedefinedIdStr)
		return -1;
	else
		return atoi(tmpTumorRedefinedIdStr.c_str());
}

void snpVec_2_matchedSamplePairId2snpSetVec(vector< pair<int, set<int> > >& matchedSamplePairId2snpSetVec, 
	vector<int>& chrNameIntVec, vector<int>& chrPosVec, vector<string>& altBaseVec_1, 
	vector<string>& altBaseVec_2, vector<string>& tumorSampeIdVec, vector<string>& normalSampleIdVec)
{
	int SNPnum = chrNameIntVec.size();
	for(int tmp = 0; tmp < SNPnum; tmp++)
	{
		string tmpTumorSampleId = tumorSampeIdVec[tmp];
		string tmpNormalSampleId = normalSampleIdVec[tmp];
		int tmpMatchedSamplePairId = tumorNormalSampleId2matchedSamplePairId(tmpTumorSampleId, tmpNormalSampleId);
		if(tmpMatchedSamplePairId < 0)
		{
			cout << "inconsistent tumorSampeId and normalSampleId" << endl;
			cout << "tmpTumorSampleId: " << tmpTumorSampleId << endl;
			cout << "tmpNormalSampleId: " << tmpNormalSampleId << endl;
			continue;
		}
		int matchedSamplePairId2snpSetVecSize = matchedSamplePairId2snpSetVec.size();
		bool existInCurrentVec_bool = false;
		for(int tmp2 = 0; tmp2 < matchedSamplePairId2snpSetVecSize; tmp2 ++)
		{
			int tmpExistMatchedSamplePairId = matchedSamplePairId2snpSetVec[tmp2].first;
			if(tmpExistMatchedSamplePairId == tmpMatchedSamplePairId)
			{
				(matchedSamplePairId2snpSetVec[tmp2].second).insert(tmp);
				existInCurrentVec_bool = true;
				break;
			}
		}
		if(!existInCurrentVec_bool)
		{
			set<int> tmpNewSet;
			tmpNewSet.insert(tmp);
			matchedSamplePairId2snpSetVec.push_back(pair<int, set<int> >(tmpMatchedSamplePairId, tmpNewSet));
		}
		else
		{}
	}

}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath TCGAmafFile outputFolder" << endl;
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

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());	

	cout << "start to parseMafFile2SNPvec ..." << endl;
	log_ofs << "start to parseMafFile2SNPvec ..." << endl;
	string TCGAmafFile = argv[2];
	string simplifiedSNPfile = outputFolderStr + "SNP.simplified.txt";
	vector<int> chrNameIntVec;
	vector<int> chrPosVec;
	vector<string> altBaseVec_1;
	vector<string> altBaseVec_2;
	vector<string> tumorSampeIdVec;
	vector<string> normalSampleIdVec;
	parseMafFile2SNPvec(TCGAmafFile, chrNameIntVec, chrPosVec, altBaseVec_1, altBaseVec_2, 
		tumorSampeIdVec, normalSampleIdVec, indexInfo, simplifiedSNPfile);

	cout << "start to do snpVec_2_matchedSamplePairId2snpSetVec ..." << endl;
	log_ofs << "start to do snpVec_2_matchedSamplePairId2snpSetVec ..." << endl;
	vector< pair<int, set<int> > > matchedSamplePairId2snpSetVec;
	snpVec_2_matchedSamplePairId2snpSetVec(matchedSamplePairId2snpSetVec, chrNameIntVec, 
		chrPosVec, altBaseVec_1, altBaseVec_2, tumorSampeIdVec, normalSampleIdVec);

	cout << "start to output somatic mutations for each samplePair ..." << endl;
	log_ofs << "start to output somatic mutations for each samplePair ..." << endl;
	int samplePairNum = matchedSamplePairId2snpSetVec.size();
	cout << "samplePairNum: " << samplePairNum << endl;
	log_ofs << "samplePairNum: " << samplePairNum << endl;
	for(int tmpSamplePair = 0; tmpSamplePair < samplePairNum; tmpSamplePair ++)
	{
		int tmpSamplePairId = matchedSamplePairId2snpSetVec[tmpSamplePair].first;
		string tmpSamplePairId_str = int_to_str(tmpSamplePairId);		
		string tmpSamplePair_snp_file_1 = outputFolderStr + tmpSamplePairId_str + "_SNP.1.txt";
		string tmpSamplePair_snp_file_2 = outputFolderStr + tmpSamplePairId_str + "_SNP.2.txt";
		ofstream SNP_1_ofs(tmpSamplePair_snp_file_1.c_str());
		ofstream SNP_2_ofs(tmpSamplePair_snp_file_2.c_str());
		for(set<int>::iterator tmpSetIter = (matchedSamplePairId2snpSetVec[tmpSamplePair].second).begin(); 
			tmpSetIter != (matchedSamplePairId2snpSetVec[tmpSamplePair].second).end(); tmpSetIter ++)
		{
			int tmpIndexInSNPvec = (*tmpSetIter);
			int tmpSNP_chrNameInt = chrNameIntVec[tmpIndexInSNPvec];
			string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
			int tmpSNP_chrPos = chrPosVec[tmpIndexInSNPvec];
			string tmpSNP_refBase = indexInfo->returnChromStrSubstr(tmpSNP_chrNameInt, tmpSNP_chrPos, 1);
			string tmpSNP_alterBase_1 = altBaseVec_1[tmpIndexInSNPvec];
			string tmpSNP_alterBase_2 = altBaseVec_2[tmpIndexInSNPvec];
			if(tmpSNP_refBase != tmpSNP_alterBase_1)
				SNP_1_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t"
					<< tmpSNP_refBase << "\t" << tmpSNP_alterBase_1 << endl;
			if(tmpSNP_refBase != tmpSNP_alterBase_2)
				SNP_2_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t"
					<< tmpSNP_refBase << "\t" << tmpSNP_alterBase_2 << endl;			
		}
		SNP_2_ofs.close();
		SNP_1_ofs.close();
	}
	log_ofs.close();
	return 0;
}