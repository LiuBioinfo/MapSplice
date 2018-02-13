// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALU_TOTAL_INFO_H
#define ALU_TOTAL_INFO_H

#include "../../../../general/index_info.h"
#include "alu_element_info.h"

using namespace std;

typedef map<string, int> Id2elementIndexMap;

class Alu_Total_Info
{
private:
	Id2elementIndexMap aluId2vecIndexMap;
	vector<Alu_Element_Info> aluElementInfoVec;
public:
	Alu_Total_Info()
	{}

	void initiate_wholeGenomePos2aluelementVecIndexVec(vector<int>& wholeGenomePos2aluelementVecIndexVec, Index_Info* indexInfo)
	{
		int aluElementInfoVecSize = aluElementInfoVec.size();
		for(int tmp = 0; tmp < aluElementInfoVecSize; tmp++)
		{
			int tmpAluElement_chrNameInt = aluElementInfoVec[tmp].return_chrNameInt();
			int tmpAluElement_startPos_inChr = aluElementInfoVec[tmp].return_startPos();
			int tmpAluElement_endPos_inChr = aluElementInfoVec[tmp].return_endPos();
			unsigned int tmpAluElement_startPos_inWholeGenome = indexInfo->getWholeGenomeLocation(
				tmpAluElement_chrNameInt, tmpAluElement_startPos_inChr);
			unsigned int tmpAluElement_endPos_inWholeGenome = indexInfo->getWholeGenomeLocation(
				tmpAluElement_chrNameInt, tmpAluElement_endPos_inChr);
			wholeGenomePos2aluelementVecIndexVec[tmpAluElement_startPos_inWholeGenome] = tmp;
			wholeGenomePos2aluelementVecIndexVec[tmpAluElement_endPos_inWholeGenome] = tmp;
		}
	}

	void parseAluInfoStr(string& tmpAluStr, string& tmpAluFamilyName, string& tmpChrName, int& tmpStartPos, int& tmpEndPos)
	{
		vector<string> tmpFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpAluStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{	
				tmpFieldStrVec.push_back(tmpAluStr.substr(startLoc));
				break;
			}
			tmpFieldStrVec.push_back(tmpAluStr.substr(startLoc, tabLoc - startLoc));
			startLoc = tabLoc + 1;
		}

		tmpAluFamilyName = tmpFieldStrVec[2];
		tmpChrName = tmpFieldStrVec[3];
		tmpStartPos = atoi((tmpFieldStrVec[0]).c_str());
		tmpEndPos = atoi((tmpFieldStrVec[1]).c_str());
	}

	void initiate_withAluInfoFile(string& aluInfoFile, Index_Info* indexInfo, string& invalidAluInfoFile, string& validAluInfoFile)
	{
		ofstream invalidAluInfo_ofs(invalidAluInfoFile.c_str());
		ofstream validAluInfo_ofs(validAluInfoFile.c_str());
		ifstream aluInfo_ifs(aluInfoFile.c_str());
		int tmpAluElementIndex = 0;
		while(!aluInfo_ifs.eof())
		{
			string tmpStr;
			getline(aluInfo_ifs, tmpStr);
			if(tmpStr == "")
				break;
			string tmpAluFamilyName, tmpChrName;
			int tmpStartPos, tmpEndPos;
			parseAluInfoStr(tmpStr, tmpAluFamilyName, tmpChrName, tmpStartPos, tmpEndPos);
			int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
			if(tmpChrNameInt < 0)
			{
				cout << "invalidChrName: " << tmpChrName << endl;
				cout << "invalidChrName.size(): " << tmpChrName.length() << endl;
				invalidAluInfo_ofs << tmpStr << endl;
				continue;
			}
			validAluInfo_ofs << tmpChrName << "\t" << tmpStartPos << "\t" << tmpEndPos << "\t" << tmpAluFamilyName << endl;
			Alu_Element_Info tmpAluElementInfo;
			string tmpAluElement_id = tmpAluFamilyName + "_" + tmpChrName + "_" + int_to_str(tmpStartPos) + "_" + int_to_str(tmpEndPos) + "_";
			tmpAluElementInfo.initiate(tmpAluElement_id, tmpAluFamilyName, tmpChrNameInt, tmpStartPos, tmpEndPos);
			aluId2vecIndexMap.insert(pair<string,int>(tmpAluElement_id, tmpAluElementIndex));
			aluElementInfoVec.push_back(tmpAluElementInfo);
			tmpAluElementIndex ++;
		}
		aluInfo_ifs.close();
		validAluInfo_ofs.close();
		invalidAluInfo_ofs.close();
	}

	void estimateAluEleAbundance(string& inputSamFile, string& outputInvalidSamFile, string& readAssignmentFile, 
		string& outputAluEleEstimatedAbundanceReadCountFile, Index_Info* indexInfo)
	{
		this->convertSam2readAssignmentFile(inputSamFile, outputInvalidSamFile, readAssignmentFile);
		this->EM_estimateAbundance(readAssignmentFile);
		this->output_aluEleAbundanceReadCount(outputAluEleEstimatedAbundanceReadCountFile, indexInfo);
	}

	void convertSam2readAssignmentFile(string& inputSamFile, string& outputInvalidSamFile, string& readAssignmentFile)
	{

	}

	void EM_estimateAbundance(string& readAssignmentFile) 
	// each line in readAssignmentFile: read_id  mapPosNum(N)  mappedAluEleIndex_1,mappedAluEleIndex_2,...mappedAluEleIndex_N,  mapQual_1,mapQual_2,...mapQual_N
	{
		this->initiate_EM(readAssignmentFile);
		this->update_EM(readAssignmentFile);
	}

	void parseReadAssignmentStr(string& tmpAluStr, string& tmpReadIdStr, int& tmpMapPosNum, vector<int>& tmpMappedAluEleIndexVec, vector<double>& tmpMapQualVec)
	// each line in readAssignmentFile: read_id  mapPosNum(N)  mappedAluEleIndex_1,mappedAluEleIndex_2,...mappedAluEleIndex_N,  mapQual_1,mapQual_2,...mapQual_N
	{
		vector<string> tmpFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpAluStr.find("\t", startLoc);
			if(tabLoc == string::npos)
			{	
				tmpFieldStrVec.push_back(tmpAluStr.substr(startLoc));
				break;
			}
			tmpFieldStrVec.push_back(tmpAluStr.substr(startLoc, tabLoc - startLoc));
			startLoc = tabLoc + 1;
		}
		tmpReadIdStr = tmpFieldStrVec[0];
		tmpMapPosNum = atoi(tmpFieldStrVec[1].c_str());
		
		// generate tmpMappedAluEleIndexVec
		string tmpMappedAluEleIndexVecStr = tmpFieldStrVec[2];
		startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpMappedAluEleIndexVecStr.find(",", startLoc);
			if(tabLoc == string::npos)
			{	
				string tmpIndexStr = tmpMappedAluEleIndexVecStr.substr(startLoc);
				int tmpIndex = atoi(tmpIndexStr.c_str());
				tmpMappedAluEleIndexVec.push_back(tmpIndex);
				break;
			}
			string tmpIndexStr = tmpMappedAluEleIndexVecStr.substr(startLoc, tabLoc - startLoc);
			int tmpIndex = atoi(tmpIndexStr.c_str());
			tmpMappedAluEleIndexVec.push_back(tmpIndex);
			startLoc = tabLoc + 1;
		}

		// generate tmpMapQualVec
		string tmpMapQualVecStr = tmpFieldStrVec[3];
		startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpMapQualVecStr.find(",", startLoc);
			if(tabLoc == string::npos)
			{	
				string tmpMapQualStr = tmpMapQualVecStr.substr(startLoc);
				int tmpMapQual = atoi(tmpMapQualStr.c_str());
				tmpMapQualVec.push_back(tmpMapQual);
				break;
			}
			string tmpMapQualStr = tmpMapQualVecStr.substr(startLoc, tabLoc - startLoc);
			int tmpMapQual = atoi(tmpMapQualStr.c_str());
			tmpMapQualVec.push_back(tmpMapQual);
			startLoc = tabLoc + 1;
		}		
	}

	void initiateReadCount2aluElementVec(int tmpMapPosNum, vector<int>& tmpMappedAluEleIndexVec, vector<double>& tmpMapQualVec)
	// average geneCount on multi candidate positions
	{
		if(tmpMapPosNum == 1)
		{
			int toAddReadCountAluElementIndex = tmpMappedAluEleIndexVec[0];
			double toAddReadCount = 1.0;
			this->addReadCount2correspondingAluElement(toAddReadCountAluElementIndex, toAddReadCount);
		}
		else // multiple mapping positions
		{
			double tmp_mapQual_sum = sumUpMapQualVec(tmpMapQualVec);
			for(int tmp = 0; tmp < tmpMapPosNum; tmp++)
			{
				int tmp_toAddReadCountAluElementIndex = tmpMappedAluEleIndexVec[tmp];
				double tmp_toAddReadCount = tmpMapQualVec[tmp] / tmp_mapQual_sum;
				this->addReadCount2correspondingAluElement(tmp_toAddReadCountAluElementIndex, tmp_toAddReadCount);
			}
		}
	}

	void addReadCount2correspondingAluElement(int tmp_toAddIndex, double tmp_toAddReadCount)
	{
		aluElementInfoVec[tmp_toAddIndex].addReadCount(tmp_toAddReadCount);
	}

	double sumUpMapQualVec(vector<double>& tmpMapQualVec)
	{
		double tmpSum = 0.0;
		for(int tmp = 0; tmp < tmpMapQualVec.size(); tmp++)
			tmpSum += tmpMapQualVec[tmp];
		return tmpSum;
	}

	void initiate_EM(string& readAssignmentFile)
	{
		ifstream readAssignment_ifs(readAssignmentFile.c_str());
		while(!readAssignment_ifs.eof())
		{
			string tmpStr;
			getline(readAssignment_ifs, tmpStr);
			if(tmpStr == "")
				break;
			string tmpReadIdStr;
			int tmpMapPosNum;
			vector<int> tmpMappedAluEleIndexVec;
			vector<double> tmpMapQualVec;
			parseReadAssignmentStr(tmpStr, tmpReadIdStr, tmpMapPosNum, tmpMappedAluEleIndexVec, tmpMapQualVec);
			this->initiateReadCount2aluElementVec(tmpMapPosNum, tmpMappedAluEleIndexVec, tmpMapQualVec);
		}
		readAssignment_ifs.close();	
	}

	void update_EM(string& readAssignmentFile)
	{

	}

	void output_aluEleAbundanceReadCount(string& outputAluEleEstimatedAbundanceReadCountFile, Index_Info* indexInfo)
	{
		ofstream aluEleEstimatedAbundanceReadCount_ofs(outputAluEleEstimatedAbundanceReadCountFile.c_str());
		int totalAluElementInfoVecSize = aluElementInfoVec.size();
		aluEleEstimatedAbundanceReadCount_ofs << "Id\tChrName\tStartPos\tEndPos\tLength\tReadCount\tAbundance" << endl;
		for(int tmp = 0; tmp < totalAluElementInfoVecSize; tmp++)
		{
			string tmpAluEle_id = this->return_id_withIndex(tmp);
			int tmpAluEle_chrNameInt = this->return_chrNameInt_withIndex(tmp);
			string tmpAluEle_chrName = indexInfo->returnChrNameStr(tmpAluEle_chrNameInt);
			int tmpAluEle_startPos = this->return_startPos_withIndex(tmp);
			int tmpAluEle_endPos = this->return_endPos_withIndex(tmp);
			int tmpAluEle_length = this->return_length_withIndex(tmp);
			double tmpAluEle_readCount = this->return_readCount_withIndex(tmp);
			double tmpAluEle_abundance = this->return_abundance_withIndex(tmp);
			aluEleEstimatedAbundanceReadCount_ofs << tmpAluEle_id << "\t" << tmpAluEle_chrName 
				<< "\t" << tmpAluEle_startPos << "\t" << tmpAluEle_endPos << "\t" << tmpAluEle_length 
				<< "\t" << tmpAluEle_readCount << "\t" << tmpAluEle_abundance << endl;
		}
		aluEleEstimatedAbundanceReadCount_ofs.close();
	}

	string return_id_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_id();
	}

	int return_chrNameInt_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_chrNameInt();
	}

	int return_startPos_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_startPos();
	}

	int return_endPos_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_endPos();
	}

	int return_length_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_length();
	}

	double return_readCount_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_readCount();
	}

	double return_abundance_withIndex(int tmpIndex)
	{
		return aluElementInfoVec[tmpIndex].return_abundance();
	}

};
#endif