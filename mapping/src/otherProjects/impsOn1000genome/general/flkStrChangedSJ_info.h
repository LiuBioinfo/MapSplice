#ifndef FLKSTRCHANGEDSJ_INFO_H
#define FLKSTRCHANGEDSJ_INFO_H

using namespace std;

typedef map<int, int> AcceptorStartPos2indexMap;
typedef map<int, map<int,int> > DonerEndPos2mapMap;

class FlkStrChangedSJ_Info
{
private:
	int chrNameInt;
	int donerEndPos;
	int acceptorStartPos;
	int supNum;
	string flkStr_ref;
	string flkStr_hap1;
	string flkStr_hap2;
public:
	FlkStrChangedSJ_Info()
	{}

	void initiate(int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos,
		int tmpSupNum, string& tmpFlkStr_ref, string& tmpFlkStr_hap1, string& tmpFlkStr_hap2)
	{
		chrNameInt = tmpChrNameInt;
		donerEndPos = tmpDonerEndPos;
		acceptorStartPos = tmpAcceptorStartPos;
		supNum = tmpSupNum;
		flkStr_ref = tmpFlkStr_ref;
		flkStr_hap1 = tmpFlkStr_hap1;
		flkStr_hap2 = tmpFlkStr_hap2;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnDonerEndPos()
	{
		return donerEndPos;
	}

	int returnAcceptorStartPos()
	{
		return acceptorStartPos;
	}

	int returnSupNum()
	{
		return supNum;
	}

	string returnFlkStr_ref()
	{
		return flkStr_ref;
	}

	string returnFlkStr_hap1()
	{
		return flkStr_hap1;
	}

	string returnFlkStr_hap2()
	{
		return flkStr_hap2;
	}	

	void addSupNum(int tmpToAddSupNum)
	{
		supNum += tmpToAddSupNum;
	}
};

class FlkStrChangedSJ_Hash_Info
{
private:
	vector<DonerEndPos2mapMap> donerEndPos2mapMapVec;
	vector<FlkStrChangedSJ_Info> flkStrChangedSJinfoVec;
public:
	FlkStrChangedSJ_Hash_Info()
	{}

	bool isCanSJ(string& tmpFlkStr)
	{
		if((tmpFlkStr == "GTAG")||(tmpFlkStr == "CTAC")
			||(tmpFlkStr == "ATAC")||(tmpFlkStr == "GTAT")
			||(tmpFlkStr == "GCAG")||(tmpFlkStr == "CTGC"))
			return true;
		else
			return false;
	}

	void search_and_return_canSJhapNum_supNum(int& output_canSJhapNum, int& output_supNum,
		int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int tmpSJvecIndex;
		bool tmpSearchBool = this->searchForSJ(tmpSJvecIndex, tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(!tmpSearchBool)
		{
			output_canSJhapNum = 0;
			output_supNum = 0;
		}
		else
		{
			output_canSJhapNum = 0;
			string tmpFlkStr_hap1 = this->return_flkStr_hap1_vecIndex(tmpSJvecIndex);
			string tmpFlkStr_hap2 = this->return_flkStr_hap2_vecIndex(tmpSJvecIndex);
			if(this->isCanSJ(tmpFlkStr_hap1))
				output_canSJhapNum ++;
			if(this->isCanSJ(tmpFlkStr_hap2))
				output_canSJhapNum ++;
			output_supNum = return_supNum_vecIndex(tmpSJvecIndex);
		}
	}

	int search_and_return_supNum(int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		int tmpSJvecIndex;
		bool tmpSearchBool = this->searchForSJ(tmpSJvecIndex, tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos);
		if(!tmpSearchBool)
			return 0;
		else
			return (this->return_supNum_vecIndex(tmpSJvecIndex));
	}	

	void initiate_chrNum(int totalChrNum)
	{
		for(int tmp = 0; tmp < totalChrNum; tmp++)
		{
			DonerEndPos2mapMap tmpNewDonerEndPos2mapMap;
			donerEndPos2mapMapVec.push_back(tmpNewDonerEndPos2mapMap);
		}
	}

	int return_flkStrChangedSJnum()
	{
		return flkStrChangedSJinfoVec.size();
	}

	int return_chrNameInt_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnChrNameInt();
	}

	int return_donerEndPos_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnDonerEndPos();
	}

	int return_acceptorStartPos_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnAcceptorStartPos();
	}

	int return_supNum_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnSupNum();
	}

	string return_flkStr_ref_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnFlkStr_ref();
	}

	string return_flkStr_hap1_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnFlkStr_hap1();
	}

	string return_flkStr_hap2_vecIndex(int tmpVecIndex)
	{
		return flkStrChangedSJinfoVec[tmpVecIndex].returnFlkStr_hap2();
	}		

	void insertFlkStrChangedSJ(int& tmpChrNameInt, int& tmpDonerEndPos, int& tmpAcceptorStartPos, 
		int& tmpSupNum, string& tmpFlkStr_ref, string& tmpFlkStr_hap1, string& tmpFlkStr_hap2)
	{
		FlkStrChangedSJ_Info tmpFlkStrChangedSJinfo;
		tmpFlkStrChangedSJinfo.initiate(tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos, 
			tmpSupNum, tmpFlkStr_ref, tmpFlkStr_hap1, tmpFlkStr_hap2);
		flkStrChangedSJinfoVec.push_back(tmpFlkStrChangedSJinfo);
		int tmpVecIndex = flkStrChangedSJinfoVec.size() - 1;
		DonerEndPos2mapMap::iterator tmpIter_1 = donerEndPos2mapMapVec[tmpChrNameInt].find(tmpDonerEndPos);
		if(tmpIter_1 != donerEndPos2mapMapVec[tmpChrNameInt].end()) // donerEndPos exists;
			(tmpIter_1->second).insert(pair<int,int>(tmpAcceptorStartPos, tmpVecIndex));
		else
		{
			AcceptorStartPos2indexMap tmpMap;
			tmpMap.insert(pair<int,int>(tmpAcceptorStartPos, tmpVecIndex));
			donerEndPos2mapMapVec[tmpChrNameInt].insert(pair<int, AcceptorStartPos2indexMap>(tmpDonerEndPos, tmpMap));
		}
	}

	bool searchForSJ(int& tmpSJvecIndex, int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos)
	{
		DonerEndPos2mapMap::iterator tmpIter_1 = donerEndPos2mapMapVec[tmpChrNameInt].find(tmpDonerEndPos);
		if(tmpIter_1 != donerEndPos2mapMapVec[tmpChrNameInt].end())
		{
			AcceptorStartPos2indexMap::iterator tmpIter_2 = (tmpIter_1->second).find(tmpAcceptorStartPos);
			if(tmpIter_2 != (tmpIter_1->second).end())
			{
				tmpSJvecIndex = tmpIter_2->second;
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}

	void loadFlkStrChangedSJfile(string& flkStrChangedSJFile, Index_Info* indexInfo)
	{
		ifstream SJ_ifs(flkStrChangedSJFile.c_str());
		while(!SJ_ifs.eof())
		{
			string tmpStr;
			getline(SJ_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos, tmpSupNum;
			string tmpFlkStr_ref, tmpFlkStr_hap1, tmpFlkStr_hap2;
			bool parse_success_bool = parseSJline(tmpStr, tmpChrNameInt, tmpDonerEndPos, 
				tmpAcceptorStartPos, tmpSupNum, tmpFlkStr_ref, tmpFlkStr_hap1, tmpFlkStr_hap2, indexInfo);
			if(parse_success_bool)
			{
				int tmpSJvecIndex;
				bool searchSJinMapBool = this->searchForSJ(tmpSJvecIndex, tmpChrNameInt, 
					tmpDonerEndPos, tmpAcceptorStartPos);
				if(searchSJinMapBool)
					flkStrChangedSJinfoVec[tmpSJvecIndex].addSupNum(tmpSupNum);
				else
					this->insertFlkStrChangedSJ(tmpChrNameInt, tmpDonerEndPos, tmpAcceptorStartPos, 
						tmpSupNum, tmpFlkStr_ref, tmpFlkStr_hap1, tmpFlkStr_hap2);
			}
		}
		SJ_ifs.close();
	}

	bool parseSJline(string& tmpStr, int& tmpChrNameInt, int& tmpDonerEndPos, int& tmpAcceptorStartPos, 
		int& tmpSupNum, string& tmpFlkStr_ref, string& tmpFlkStr_hap1, string& tmpFlkStr_hap2, Index_Info* indexInfo)
	{
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		int tabLoc_4 = tmpStr.find("\t", tabLoc_3 + 1);
		int tabLoc_5 = tmpStr.find("\t", tabLoc_4 + 1);
		int tabLoc_6 = tmpStr.find("\t", tabLoc_5 + 1);
		int tabLoc_7 = tmpStr.find("\t", tabLoc_6 + 1);
		string tmpChrNameStr = tmpStr.substr(0, tabLoc_1);
		string tmpDonerEndPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpAcceptorStartPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpSupNumStr = tmpStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		tmpFlkStr_ref = tmpStr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
		tmpFlkStr_hap1 = tmpStr.substr(tabLoc_5 + 1, tabLoc_6 - tabLoc_5 - 1);
		tmpFlkStr_hap2 = tmpStr.substr(tabLoc_6 + 1, tabLoc_7 - tabLoc_6 - 1);
		tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		tmpDonerEndPos = atoi(tmpDonerEndPosStr.c_str());
		tmpAcceptorStartPos = atoi(tmpAcceptorStartPosStr.c_str());
		tmpSupNum = atoi(tmpSupNumStr.c_str());
		if(tmpChrNameInt < 0)
			return false;
		else
			return true;
	}
};
#endif