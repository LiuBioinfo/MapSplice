// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef REISSUEDGENOMEID2TAXOID_INFO_H
#define REISSUEDGENOMEID2TAXOID_IFNO_H
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
#include "../../mps3Lib/read_block_test.h"
#include "../../mps3Lib/otherFunc.h"
#include "NCBIfullTaxoID2Name_info.h"
using namespace std;

class ReissuedGenomeID2TaxoID_Info
{
private:
	vector<int> oriGenomeIdVec;

	vector<int> speciesIdVec;
	vector<int> genusIdVec;
	vector<int> phylumIdVec;

	vector<string> speciesNameVec;
	vector<string> genusNameVec;
	vector<string> phylumNameVec;	
public:
	ReissuedGenomeID2TaxoID_Info()
	{}

	int returnTaxoId_fromReissuedId_species(int tmpReissuedId)
	{
		return speciesIdVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId_genus(int tmpReissuedId)
	{
		return genusIdVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId_phylum(int tmpReissuedId)
	{
		return phylumIdVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_species(int tmpReissuedId)
	{
		return speciesNameVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_genus(int tmpReissuedId)
	{
		return genusNameVec[tmpReissuedId - 1];
	}

	string returnTaxoName_fromReissuedId_phylum(int tmpReissuedId)
	{
		return phylumNameVec[tmpReissuedId - 1];
	}

	int returnTaxoId_fromReissuedId(int tmpReissuedId, int tmpRank)
	//1-domain, 2-kindom, 3-phylum, 4-class, 
	//5-order, 6-family, 7-genus, 8-species
	{
		if(tmpRank == 3)
			return this->returnTaxoId_fromReissuedId_phylum(tmpReissuedId);
		else if(tmpRank = 7)
			return this->returnTaxoId_fromReissuedId_genus(tmpReissuedId);
		else if(tmpRank = 8)
			return this->returnTaxoId_fromReissuedId_species(tmpReissuedId);		
		else
		{
			cout << "invalid rank: " << tmpRank << endl;
			exit(1);
		}
	}

	string returnTaxoName_fromReissuedId(int tmpReissuedId, int tmpRank)
	//1-domain, 2-kindom, 3-phylum, 4-class, 
	//5-order, 6-family, 7-genus, 8-species
	{
		if(tmpRank == 3)
			return this->returnTaxoName_fromReissuedId_phylum(tmpReissuedId);
		else if(tmpRank = 7)
			return this->returnTaxoName_fromReissuedId_genus(tmpReissuedId);
		else if(tmpRank = 8)
			return this->returnTaxoName_fromReissuedId_species(tmpReissuedId);		
		else
		{
			cout << "invalid rank: " << tmpRank << endl;
			exit(1);
		}
	}

	void initiate_reissuedId2oriGenomeIdFile_oriGenomeId2taxoIdFile_NCBIfullTaxoId2NameFile(
		string& reissuedId2oriGenomeIdFile, string& oriGenomeid2taxoIdFile)
	{
		ifstream oriGenomeid2taxoId_ifs(oriGenomeid2taxoIdFile.c_str());
		vector<int> interOriGenomeIdVec;
		vector<int> interSpeciesIdVec;
		vector<int> interGenusIdVec;
		vector<int> interPhylumIdVec;
		while(!oriGenomeid2taxoId_ifs.eof())
		{
			string tmpStr;
			getlien(oriGenomeid2taxoId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			vector<string> tmpFieldVec;
			parseStr2fieldVec(tmpFieldVec, tmpStr);
			int tmpOriGenomeId = atoi(tmpFieldVec[0].c_str());
			int tmpSpeciesId = atoi(tmpFieldVec[1].c_str());
			int tmpGenusId = atoi(tmpFieldVec[2].c_str());
			int tmpPhylumId = atoi(tmpFieldVec[6].c_str());
			interOriGenomeIdVec.push_back(tmpOriGenomeId);
			interSpeciesIdVec.push_back(tmpSpeciesId);
			interGenusIdVec.push_back(tmpGenusId);
			interPhylumIdVec.push_back(tmpPhylumId);
		}
		oriGenomeid2taxoId_ifs.close();
	
		ifstream reissuedId2oriGenomeId_ifs(reissuedId2oriGenomeIdFile.c_str());
		while(!reissuedId2oriGenomeId_ifs.eof())
		{
			string tmpStr;
			getline(reissuedId2oriGenomeId_ifs, tmpStr);
			if(tmpStr == "")
				break;
			int tabLoc = tmpStr.find("\t");
			string tmpReissuedIdStr = tmpStr.substr(0, tabLoc);
			string tmpOriGenomeIdStr = tmpStr.substr(tabLoc + 1);
			int tmpReissuedId = atoi(tmpReissuedIdStr.c_str());
			int tmpOriGenomeId = atoi(tmpOriGenomeIdStr.c_str());
		}
		reissuedId2oriGenomeId_ifs.close();

		this->initiateNameVec_idVec_NCBIfullTaxoId2NameFile(NCBIfullTaxoId2NameFile);
	}

	void initiateNameVec_idVec_NCBIfullTaxoId2NameFile(string& NCBIfullTaxoId2NameFile)
	{
		NCBIfullTaxoID2Name_Info fullTaxoId2NameInfo;
		fullTaxoId2NameInfo.initiate_taxoID2NameFile(NCBIfullTaxoId2NameFile);
		int NCBIfullTaxoIdMax = fullTaxoId2NameInfo.return_NCBIfullTaxoIdMax();
		int genomeIdVecSize = oriGenomeIdVec.size();
		for(int tmp = 0; tmp < genomeIdVecSize; tmp++)
		{
			int tmpSpeciesId = speciesIdVec[tmp];
			int tmpGenusId = genusIdVec[tmp];
			int tmpPhylumId = phylumIdVec[tmp];
			if((tmpSpeciesId < 0)||(tmpSpeciesId > NCBIfullTaxoIdMax))
				speciesNameVec.push_back("NULL");
			else
				speciesNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpSpeciesId));

			if((tmpGenusId < 0)||(tmpGenusId > NCBIfullTaxoIdMax))
				genusNameVec.push_back("NULL");
			else
				genusNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpGenusId));

			if((tmpPhylumId < 0)||(tmpPhylumId > NCBIfullTaxoIdMax))
				phylumNameVec.push_back("NULL");
			else
				phylumNameVec.push_back(fullTaxoId2NameInfo.return_taxoName(tmpPhylumId));			
		}
	}

	void parseStr2fieldVec(vector<string>& tmpFieldVec, string& tmpStr)
	{
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpField = tmpStr.substr(startLoc, tabLoc-startLoc);
			tmpFieldVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		tmpFieldVec.push_back(tmpStr.substr(startLoc));
	}
};
#endif