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
#include <sstream>

using namespace std;

void parsePsSJexpLine(string& tmpStr, string& tmpChrName, int& tmpStartPos, 
	int& tmpEndPos, int& tmpRefCanSJhapNum, vector< pair<int,int> >& tmpSampleHapNumExpPairVec)
{
	vector<string> fieldVec;
	int startLoc = 0;
	for(int tmp = 0; ; tmp++)
	{
		int tabLoc = tmpStr.find("\t", startLoc);
		if(tabLoc == string::npos)
		{
			fieldVec.push_back(tmpStr.substr(startLoc));
			break;
		}
		string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
		fieldVec.push_back(tmpField);
		startLoc = tabLoc + 1;
	}
	tmpChrName = fieldVec[0];
	string tmpStartPosStr = fieldVec[1];
	string tmpEndPosStr = fieldVec[2];
	string tmpRefCanSJhapNumStr = fieldVec[3];
	tmpStartPos = atoi(tmpStartPosStr.c_str());
	tmpEndPos = atoi(tmpEndPosStr.c_str());
	tmpRefCanSJhapNum = atoi(tmpRefCanSJhapNumStr.c_str());
	int tmpSampleNum = (fieldVec.size() - 4)/2;
	for(int tmp = 0; tmp < tmpSampleNum; tmp++)
	{
		string tmpSample_hapNumStr = fieldVec[tmp*2 + 4];
		string tmpSample_supNumStr = fieldVec[tmp*2 + 5];
		int tmpSample_hapNum = atoi(tmpSample_hapNumStr.c_str());
		int tmpSample_supNum = atoi(tmpSample_supNumStr.c_str());
		tmpSampleHapNumExpPairVec.push_back(pair<int,int>(tmpSample_hapNum, tmpSample_supNum));
	}
}

void get_sampleNum_supNum_varyCanSJhapNum(
	vector< pair<int,int> >& sampleCanSJhapNumExpPairVec, 
	int& canSJhapNum_0_sampleNum, int& canSJhapNum_0_supNum,
	int& canSJhapNum_1_sampleNum, int& canSJhapNum_1_supNum, 
	int& canSJhapNum_2_sampleNum, int& canSJhapNum_2_supNum)
{
	canSJhapNum_0_sampleNum = 0; 
	canSJhapNum_0_supNum = 0;
	canSJhapNum_1_sampleNum = 0;
	canSJhapNum_1_supNum = 0;
	canSJhapNum_2_sampleNum = 0;
	canSJhapNum_2_supNum = 0;
	for(int tmp = 0; tmp < sampleCanSJhapNumExpPairVec.size(); tmp++)
	{
		int tmpCanSJhapNum = sampleCanSJhapNumExpPairVec[tmp].first;
		int tmpSupNum = sampleCanSJhapNumExpPairVec[tmp].second;
		if(tmpCanSJhapNum == 0)
		{
			canSJhapNum_0_sampleNum ++;
			canSJhapNum_0_supNum += tmpSupNum;
		}
		else if(tmpCanSJhapNum == 1)
		{
			canSJhapNum_1_sampleNum ++;
			canSJhapNum_1_supNum += tmpSupNum;
		}
		else if(tmpCanSJhapNum == 2)
		{
			canSJhapNum_2_sampleNum ++;
			canSJhapNum_2_supNum += tmpSupNum;
		}
		else
		{
			cout << "invalid tmpCanSJhapNum: " << tmpCanSJhapNum << endl;
			exit(1);
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputPsSJexpression" << endl;
		cout << "#2 outputDir" << endl;
		cout << "#3 minSupNum" << endl;
		exit(1);
	}
	string inputPsSJexpression = argv[1];
	string outputDir = argv[2];
	string minSupNumStr = argv[3];
	int minSupNum = atoi(minSupNumStr.c_str());
	cout << "minSupNum: " << minSupNum << endl;
	outputDir += "/";
	string mkdir_cmd = "mkdir " + outputDir;
	system(mkdir_cmd.c_str());	

	vector<string> chrNameVec;
	vector< pair<int,int> > posPairVec;
	vector<int> refCanSJhapNumVec;
	vector< vector< pair<int,int> > > sampleCanSJhapNumExpPairVecVec;
	
	ifstream psSJexp_ifs(inputPsSJexpression.c_str());
	string tmpHeader;
	getline(psSJexp_ifs, tmpHeader);
	while(!psSJexp_ifs.eof())
	{
		string tmpStr;
		getline(psSJexp_ifs, tmpStr);
		if(tmpStr == "")
			break;
		string tmpChrName; 
		int tmpStartPos, tmpEndPos, tmpRefCanSJhapNum;
		vector< pair<int,int> > tmpSampleHapNumExpPairVec;
		parsePsSJexpLine(tmpStr, tmpChrName, tmpStartPos, tmpEndPos, tmpRefCanSJhapNum, tmpSampleHapNumExpPairVec);
		chrNameVec.push_back(tmpChrName);
		posPairVec.push_back(pair<int,int>(tmpStartPos, tmpEndPos));
		refCanSJhapNumVec.push_back(tmpRefCanSJhapNum);
		sampleCanSJhapNumExpPairVecVec.push_back(tmpSampleHapNumExpPairVec);
	}
	psSJexp_ifs.close();

	string outputFile_master = outputDir + "/master.txt";
	string outputFile_non2can = outputDir + "/non2can_master.txt";
	string outputFile_non2can_het_hom_merge = outputDir + "/non2can_het_hom_merge.txt";
	string outputFile_non2can_het = outputDir + "/non2can_het.txt";
	string outputFile_non2can_hom = outputDir + "/non2can_hom.txt";
	string outputFile_can2non = outputDir + "/can2non_master.txt";
	string outputFile_can2non_het_hom_merge = outputDir + "/can2non_het_hom_merge.txt";
	string outputFile_can2non_het = outputDir + "/can2non_het.txt";
	string outputFile_can2non_hom = outputDir + "/can2non_hom.txt";

	ofstream master_ofs(outputFile_master.c_str());

	ofstream non2can_ofs(outputFile_non2can.c_str());
	ofstream non2can_ofs_het(outputFile_non2can_het.c_str());
	ofstream non2can_ofs_hom(outputFile_non2can_hom.c_str());
	ofstream non2can_ofs_het_hom_merge(outputFile_non2can_het_hom_merge.c_str());	

	ofstream can2non_ofs(outputFile_can2non.c_str());
	ofstream can2non_ofs_het(outputFile_can2non_het.c_str());
	ofstream can2non_ofs_hom(outputFile_can2non_hom.c_str());
	ofstream can2non_ofs_het_hom_merge(outputFile_can2non_het_hom_merge.c_str());

	// master_ofs << "chrName\tpos_1\tpos_2\trefCanSJlabel";
	// master_ofs << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0";
	// master_ofs << "\tsampleNum_canSJlabel_1\tsupNum_canSJlabel_1\tavgSup_canSJlabel_1";
	// master_ofs << "\tsampleNum_canSJlabel_2\tsupNum_canSJlabel_2\tavgSup_canSJlabel_2";
	// master_ofs << "\tsampleNum_canSJlabel_non0\tsupNum_canSJlabel_non0\tavgSup_canSJlabel_non0" << endl;

	// non2can_ofs << "chrName\tpos_1\tpos_2";
	// non2can_ofs << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0";
	// non2can_ofs << "\tsampleNum_canSJlabel_1\tsupNum_canSJlabel_1\tavgSup_canSJlabel_1";
	// non2can_ofs << "\tsampleNum_canSJlabel_2\tsupNum_canSJlabel_2\tavgSup_canSJlabel_2";
	// non2can_ofs << "\tsampleNum_canSJlabel_1/2\tsupNum_canSJlabel_1/2\tavgSup_canSJlabel_1/2" << endl;
	// non2can_ofs_het << "chrName\tpos_1\tpos_2";
	// non2can_ofs_het << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0";
	// non2can_ofs_het << "\tsampleNum_canSJlabel_1\tsupNum_canSJlabel_1\tavgSup_canSJlabel_1" << endl;
	// non2can_ofs_hom << "chrName\tpos_1\tpos_2";
	// non2can_ofs_hom << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0";
	// non2can_ofs_hom << "\tsampleNum_canSJlabel_2\tsupNum_canSJlabel_2\tavgSup_canSJlabel_2" << endl;
	// non2can_ofs_het_hom_merge << "chrName\tpos_1\tpos_2";
	// non2can_ofs_het_hom_merge << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0";
	// non2can_ofs_het_hom_merge << "\tsampleNum_canSJlabel_1/2\tsupNum_canSJlabel_1/2\tavgSup_canSJlabel_1/2\tCategory" << endl;


	// can2non_ofs << "chrName\tpos_1\tpos_2";
	// can2non_ofs << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0"; 
	// can2non_ofs << "\tsampleNum_canSJlabel_1\tsupNum_canSJlabel_1\tavgSup_canSJlabel_1";
	// can2non_ofs << "\tsampleNum_canSJlabel_2\tsupNum_canSJlabel_2\tavgSup_canSJlabel_2";
	// can2non_ofs << "\tsampleNum_canSJlabel_1/2\tsupNum_canSJlabel_1/2\tavgSup_canSJlabel_1/2" << endl;
	// can2non_ofs_het << "chrName\tpos_1\tpos_2";
	// can2non_ofs_het << "\tsampleNum_canSJlabel_2\tsupNum_canSJlabel_2\tavgSup_canSJlabel_2";
	// can2non_ofs_het << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0" << endl;
	// can2non_ofs_hom << "chrName\tpos_1\tpos_2";
	// can2non_ofs_hom << "\tsampleNum_canSJlabel_1\tsupNum_canSJlabel_1\tavgSup_canSJlabel_1";
	// can2non_ofs_hom << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0" << endl;
	// can2non_ofs_het_hom_merge << "chrName\tpos_1\tpos_2";
	// can2non_ofs_het_hom_merge << "\tsampleNum_canSJlabel_2/1\tsupNum_canSJlabel_2/1\tavgSup_canSJlabel_2/1";
	// can2non_ofs_het_hom_merge << "\tsampleNum_canSJlabel_0\tsupNum_canSJlabel_0\tavgSup_canSJlabel_0\tCategory" << endl;


	cout << "refCanSJhapNumVec.size(): " << refCanSJhapNumVec.size() << endl;
	for(int tmp = 0; tmp < refCanSJhapNumVec.size(); tmp++)
	{
		int canSJhapNum_0_sampleNum = 0, canSJhapNum_1_sampleNum = 0, canSJhapNum_2_sampleNum = 0,
			canSJhapNum_0_supNum = 0, canSJhapNum_1_supNum = 0, canSJhapNum_2_supNum = 0;
		get_sampleNum_supNum_varyCanSJhapNum(
			sampleCanSJhapNumExpPairVecVec[tmp], 
			canSJhapNum_0_sampleNum, canSJhapNum_0_supNum,
			canSJhapNum_1_sampleNum, canSJhapNum_1_supNum, 
			canSJhapNum_2_sampleNum, canSJhapNum_2_supNum);
		int canSJhapNum_non0_sampleNum = canSJhapNum_1_sampleNum + canSJhapNum_2_sampleNum;
		int canSJhapNum_non0_supNum = canSJhapNum_1_supNum + canSJhapNum_2_supNum;

		double tmpSJcov_0_avg, tmpSJcov_1_avg, tmpSJcov_2_avg, tmpSJcov_non0_avg;// tmpSJcov_non2_avg;
		if(canSJhapNum_0_sampleNum == 0)
			tmpSJcov_0_avg = 0;
		else
			tmpSJcov_0_avg = (double)canSJhapNum_0_supNum/(double)canSJhapNum_0_sampleNum;
		if(canSJhapNum_1_sampleNum == 0)
			tmpSJcov_1_avg = 0;
		else
			tmpSJcov_1_avg = (double)canSJhapNum_1_supNum/(double)canSJhapNum_1_sampleNum;
		if(canSJhapNum_2_sampleNum == 0)
			tmpSJcov_2_avg = 0;
		else
			tmpSJcov_2_avg = (double)canSJhapNum_2_supNum/(double)canSJhapNum_2_sampleNum;	
		if(canSJhapNum_non0_sampleNum == 0)
			tmpSJcov_non0_avg = 0;
		else
			tmpSJcov_non0_avg = (double)canSJhapNum_non0_supNum/(double)canSJhapNum_non0_sampleNum;

		if((refCanSJhapNumVec[tmp] == 0)&&(tmpSJcov_non0_avg >= minSupNum)) // non2can
		{	
			if((canSJhapNum_0_sampleNum > 0)&&(canSJhapNum_non0_sampleNum > 0))
				non2can_ofs << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg
					<< "\t" << canSJhapNum_non0_sampleNum << "\t" << canSJhapNum_non0_supNum << "\t" << tmpSJcov_non0_avg << endl;
			if((canSJhapNum_0_sampleNum > 0)&&(canSJhapNum_1_sampleNum > 0))
				non2can_ofs_het << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg 
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg << endl;
			if((canSJhapNum_0_sampleNum > 0)&&(canSJhapNum_2_sampleNum > 0))
				non2can_ofs_hom << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg 
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg << endl;
			if((canSJhapNum_0_sampleNum > 0)&&(canSJhapNum_1_sampleNum > 0)&&(canSJhapNum_2_sampleNum > 0))
			{
				non2can_ofs_het_hom_merge << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg 
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg << "\tHet" << endl;
				non2can_ofs_het_hom_merge << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg 
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg << "\tHom" << endl;									
			}
		}
		else if((refCanSJhapNumVec[tmp] == 1)&&(tmpSJcov_non0_avg >= minSupNum)) // can2non
		{
			if((canSJhapNum_non0_sampleNum > 0)&&(canSJhapNum_0_sampleNum > 0))
				can2non_ofs << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second  
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg 
					<< "\t" << canSJhapNum_non0_sampleNum << "\t" << canSJhapNum_non0_supNum << "\t" << tmpSJcov_non0_avg << endl;
			if((canSJhapNum_1_sampleNum > 0)&&(canSJhapNum_0_sampleNum > 0))
				can2non_ofs_het << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg << endl;
			if((canSJhapNum_2_sampleNum > 0)&&(canSJhapNum_0_sampleNum > 0))
				can2non_ofs_hom << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg << endl;
			if((canSJhapNum_2_sampleNum > 0)&&(canSJhapNum_1_sampleNum > 0)&&(canSJhapNum_0_sampleNum > 0))
			{
				can2non_ofs_het_hom_merge << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg << "\tHet" << endl;
				can2non_ofs_het_hom_merge << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" << posPairVec[tmp].second 
					<< "\t" << canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg 
					<< "\t" << canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg << "\tHom" << endl;
			}
		}
		else if((refCanSJhapNumVec[tmp] == 0)||(refCanSJhapNumVec[tmp] == 1))
		{}
		else
		{
			cout << "invalid refCanSJhapNumVec[tmp]: " << refCanSJhapNumVec[tmp] << endl;
			exit(1);
		}

		master_ofs << chrNameVec[tmp] << "\t" << posPairVec[tmp].first << "\t" 
			<< posPairVec[tmp].second << "\t" << refCanSJhapNumVec[tmp] << "\t"
			<< canSJhapNum_0_sampleNum << "\t" << canSJhapNum_0_supNum << "\t" << tmpSJcov_0_avg << "\t"
			<< canSJhapNum_1_sampleNum << "\t" << canSJhapNum_1_supNum << "\t" << tmpSJcov_1_avg << "\t"
			<< canSJhapNum_2_sampleNum << "\t" << canSJhapNum_2_supNum << "\t" << tmpSJcov_2_avg << "\t"
			<< canSJhapNum_non0_sampleNum << "\t" << canSJhapNum_non0_supNum << "\t" << tmpSJcov_non0_avg << endl;			
	}
	master_ofs.close();
	non2can_ofs.close();
	non2can_ofs_het.close();
	non2can_ofs_hom.close();
	can2non_ofs.close();
	can2non_ofs_het.close();
	can2non_ofs_hom.close();
	return 0;
}