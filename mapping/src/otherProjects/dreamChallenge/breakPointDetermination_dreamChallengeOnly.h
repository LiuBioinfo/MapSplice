#ifndef BreakPointDetermination_dreamChallengeOnly_H
#define BreakPointDetermination_dreamChallengeOnly_H

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
#include "exonBoundaryHash_info.h"
using namespace std;

class BreakPointDetermination_dreamChallengeOnly_Info
{
private:
	int chrNameInt_gene1;
	int chrNameInt_gene2;
	int breakPointPos_gene1;
	int breakPointPos_gene2;
	string strand_gene1;
	string strand_gene2;
	int fusionCase;
public:
	BreakPointDetermination_dreamChallengeOnly_Info()
	{}

	bool determine_breakPoint(PeSam_Info& peSamInfo_body, PE_Read_Alignment_Info& peAlignInfo_end, 
		bool upstreamHead_or_downstreamTail_bool, int bufferSize, 
		ExonBoundaryHash_Info& exonBoundaryHashInfo, Index_Info* indexInfo)
	{
		if(!peAlignInfo_end.return_validForFusionPostAnalysis_bool())
			return false;
		bool peAlignInfo_end_forOrRev_bool = peAlignInfo_end.returnForMapOrRcmMap_SE_uniqAlign();
		if(upstreamHead_or_downstreamTail_bool) // fusion at upstream head
		{
			if(peAlignInfo_end_forOrRev_bool) // for for upstream end, case 2,5
			{
				if(peAlignInfo_end.returnLastJumpCodeType_SE_uniqAlign() == "S")
					return false;
				int leftAnchor_chrNameInt = peAlignInfo_end.returnChrNameInt_SE_uniqAlign(indexInfo);
				int leftAnchor_endPos = peAlignInfo_end.returnOnlySeAlign_endPos();
				int rightAnchor_chrNameInt = peSamInfo_body.returnChrNameInt();
				int rightAnchor_startPos = peSamInfo_body.returnStartPos_upstreamRead();
				for(int tmp = 0; tmp <= bufferSize; tmp++)
				{
					// positive buffer
					int tmp_breakPoint_left_candi = leftAnchor_endPos + tmp;
					int tmp_breakPoint_right_candi = rightAnchor_startPos + tmp;
					// case 2
					if(exonBoundaryHashInfo.confirm(false, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "+";
						fusionCase = 2;
						return true;
					}
					// case 5
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "-";
						fusionCase = 5;
						return true;
					}					

					// negative buffer
					tmp_breakPoint_left_candi = leftAnchor_endPos - tmp;
					tmp_breakPoint_right_candi = rightAnchor_startPos - tmp;					
					// case 2
					if(exonBoundaryHashInfo.confirm(false, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "+";
						fusionCase = 2;
						return true;
					}
					// case 5
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "-";
						fusionCase = 5;
						return true;
					}
				}
			}
			else // rev for upstream end, case 10,11
			{
				if(peAlignInfo_end.returnFirstJumpCodeType_SE_uniqAlign() == "S")
					return false;
				int leftAnchor_chrNameInt = peAlignInfo_end.returnChrNameInt_SE_uniqAlign(indexInfo);
				int leftAnchor_startPos = peAlignInfo_end.returnOnlySeAlign_startPos();
				int rightAnchor_chrNameInt = peSamInfo_body.returnChrNameInt();
				int rightAnchor_startPos = peSamInfo_body.returnStartPos_upstreamRead();
				for(int tmp = 0; tmp <= bufferSize; tmp++)
				{
					// positive buffer
					int tmp_breakPoint_left_candi = leftAnchor_startPos - tmp;
					int tmp_breakPoint_right_candi = rightAnchor_startPos + tmp;
					// case 10
					if(exonBoundaryHashInfo.confirm(true, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "-";
						strand_gene2 = "+";
						fusionCase = 10;
						return true;
					}
					// case 11
					if(exonBoundaryHashInfo.confirm(true, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "+";
						fusionCase = 11;
						return true;
					}

					// negative buffer
					tmp_breakPoint_left_candi = leftAnchor_startPos + tmp;
					tmp_breakPoint_right_candi = rightAnchor_startPos - tmp;
					// case 10
					if(exonBoundaryHashInfo.confirm(true, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "-";
						strand_gene2 = "+";
						fusionCase = 10;
						return true;
					}
					// case 11
					if(exonBoundaryHashInfo.confirm(true, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "+";
						fusionCase = 11;
						return true;
					}					
				}	
			}
		}
		else // fusion at downstream tail
		{
			if(peAlignInfo_end_forOrRev_bool) // rev rev downstream end, case 1,4
			{
				if(peAlignInfo_end.returnFirstJumpCodeType_SE_uniqAlign() == "S")
					return false;
				int leftAnchor_chrNameInt = peSamInfo_body.returnChrNameInt();
				int leftAnchor_endPos = peSamInfo_body.returnEndPos_downstreamRead();				
				int rightAnchor_chrNameInt = peAlignInfo_end.returnChrNameInt_SE_uniqAlign(indexInfo);
				int rightAnchor_startPos = peAlignInfo_end.returnOnlySeAlign_startPos();
				for(int tmp = 0; tmp <= bufferSize; tmp++)
				{
					// positive buffer
					int tmp_breakPoint_left_candi = leftAnchor_endPos + tmp;
					int tmp_breakPoint_right_candi = rightAnchor_startPos + tmp;
					// case 1
					if(exonBoundaryHashInfo.confirm(false, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "+";
						fusionCase = 1;
						return true;
					}
					// case 4
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "-";
						fusionCase = 4;
						return true;
					}					

					// negative buffer
					tmp_breakPoint_left_candi = leftAnchor_endPos - tmp;
					tmp_breakPoint_right_candi = rightAnchor_startPos - tmp;					
					// case 2
					if(exonBoundaryHashInfo.confirm(false, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "+";
						fusionCase = 1;
						return true;
					}
					// case 5
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi))
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "-";
						strand_gene2 = "-";
						fusionCase = 4;
						return true;
					}
				}
			}
			else // rev for downstream end, case 7,8
			{
				if(peAlignInfo_end.returnLastJumpCodeType_SE_uniqAlign() == "S")
					return false;
				int leftAnchor_chrNameInt = peSamInfo_body.returnChrNameInt();
				int leftAnchor_endPos = peSamInfo_body.returnEndPos_downstreamRead();
				int rightAnchor_chrNameInt = peAlignInfo_end.returnChrNameInt_SE_uniqAlign(indexInfo);
				int rightAnchor_endPos = peAlignInfo_end.returnOnlySeAlign_endPos();

				for(int tmp = 0; tmp <= bufferSize; tmp++)
				{
					// positive buffer
					int tmp_breakPoint_left_candi = leftAnchor_endPos + tmp;
					int tmp_breakPoint_right_candi = rightAnchor_endPos - tmp;

					// case 7
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(false, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "+";
						strand_gene2 = "-";
						fusionCase = 7;
						return true;
					}
					// case 8
					if(exonBoundaryHashInfo.confirm(true, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "-";
						fusionCase = 8;
						return true;
					}

					// negative buffer				
					tmp_breakPoint_left_candi = leftAnchor_endPos - tmp;
					tmp_breakPoint_right_candi = rightAnchor_endPos + tmp;
			
					// case 7
					if(exonBoundaryHashInfo.confirm(false, false, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(false, true, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = rightAnchor_chrNameInt;
						chrNameInt_gene2 = leftAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_right_candi;
						breakPointPos_gene2 = tmp_breakPoint_left_candi;
						strand_gene1 = "+";
						strand_gene2 = "-";
						fusionCase = 7;
						return true;
					}
					// case 8
					if(exonBoundaryHashInfo.confirm(true, true, leftAnchor_chrNameInt, tmp_breakPoint_left_candi)
						&&exonBoundaryHashInfo.confirm(true, false, rightAnchor_chrNameInt, tmp_breakPoint_right_candi)) 
					{
						chrNameInt_gene1 = leftAnchor_chrNameInt;
						chrNameInt_gene2 = rightAnchor_chrNameInt;
						breakPointPos_gene1 = tmp_breakPoint_left_candi;
						breakPointPos_gene2 = tmp_breakPoint_right_candi;
						strand_gene1 = "+";
						strand_gene2 = "-";
						fusionCase = 8;
						return true;
					}
				}
			}
		}
		return false;
	}

	string returnBreakPointInfoStr(Index_Info* indexInfo)
	{
		string chrNameStr_gene1 = indexInfo->returnChrNameStr(chrNameInt_gene1);
		string chrNameStr_gene2 = indexInfo->returnChrNameStr(chrNameInt_gene2);
		string tmpBreakPointInfoStr = chrNameStr_gene1 + "\t" + int_to_str(breakPointPos_gene1)
			+ "\t" + chrNameStr_gene2 + "\t" + int_to_str(breakPointPos_gene2) + "\t" 
			+ strand_gene1 + "\t" + strand_gene2 + "\tCASE_" + int_to_str(fusionCase);  
		return tmpBreakPointInfoStr;
	}

	int return_chrNameInt_gene1()
	{
		return chrNameInt_gene1;
	}

	int return_chrNameInt_gene2()
	{
		return chrNameInt_gene2;
	}

	int return_breakPointPos_gene1()
	{
		return breakPointPos_gene1;
	}

	int return_breakPointPos_gene2()
	{
		return breakPointPos_gene2;
	}

	string return_strand_gene1()
	{
		return strand_gene1;
	}

	string return_strand_gene2()
	{
		return strand_gene2;
	}

	int return_fusion_case()
	{
		return fusionCase;
	}
};
#endif