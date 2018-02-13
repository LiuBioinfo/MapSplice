#ifndef FIXDOUBLEANCHORSPLICE_MULTIREF_INFO_H
#define FIXDOUBLEANCHORSPLICE_MULTIREF_INFO_H

//#include "fixDoubleAnchorSplice_complicate_info.h"
#include "nw_DP.h"

class FixDoubleAnchor_Splice_MultiRef_Info
{
private:
		vector<int> donerMapCumulativeMismatchNumVec_sd;
		vector<int> donerMapMismatchPosVec_sd;
		vector<char> donerMapMismatchCharVec_sd;
		vector<int> acceptorMapCumulativeMismatchNumVec_sd;
		vector<int> acceptorMapMismatchPosVec_sd;
		vector<char> acceptorMapMismatchCharVec_sd;
		int max_extension_forward_doner_sd;
		int max_extension_backward_acceptor_sd;
		vector< pair< int, pair<int, int> > > spliceSiteVec_sd; // < < splice_site = forwared_doner_# , mismatch # >, flank_string_case >
		int bestSpliceSiteIndex_sd;
		vector<int> bestSpliceMismatchPosVec_sd;
		vector<char> bestSpliceMismatchCharVec_sd;

		vector< vector<int> > donerMapCumulativeMismatchNumVec_altVec;
		vector< vector<int> > donerMapMismatchPosVec_altVec;
		vector< vector<char> > donerMapMismatchCharVec_altVec;
		vector< vector<int> > acceptorMapCumulativeMismatchNumVec_altVec;
		vector< vector<int> > acceptorMapMismatchPosVec_altVec;
		vector< vector<char> > acceptorMapMismatchCharVec_altVec;
		vector<int> max_extension_forward_doner_altVec;
		vector<int> max_extension_backward_acceptor_altVec;		
		vector< vector< pair< int, pair<int, int> > > > spliceSiteVec_altVec;
		vector<int> bestSpliceSiteIndex_altVec;
		vector< vector<int> > bestSpliceMismatchPosVec_altVec;
		vector< vector<char> > bestSpliceMismatchCharVec_altVec;		

		vector<int> donerMapCumulativeMismatchNumVec;
		vector<int> donerMapMismatchPosVec;
		vector<char> donerMapMismatchCharVec;
		vector<int> acceptorMapCumulativeMismatchNumVec;
		vector<int> acceptorMapMismatchPosVec;
		vector<char> acceptorMapMismatchCharVec;
		int max_extension_forward_doner;
		int max_extension_backward_acceptor;
		vector< pair< int, pair<int, int> > > spliceSiteVec; // < < splice_site = forwared_doner_# , mismatch # >, flank_string_case >
		int bestSpliceSiteIndex;
		vector<int> bestSpliceMismatchPosVec;
		vector<char> bestSpliceMismatchCharVec;

		int bestSpliceSite_sdRef_altRef_index;
		int mismatchPos_interval_min;
public:
		// void copyMismatchPos2TargetVec(vector<int>& targetMismatchPosVec)
		// {
		// 	for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
		// 	{
		// 		targetMismatchPosVec.push_back(bestSpliceMismatchPosVec[tmp]);
		// 	}
		// }
		// void copyMismatchChar2TargetVec(vector<char>& targetMismatchCharVec)
		// {
		// 	for(int tmp = 0; tmp < bestSpliceMismatchCharVec.size(); tmp++)
		// 	{
		// 		targetMismatchCharVec.push_back(bestSpliceMismatchCharVec[tmp]);
		// 	}
		// }

		int getMinMismatchInterval(vector<int>& bestSpliceMismatchPosVec_sorted)
		{
				// cout << "before sorting the vector<int> ..." << endl;
				// for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				// {
				// 	cout << bestSpliceMismatchPosVec[tmp] << ",";
				// }
				// cout << endl;
				// sort(bestSpliceMismatchPosVec.begin(), bestSpliceMismatchPosVec.end());
				// cout << "after sorting ..." << endl;
				// for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				// {
				// 	cout << bestSpliceMismatchPosVec[tmp] << ",";
				// }
				// cout << endl;
				for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				{
					bestSpliceMismatchPosVec_sorted.push_back(bestSpliceMismatchPosVec[tmp]);
				}				
				sort(bestSpliceMismatchPosVec_sorted.begin(), bestSpliceMismatchPosVec_sorted.end());

				int tmpMinMismatchPosInterval = 100;
				for(int tmp = 0; tmp < bestSpliceMismatchPosVec_sorted.size()-1; tmp++)
				{
					int tmpIndex_1 = tmp;
					int tmpIndex_2 = tmp + 1;
					int tmpMismatchPos_1 = bestSpliceMismatchPosVec_sorted[tmpIndex_1];
					int tmpMismatchPos_2 = bestSpliceMismatchPosVec_sorted[tmpIndex_2];
					int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
						tmpMinMismatchPosInterval = tmpMismatchPosInterval;
				}
				return tmpMinMismatchPosInterval;
		}

		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> bestSpliceMismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(bestSpliceMismatchPosVec_sorted);
			//cout << "minimumMismatchGap: " << minimumMismatchGap << endl;
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < bestSpliceMismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = bestSpliceMismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = bestSpliceMismatchPosVec_sorted[tmpIndex_2];
				int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
				if(tmpMismatchPosInterval < minimumMismatchGap)
				{
					cout << "some thing wrong in getMinMismatchInterval_2" << endl;
					exit(1);
				}
				else if(tmpMismatchPosInterval == minimumMismatchGap)
				{	
					if(!checkedMinimumGap_bool)
						checkedMinimumGap_bool = true;
					else
						return tmpMismatchPosInterval;
				}
				else
				{
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval_2)
						tmpMinMismatchPosInterval_2 = tmpMismatchPosInterval;
				}
			}
			return tmpMinMismatchPosInterval_2;
		}

		bool filterInterResults(int maxMismatchNum_semi_noncanonical,
			int bestSpliceSiteIndex)
		{
			//cout << "start to filter interResults ..." << endl;
			int mismatchPosVecSize = bestSpliceMismatchPosVec.size();
			//cout << "mismatchPosVecSize: " << mismatchPosVecSize << endl;
			if(mismatchPosVecSize >= 3)
			{
				int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
				//cout << "tmpMinMismatchPosInterval_2: " << tmpMinMismatchPosInterval << endl;
				if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
					return false;
			}

			int tmpSJcase = (spliceSiteVec[bestSpliceSiteIndex].second).second;
			if(tmpSJcase < 5)
			{
				if(mismatchPosVecSize > maxMismatchNum_semi_noncanonical)
					return false;
			}
			return true;
		}

		FixDoubleAnchor_Splice_MultiRef_Info()
		{
			bestSpliceSiteIndex = -1;
			mismatchPos_interval_min = 3;
			bestSpliceSite_sdRef_altRef_index = -1;
			//indel_length_max = INDEL_BESIDE_SJ_LEN_MAX;
		}
		
		int returnBestSpliceMismatchPosVecSize()
		{
			return bestSpliceMismatchPosVec.size();
		}
		int returnBestSpliceMismatchCharVecSize()
		{
			return bestSpliceMismatchCharVec.size();
		}

		int returnMismatchPosInRead(int index_bestSpliceMismatchPosVec)
		{
			return bestSpliceMismatchPosVec[index_bestSpliceMismatchPosVec];
		}
		char returnMismatchChar(int index_bestSpliceMismatchCharVec)
		{
			return bestSpliceMismatchCharVec[index_bestSpliceMismatchCharVec];
		}

		void generateBestSpliceMismatchVec(int toFixSeqLength, int startLocInRead)
		{
			if(STORE_MISMATCH_CHA)
				this->generateBestSpliceMismatchVec_Pos_Char(toFixSeqLength, startLocInRead);
			else
				this->generateBestSpliceMismatchVec_Pos(toFixSeqLength, startLocInRead);
		}

		void generateBestSpliceMismatchVec_Pos(int toFixSeqLength, int startLocInRead)
		{
			int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
			//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
			for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = donerMapMismatchPosVec[tmp];
				if(tmpMismatchPos <= bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				}
			}

			for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
				if(tmpMismatchPos > bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				}
			}
		}

		void generateBestSpliceMismatchVec_Pos_Char(int toFixSeqLength, int startLocInRead)
		{
			int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
			//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
			for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = donerMapMismatchPosVec[tmp];
				if(tmpMismatchPos <= bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
					bestSpliceMismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
				}
			}

			for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
				if(tmpMismatchPos > bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
					bestSpliceMismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
				}
			}
		}

		string returnDonerAndAcceptorMapCumulativeMismatchNumVecStr()
		{
			string tmpStr = "donerMapCumulativeMismatchNumVec:\n";
			for(int tmp = 0; tmp < donerMapCumulativeMismatchNumVec.size(); tmp++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(donerMapCumulativeMismatchNumVec[tmp]) + ", ";
			}
			tmpStr = tmpStr + "\nacceptorMapCumulativeMismatchNumVec:\n";
			for(int tmp = 0; tmp < acceptorMapCumulativeMismatchNumVec.size(); tmp++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(acceptorMapCumulativeMismatchNumVec[tmp]) + ", ";
			}
			tmpStr += "\n";
			return tmpStr;
		}

		string returnSpliceSiteVecStr()
		{
			string tmpStr = "spliceSiteVec:\n";
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				tmpStr = tmpStr + "doner_length: " + int_to_str(spliceSiteVec[tmp].first) + " mismatch#: " 
					+ int_to_str((spliceSiteVec[tmp].second).first) + " flank_string_case: " 
					+ int_to_str((spliceSiteVec[tmp].second).second) + "\n";
			}
			return tmpStr;
		}
		int returnComparedPenalty(int mismatchPenalty,
					int semiCanonical_penalty, int nonCanonical_penalty)
		{
			int tmpMismatchNum_bestSplice = this->returnBestSplice_mismatchNum();
			int tmpFlankStringCase_bestSplice = this->returnBestSplice_flankString();
			int tmpMismatch_penalty = tmpMismatchNum_bestSplice *  mismatchPenalty;
			int tmpSpliceTypePenalty;
			if(tmpFlankStringCase_bestSplice >= 5)
			{
				tmpSpliceTypePenalty = 0;
			}
			else if(tmpFlankStringCase_bestSplice == 0)
			{
				tmpSpliceTypePenalty = nonCanonical_penalty;
			}
			else
				tmpSpliceTypePenalty = semiCanonical_penalty;
			return tmpMismatch_penalty + tmpSpliceTypePenalty;
		}
		int returnBestSplice_flankString()
		{
			int flankString_case = spliceSiteVec[bestSpliceSiteIndex].second.second;
			return flankString_case;
		}
		bool returnBestSplice_canonicalOrNot()
		{
			return this->SJ_canonicalOrNot(this->returnBestSplice_flankString());
		}
		int returnBestSplice_mismatchNum()
		{
			int bestSJ_mismatchNum = spliceSiteVec[bestSpliceSiteIndex].second.first;
			return bestSJ_mismatchNum;
		}

		int returnBestSplice_prefixMatchLength()
		{
			return spliceSiteVec[bestSpliceSiteIndex].first;
		}

		int returnFlankStringCase(const string& flank_string)   
		{
			if(flank_string == "ATAC")
			{
				return 1;
			}
			else if(flank_string == "GTAT")
			{
				return 2;
			}
			else if(flank_string == "CTGC")
			{
				return 3;
			}
			else if(flank_string == "GCAG")
			{
				return 4;
			}
			else if(flank_string == "GTAG")
			{
				return 5;
			}
			else if(flank_string == "CTAC")
			{
				return 6;
			}
			else
			{
				return 0;
			}
		}

		bool fixSpliceResultConfident()
		{
			if((this->returnBestSplice_mismatchNum() == 0)
				||((this->returnBestSplice_canonicalOrNot())&&(this->returnBestSplice_mismatchNum() == 1)))
			{
				return true;
			}
			else
			{
				return false;
			}
		}		

		bool detectBestSpliceSite_sd_includeSNPhashInfoVec(//_prefer_canonical_lessMismatch
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum,
			bool annotation_provided_bool, 
			bool Do_annotation_only_bool, Annotation_Info* annotationInfo, 
			vector< SNPhash_Info >& snpHashInfoVec)	
		{
			// apply to reference genome
			this->scanGenomeAndReadSeq_sd(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
				readSeq_inProcess, indexInfo, chrNameInt, max_allowed_mismatchNum);				
			this->generateSpliceSiteVec_sd(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				indexInfo, chrNameInt, max_allowed_mismatchNum, annotation_provided_bool, Do_annotation_only_bool, annotationInfo);			
			double bestSpliceSite_penalty_sd;
			bestSpliceSiteIndex_sd = this->selectBestSpliceSite_returnPenalty_sd(bestSpliceSite_penalty);

			// apply to multiple alternative references
			this->scanGenomeAndReadSeq_includeSNPhashInfoVec(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
				readSeq_inProcess, indexInfo, chrNameInt, max_allowed_mismatchNum, snpHashInfoVec);	
			this->generateSpliceSiteVec_includeSNPhashInfoVec(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				indexInfo, chrNameInt, max_allowed_mismatchNum, annotation_provided_bool, Do_annotation_only_bool, 
				annotationInfo, snpHashInfoVec);
			int tmpAltReferenceNum = snpHashInfoVec.size();
			vector<double> bestSpliceSite_penalty_altVec;
			for(int tmpAltReferenceIndex = 0; tmpAltReferenceIndex < tmpAltReferenceNum; tmpAltReferenceIndex ++)
			{	
				double tmpBestSpliceSite_penalty_alt
				int tmpBestSpliceSiteIndex_alt = this->selectBestSpliceSite_returnPenalty_alt(
					tmpAltReferenceIndex, tmpBestSpliceSite_penalty_alt);
				bestSpliceSiteIndex_altVec.push_back(tmpBestSpliceSiteIndex_altReferene);
				bestSpliceSite_penalty_altVec.push_back(tmpBestSpliceSite_penalty_altReference);
			}
			
			// do selection between sd and altVec
			vector<int> bestSpliceSiteIndexVec_candi;
			vector<double> bestSpliceSitePenaltyVec_candi;
			bestSpliceSiteIndexVec_candi.push_back(bestSpliceSiteIndex_sd);
			bestSpliceSitePenaltyVec_candi.push_back(bestSpliceSite_penalty_sd);
			for(int tmp = 0; tmp < tmpAltReferenceNum; tmp++)
			{
				bestSpliceSiteIndexVec_candi.push_back(bestSpliceSiteIndex_altVec[tmp]);
				bestSpliceSitePenaltyVec_candi.push_back(bestSpliceSite_penalty_altVec[tmp]);
			}

			int bestSpliceSiteIndexVec_candi_size = bestSpliceSiteIndexVec_candi.size();
			double currentBestSpliceSite_penalty = 100000.0;
			int currentBestSpliceSite_index = -1;
			bestSpliceSite_sdRef_altRef_index = -1; // standard reference == 0; tmpAltReferenceNum >= alterReference > 0
			for(int tmp = 0; tmp < tmpAltReferenceNum + 1; tmp++)
			{
				int tmpBestSpliceIndex_candi = bestSpliceSiteIndexVec_candi[tmp];
				if(tmpBestSpliceIndex_candi >= 0)
				{
					double tmpPenalty_candi = bestSpliceSitePenaltyVec_candi[tmp];
					if(tmpPenalty_candi < currentBestSpliceSite_penalty)
					{
						currentBestSpliceSite_index = tmpBestSpliceIndex_candi;
						currentBestSpliceSite_penalty = tmpPenalty_candi;
						bestSpliceSite_sdRef_altRef_index = tmp;
					}
				}
			}

			int maxMismatchNum_semi_noncanonical = max_allowed_mismatchNum;
			if((currentBestSpliceSite_index >= 0)&&(bestSpliceSite_sdRef_altRef_index >= 0))
			{
				if(bestSpliceSite_sdRef_altRef_index == 0) // cp sd
				{
					// vector<int> donerMapCumulativeMismatchNumVec;
					int donerMapCumulativeMismatchNumVec_sd_size = donerMapCumulativeMismatchNumVec_sd.size();
					for(int tmp = 0; tmp < donerMapCumulativeMismatchNumVec_sd_size; tmp ++)
						donerMapCumulativeMismatchNumVec.push_back(donerMapCumulativeMismatchNumVec_sd[tmp]);
					// vector<int> donerMapMismatchPosVec;
					int donerMapMismatchPosVec_sd_size = donerMapMismatchPosVec_sd.size();
					for(int tmp = 0; tmp < donerMapMismatchPosVec_sd_size; tmp++)
						donerMapMismatchPosVec.push_back(donerMapMismatchPosVec_sd[tmp]);
					// vector<char> donerMapMismatchCharVec;
					int donerMapMismatchCharVec_sd_size = donerMapMismatchCharVec_sd.size();
					for(int tmp = 0; tmp < donerMapMismatchCharVec_sd_size; tmp++)
						donerMapMismatchCharVec.push_back(donerMapMismatchCharVec_sd[tmp]);
					// vector<int> acceptorMapCumulativeMismatchNumVec;
					int acceptorMapCumulativeMismatchNumVec_sd_size = acceptorMapCumulativeMismatchNumVec_sd.size();
					for(int tmp = 0; tmp < acceptorMapCumulativeMismatchNumVec_sd_size; tmp++)
						acceptorMapCumulativeMismatchNumVec.push_back(acceptorMapCumulativeMismatchNumVec_sd[tmp]);
					// vector<int> acceptorMapMismatchPosVec;
					int acceptorMapMismatchPosVec_sd_size = acceptorMapMismatchPosVec_sd.size();
					for(int tmp = 0; tmp < acceptorMapMismatchPosVec_sd_size; tmp++)
						acceptorMapMismatchPosVec.push_back(acceptorMapMismatchPosVec_sd[tmp]);
					// vector<char> acceptorMapMismatchCharVec;
					int acceptorMapMismatchCharVec_sd_size = acceptorMapMismatchCharVec_sd.size();
					for(int tmp = 0; tmp < acceptorMapMismatchCharVec_sd_size; tmp++)
						acceptorMapMismatchCharVec.push_back(acceptorMapMismatchCharVec_sd[tmp]);
					// int max_extension_forward_doner;
					max_extension_forward_doner = max_extension_forward_doner_sd;
					// int max_extension_backward_acceptor;
					max_extension_backward_acceptor = max_extension_backward_acceptor_sd;
					// vector< pair< int, pair<int, int> > > spliceSiteVec; 
					int spliceSiteVec_sd_size = spliceSiteVec_sd.size();
					for(int tmp = 0; tmp < spliceSiteVec_sd_size; tmp++)
						spliceSiteVec.push_back(spliceSiteVec_sd[tmp]);
					// int bestSpliceSiteIndex;
					bestSpliceSiteIndex = bestSpliceSiteIndex_sd;
					// vector<int> bestSpliceMismatchPosVec;
					int bestSpliceMismatchPosVec_sd_size = bestSpliceMismatchPosVec_sd.size();
					for(int tmp = 0; tmp < bestSpliceMismatchPosVec_sd_size; tmp++)
						bestSpliceMismatchPosVec.push_back(bestSpliceMismatchPosVec_sd[tmp]);
					// vector<char> bestSpliceMismatchCharVec;
					int bestSpliceMismatchCharVec_sd_size = bestSpliceMismatchCharVec_sd.size();
					for(int tmp = 0; tmp < bestSpliceMismatchCharVec_sd_size; tmp++)
						bestSpliceMismatchCharVec.push_back(bestSpliceMismatchCharVec_sd[tmp]);
				}
				else // cp alt
				{
					int tmpIndexInAltVec = bestSpliceSite_sdRef_altRef_index - 1;
					// vector<int> donerMapCumulativeMismatchNumVec;
					int donerMapCumulativeMismatchNumVec_alt_size = donerMapCumulativeMismatchNumVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < donerMapCumulativeMismatchNumVec_alt_size; tmp ++)
						donerMapCumulativeMismatchNumVec.push_back((donerMapCumulativeMismatchNumVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<int> donerMapMismatchPosVec;
					int donerMapMismatchPosVec_alt_size = donerMapMismatchPosVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < donerMapMismatchPosVec_alt_size; tmp++)
						donerMapMismatchPosVec.push_back((donerMapMismatchPosVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<char> donerMapMismatchCharVec;
					int donerMapMismatchCharVec_alt_size = donerMapMismatchCharVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < donerMapMismatchCharVec_sd_size; tmp++)
						donerMapMismatchCharVec.push_back((donerMapMismatchCharVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<int> acceptorMapCumulativeMismatchNumVec;
					int acceptorMapCumulativeMismatchNumVec_alt_size = acceptorMapCumulativeMismatchNumVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < acceptorMapCumulativeMismatchNumVec_alt_size; tmp++)
						acceptorMapCumulativeMismatchNumVec.push_back((acceptorMapCumulativeMismatchNumVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<int> acceptorMapMismatchPosVec;
					int acceptorMapMismatchPosVec_alt_size = acceptorMapMismatchPosVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < acceptorMapMismatchPosVec_alt_size; tmp++)
						acceptorMapMismatchPosVec.push_back((acceptorMapMismatchPosVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<char> acceptorMapMismatchCharVec;
					int acceptorMapMismatchCharVec_alt_size = acceptorMapMismatchCharVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < acceptorMapMismatchCharVec_alt_size; tmp++)
						acceptorMapMismatchCharVec.push_back((acceptorMapMismatchCharVec_alt[tmpIndexInAltVec])[tmp]);
					// int max_extension_forward_doner;
					max_extension_forward_doner = max_extension_forward_doner_alt[tmpIndexInAltVec];
					// int max_extension_backward_acceptor;
					max_extension_backward_acceptor = max_extension_backward_acceptor_alt[tmpIndexInAltVec];
					// vector< pair< int, pair<int, int> > > spliceSiteVec; 
					int spliceSiteVec_alt_size = spliceSiteVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < spliceSiteVec_alt_size; tmp++)
						spliceSiteVec.push_back((spliceSiteVec_alt[tmpIndexInAltVec])[tmp]);
					// int bestSpliceSiteIndex;
					bestSpliceSiteIndex = bestSpliceSiteIndex_alt[tmpIndexInAltVec];
					// vector<int> bestSpliceMismatchPosVec;
					int bestSpliceMismatchPosVec_alt_size = bestSpliceMismatchPosVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < bestSpliceMismatchPosVec_alt_size; tmp++)
						bestSpliceMismatchPosVec.push_back((bestSpliceMismatchPosVec_alt[tmpIndexInAltVec])[tmp]);
					// vector<char> bestSpliceMismatchCharVec;
					int bestSpliceMismatchCharVec_alt_size = bestSpliceMismatchCharVec_alt[tmpIndexInAltVec].size();
					for(int tmp = 0; tmp < bestSpliceMismatchCharVec_alt_size; tmp++)
						bestSpliceMismatchCharVec.push_back((bestSpliceMismatchCharVec_alt[tmpIndexInAltVec])[tmp]);
				}

				this->generateBestSpliceMismatchVec(toFixSeqLocInRead_end-toFixSeqLocInRead_start+1, toFixSeqLocInRead_start);
				bool filter_bool = filterInterResults(maxMismatchNum_semi_noncanonical, currentBestSpliceSite_index);
				if(filter_bool)
					return true;
				else
					return false;
			}
			else
				return false;
		}

		void generateBestSpliceMismatchVec_withSdAltRefIndex(int toFixSeqLength, int startLocInRead, int sdRef_altRef_index)
		{
			if(sdRef_altRef_index)
			{	
				int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
				//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
				for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = donerMapMismatchPosVec[tmp];
					if(tmpMismatchPos <= bestSplice_doner_length + startLocInRead - 1)
					{
						bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
						bestSpliceMismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
					}
				}

				for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
					if(tmpMismatchPos > bestSplice_doner_length + startLocInRead - 1)
					{
						bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
						bestSpliceMismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
					}
				}
			}
			else
			{

			}
		}


		double returnSJpenalty(int flank_string_case)
		{
			if((flank_string_case == 5)||(flank_string_case == 6)) 
				return CANONICAL_SJ_PENALTY;
			else if((flank_string_case >= 1)&&(flank_string_case <= 4))
				return SEMICANONICAL_SJ_PENALTY;
			else if(flank_string_case == 0)
				return NONCANONICAL_SJ_PENALTY;
			else
			{
				cout << "incorrect flank_string_case in FixDoubleAnchorSplice_Info.h" << endl;
				exit(1);
			}
		}

		
		bool SJ_canonicalOrNot(int flank_string_case)
		{
			if((flank_string_case == 5)||(flank_string_case == 6)) 
				return true;
			else
				return false;
		}

		bool SJ_semicanonicalOrNot(int flank_string_case)
		{
			if((flank_string_case >= 1)&&(flank_string_case <= 4)) 
				return true;
			else
				return false;
		}

		bool SJ_noncanonicalOrNot(int flank_string_case)
		{
			if(flank_string_case == 0) 
				return true;
			else
				return false;
		}

		int selectBestSpliceSite() 
		// penalty = mismatchNum + SJ_penalty(canon=0, semi=1.5, noncanon=2);
		{
			double currentBestSJ_penalty = 100000.0;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				double tmpMismatchNum = (double)spliceSiteVec[tmp].second.first;
				double tmpSJpenalty = this->returnSJpenalty(spliceSiteVec[tmp].second.second);
				double tmpPenalty = tmpMismatchNum + tmpSJpenalty;
				if(tmpPenalty < currentBestSJ_penalty)
				{
					currentBestSJ_penalty = tmpPenalty;
					currentBestSJ_index = tmp;
				}
			}			
			return currentBestSJ_index;
		}

		int selectBestSpliceSite_returnPenalty_sd(double& tmpPenalty) 
		// penalty = mismatchNum + SJ_penalty(canon=0, semi=1.5, noncanon=2);
		{
			double currentBestSJ_penalty = 100000.0;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec_sd.size(); tmp++)
			{
				double tmpMismatchNum = (double)spliceSiteVec_sd[tmp].second.first;
				double tmpSJpenalty = this->returnSJpenalty(spliceSiteVec_sd[tmp].second.second);
				double tmpPenalty = tmpMismatchNum + tmpSJpenalty;
				if(tmpPenalty < currentBestSJ_penalty)
				{
					currentBestSJ_penalty = tmpPenalty;
					currentBestSJ_index = tmp;
				}
			}
			tmpPenalty = currentBestSJ_penalty;
			return currentBestSJ_index;
		}

		int selectBestSpliceSite_returnPenalty_alt(int tmpAltReferenceIndex, double& tmpPenalty) 
		// penalty = mismatchNum + SJ_penalty(canon=0, semi=1.5, noncanon=2);
		{
			double currentBestSJ_penalty = 100000.0;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec_altVec[tmpAltReferenceIndex].size(); tmp++)
			{
				double tmpMismatchNum = (double)(spliceSiteVec_altVec[tmpAltReferenceIndex])[tmp].second.first;
				double tmpSJpenalty = this->returnSJpenalty((spliceSiteVec_altVec[tmpAltReferenceIndex])[tmp].second.second);
				double tmpPenalty = tmpMismatchNum + tmpSJpenalty;
				if(tmpPenalty < currentBestSJ_penalty)
				{
					currentBestSJ_penalty = tmpPenalty;
					currentBestSJ_index = tmp;
				}
			}
			tmpPenalty = currentBestSJ_penalty;
			return currentBestSJ_index;
		}

		void generateSpliceSiteVec_noAnnotationProvided_sd(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner_sd; tmp_forwared_doner ++ )
			{
				int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
				if(tmp_backward_acceptor > max_extension_backward_acceptor_sd)
					continue;
				int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
				int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
				int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec_sd[index_donerMapCumulativeMismatchNumVec];
				int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec_sd[index_acceptorMapCumulativeMismatchNumVec];
				int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;
				if(tmp_mismatch_sum <= max_allowed_mismatchNum)
				{
					int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
					int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;
					string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
					int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
					spliceSiteVec_sd.push_back(pair< int, pair<int, int> >(
						tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
				}
			}
			return;
		}

		void generateSpliceSiteVec_noAnnotationProvided_includeSNPhashInfoVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum,
			vector<SNPhash_Info>& snpHashInfoVec)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			int tmpAltReferenceNum = snpHashInfoVec.size();
			for(int tmpAltReferenceIndex = 0; tmpAltReferenceIndex < tmpAltReferenceNum; tmpAltReferenceIndex ++)
			{
				vector< pair< int, pair<int, int> > > tmpSpliceSiteVec_alt;
				spliceSiteVec_altVec.push_back(tmpSpliceSiteVec_alt);
				for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner_altVec[tmpAltReferenceIndex]; tmp_forwared_doner ++ )
				{
					int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
					if(tmp_backward_acceptor > max_extension_backward_acceptor_altVec[tmpAltReferenceIndex])
						continue;
					int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
					int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
					int tmp_doner_mismatch = (donerMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex])[index_donerMapCumulativeMismatchNumVec];
					int tmp_acceptor_mismatch = (acceptorMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex])[index_acceptorMapCumulativeMismatchNumVec];
					int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;

					if(tmp_mismatch_sum <= max_allowed_mismatchNum)
					{
						int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
						int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;
						string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
						string tmpSNPbaseInFlankString_1;
						string tmpSNPbaseInFlankString_2;
						string tmpSNPbaseInFlankString_3;
						string tmpSNPbaseInFlankString_4;
						bool tmpSearchBaseInFlankString_1_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_doner + 1, tmpSNPbaseInFlankString_1);
						if(tmpSearchBaseInFlankString_1_bool)
							tmp_flank_string.replace(0, 1, tmpSNPbaseInFlankString_1);
						bool tmpSearchBaseInFlankString_2_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_doner + 2, tmpSNPbaseInFlankString_2);
						if(tmpSearchBaseInFlankString_2_bool)
							tmp_flank_string.replace(1, 1, tmpSNPbaseInFlankString_2);
						bool tmpSearchBaseInFlankString_3_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_acceptor - 2, tmpSNPbaseInFlankString_3);
						if(tmpSearchBaseInFlankString_3_bool)
							tmp_flank_string.replace(2, 1, tmpSNPbaseInFlankString_3);
						bool tmpSearchBaseInFlankString_4_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_acceptor - 1, tmpSNPbaseInFlankString_4);
						if(tmpSearchBaseInFlankString_4_bool)
							tmp_flank_string.replace(3, 1, tmpSNPbaseInFlankString_4);
						int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
						spliceSiteVec_altVec[tmpAltReferenceIndex].push_back(pair< int, pair<int, int> >(tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
					}
				}
			}
			return;
		}

		void generateSpliceSiteVec_doAnnotationOnly_sd(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			Annotation_Info* annotationInfo)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner_sd; tmp_forwared_doner ++ )
			{
				int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
				if(tmp_backward_acceptor > max_extension_backward_acceptor_sd)
					continue;
				// search in annotated SJs
				int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
				int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;				
				bool foundInAnnotation_bool = annotationInfo->SJfoundInAnnotation(chrNameInt, tmpChromPos_doner,
					tmpChromPos_acceptor);
				if(!foundInAnnotation_bool)
					continue;
				int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
				int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
				int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec_sd[index_donerMapCumulativeMismatchNumVec];
				int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec_sd[index_acceptorMapCumulativeMismatchNumVec];
				int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;
				if(tmp_mismatch_sum <= max_allowed_mismatchNum)
				{
					string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
					int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
					spliceSiteVec_sd.push_back(pair< int, pair<int, int> >(
						tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
				}
			}
			return;
		}

		void generateSpliceSiteVec_doAnnotationOnly_includeSNPhashInfoVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			Annotation_Info* annotationInfo, vector<SNPhash_Info>& snpHashInfoVec)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			int altReferenceSeqNum = snpHashInfoVec.size();
			for(int tmpAltReferenceIndex = 0; tmpAltReferenceIndex < altReferenceSeqNum; tmpAltReferenceIndex ++)
			{
				vector< pair< int, pair<int, int> > > tmpSpliceSiteVec_alt;
				spliceSiteVec_altVec.push_back(tmpSpliceSiteVec_alt);
				for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner_altVec[tmpAltReferenceIndex]; tmp_forwared_doner ++ )
				{
					int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
					if(tmp_backward_acceptor > max_extension_backward_acceptor_altVec[tmpAltReferenceIndex])
						continue;
					// search in annotated SJs
					int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
					int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;				
					bool foundInAnnotation_bool = annotationInfo->SJfoundInAnnotation(chrNameInt, tmpChromPos_doner,
						tmpChromPos_acceptor);
					if(!foundInAnnotation_bool)
						continue;
					int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
					int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
					int tmp_doner_mismatch = (donerMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex])[index_donerMapCumulativeMismatchNumVec];
					int tmp_acceptor_mismatch = (acceptorMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex])[index_acceptorMapCumulativeMismatchNumVec];
					int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;
					if(tmp_mismatch_sum <= max_allowed_mismatchNum)
					{
						string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
						string tmpSNPbaseInFlankString_1;
						string tmpSNPbaseInFlankString_2;
						string tmpSNPbaseInFlankString_3;
						string tmpSNPbaseInFlankString_4;
						bool tmpSearchBaseInFlankString_1_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_doner + 1, tmpSNPbaseInFlankString_1);
						if(tmpSearchBaseInFlankString_1_bool)
							tmp_flank_string.replace(0, 1, tmpSNPbaseInFlankString_1);
						bool tmpSearchBaseInFlankString_2_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_doner + 2, tmpSNPbaseInFlankString_2);
						if(tmpSearchBaseInFlankString_2_bool)
							tmp_flank_string.replace(1, 1, tmpSNPbaseInFlankString_2);
						bool tmpSearchBaseInFlankString_3_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_acceptor - 2, tmpSNPbaseInFlankString_3);
						if(tmpSearchBaseInFlankString_3_bool)
							tmp_flank_string.replace(2, 1, tmpSNPbaseInFlankString_3);
						bool tmpSearchBaseInFlankString_4_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
							chrNameInt, tmpChromPos_acceptor - 1, tmpSNPbaseInFlankString_4);
						if(tmpSearchBaseInFlankString_4_bool)
							tmp_flank_string.replace(3, 1, tmpSNPbaseInFlankString_4);
						int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
						spliceSiteVec_altVec[tmpAltReferenceIndex].push_back(pair< int, pair<int, int> >(
							tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
					}
				}
			}
			return;
		}		

		void generateSpliceSiteVec_sd(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			bool annotation_provided_bool, bool Do_annotation_only_bool, 
			Annotation_Info* annotationInfo)
		{
			if(annotation_provided_bool)
			{
				if(Do_annotation_only_bool)
				{
					this->generateSpliceSiteVec_doAnnotationOnly_sd(
						toFixSeqLocInRead_start, toFixSeqLocInRead_end,
						toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
						indexInfo, chrNameInt, max_allowed_mismatchNum, annotationInfo);
				}
				else
				{}
			}	
			else
			{
				this->generateSpliceSiteVec_noAnnotationProvided_sd(
					toFixSeqLocInRead_start, toFixSeqLocInRead_end,
					toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
					indexInfo, chrNameInt, max_allowed_mismatchNum);
			}
		}

		void generateSpliceSiteVec_includeSNPhashInfoVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			bool annotation_provided_bool, bool Do_annotation_only_bool, 
			Annotation_Info* annotationInfo, vector<SNPhash_Info>& snpHashInfoVec)
		{
			if(annotation_provided_bool)
			{
				if(Do_annotation_only_bool)
				{
					// for now, does not include SNPhashInfo for generateSpliceSiteVec_doAnnotationOnly
					this->generateSpliceSiteVec_doAnnotationOnly_includeSNPhashInfoVec(
						toFixSeqLocInRead_start, toFixSeqLocInRead_end,
						toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
						indexInfo, chrNameInt, max_allowed_mismatchNum, annotationInfo);
				}
				else
				{}
			}	
			else
			{
				this->generateSpliceSiteVec_noAnnotationProvided_includeSNPhashInfoVec(
					toFixSeqLocInRead_start, toFixSeqLocInRead_end,
					toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
					indexInfo, chrNameInt, max_allowed_mismatchNum, snpHashInfoVec);
			}
		}

		void scanGenomeAndReadSeq_sd(int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)	
		{
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			max_extension_forward_doner_sd = toFixSeqLength;
			max_extension_backward_acceptor_sd = toFixSeqLength;

			int tmpDonerMapCumulativeMismatchNum_sd = 0;
			donerMapCumulativeMismatchNumVec_sd.push_back(0);
			for(int tmp = 0; tmp < toFixSeqLength; tmp++)
			{
				char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_doner + tmp);
				if(readSeq_inProcess.at(toFixSeqLocInRead_start - 1 + tmp) != charInRef)
				{
					tmpDonerMapCumulativeMismatchNum_sd ++;
					donerMapMismatchPosVec_sd.push_back(tmp+toFixSeqLocInRead_start);
					donerMapMismatchCharVec_sd.push_back(charInRef);
				}
				if(tmpDonerMapCumulativeMismatchNum_sd > max_allowed_mismatchNum)
				{
					max_extension_forward_doner_sd = tmp;
					break;
				}
				donerMapCumulativeMismatchNumVec_sd.push_back(tmpDonerMapCumulativeMismatchNum_sd);
			}

			int tmpAcceptorMapCumulativeMismatchNum_sd = 0;
			acceptorMapCumulativeMismatchNumVec_sd.push_back(0);
			for(int tmp = 0; tmp < toFixSeqLength; tmp++)
			{
				char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_acceptor - tmp);
				if(readSeq_inProcess.at(toFixSeqLocInRead_end - tmp - 1) != charInRef)
				{
					tmpAcceptorMapCumulativeMismatchNum_sd ++;
					acceptorMapMismatchPosVec_sd.push_back(toFixSeqLocInRead_end - tmp);
					acceptorMapMismatchCharVec_sd.push_back(charInRef);
				}
				if(tmpAcceptorMapCumulativeMismatchNum_sd > max_allowed_mismatchNum)
				{
					max_extension_backward_acceptor_sd = tmp;
					break;
				}
				acceptorMapCumulativeMismatchNumVec_sd.push_back(tmpAcceptorMapCumulativeMismatchNum_sd);
			}
			return;
		}

		void scanGenomeAndReadSeq_includeSNPhashInfoVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end, int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum,
			vector<SNPhash_Info>& snpHashInfoVec)	
		{
			int altReferenceSeqNum = snpHashInfoVec.size();
			//cout << endl << "start to scan scanGenomeAndReadSeq " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;

			for(int tmpAltReferenceIndex = 0; tmpAltReferenceIndex < altReferenceSeqNum; tmpAltReferenceIndex ++)
			{	
				max_extension_forward_doner_altVec.push_back(toFixSeqLength);				
				vector<int> tmpCumulativeMismatchNumVec_alt;
				donerMapCumulativeMismatchNumVec_altVec.push_back(tmpCumulativeMismatchNumVec_alt);
				int tmpDonerMapCumulativeMismatchNum_tmpAltReference = 0;
				donerMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex].push_back(0);
				for(int tmp = 0; tmp < toFixSeqLength; tmp++)
				{
					char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_doner + tmp);
					string tmpSNPbaseStr;
					bool searchInTmpSNPhashInfo_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
						chrNameInt, toFixSeqMapPos_doner + tmp, tmpSNPbaseStr);
					if(searchInTmpSNPhashInfo_bool)
						charInRef = tmpSNPbaseStr.at[0];

					if(readSeq_inProcess.at(toFixSeqLocInRead_start - 1 + tmp) != charInRef)
					{
						tmpDonerMapCumulativeMismatchNum_tmpAltReference ++;
						donerMapMismatchPosVec_altVec[tmpAltReferenceIndex].push_back(tmp+toFixSeqLocInRead_start);
						donerMapMismatchCharVec_altVec[tmpAltReferenceIndex].push_back(charInRef);
					}
					if(tmpDonerMapCumulativeMismatchNum_tmpAltReference > max_allowed_mismatchNum)
					{
						max_extension_forward_doner_altVec[tmpAltReferenceIndex] = tmp;
						break;
					}
					donerMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex].push_back(tmpDonerMapCumulativeMismatchNum_tmpAltReference);
				}
			}

			for(int tmpAltReferenceIndex = 0; tmpAltReferenceIndex < altReferenceSeqNum; tmpAltReferenceIndex ++)
			{	
				max_extension_backward_acceptor_altVec.push_back(toFixSeqLength);
				vector<int> tmpCumulativeMismatchNumVec_alt;
				acceptorMapCumulativeMismatchNumVec_altVec.push_back(tmpCumulativeMismatchNumVec_alt);
				int tmpAcceptorMapCumulativeMismatchNum_tmpAltReference = 0;
				acceptorMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex].push_back(0);
				for(int tmp = 0; tmp < toFixSeqLength; tmp++)
				{
					char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_acceptor - tmp);
					string tmpSNPbaseStr;
					bool searchInTmpSNPhashInfo_bool = snpHashInfoVec[tmpAltReferenceIndex].searchAndReturnSNPbase(
						chrNameInt, toFixSeqMapPos_acceptor - tmp, tmpSNPbaseStr);
					if(searchInTmpSNPhashInfo_bool)
						charInRef = tmpSNPbaseStr.at[0];

					if(readSeq_inProcess.at(toFixSeqLocInRead_end - tmp - 1) != charInRef)
					{
						tmpAcceptorMapCumulativeMismatchNum_tmpAltReference ++;
						acceptorMapMismatchPosVec_altVec[tmpAltReferenceIndex].push_back(toFixSeqLocInRead_end - tmp);
						acceptorMapMismatchCharVec_altVec[tmpAltReferenceIndex].push_back(charInRef);
					}
					if(tmpAcceptorMapCumulativeMismatchNum_tmpAltReference > max_allowed_mismatchNum)
					{
						max_extension_backward_acceptor_altVec[tmpAltReferenceIndex] = tmp;
						break;
					}
					acceptorMapCumulativeMismatchNumVec_altVec[tmpAltReferenceIndex].push_back(tmpAcceptorMapCumulativeMismatchNum_tmpAltReference);
					//cout << "tmp: " << tmp << "tmpAcceptorMapCumulativeMismatchNum: " << tmpAcceptorMapCumulativeMismatchNum << endl;
				}
			}
			return;
			//int tmpSpliceSite = 0;
		}
};

#endif