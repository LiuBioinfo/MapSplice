all: alignAll annotation

OPTFLAGS = -O0 -g -rdynamic -std=c++0x
OUTPUTDIR = ./bin/

CFLAGS += $(OPTFLAGS) -fopenmp

mps: sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify \
	detection_globalMap generateFusionJuncInfoFromAdjustedBreakPointFile separateFusionJuncStrandedOrNot \
	filterFusionJunc_anchorSeqSimilarity filterNonStrandedFusionJunc_anchorSeqSimilarity \
	compareFusionJuncResultsWithGTFannotation filtering_supNum_repeat
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps_regular src/mps_arrayList.cpp
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_regular_circRNA src/mps_arrayList.cpp	
	g++ $(CFLAGS) -D MPS_FUSION_POST_NEW -o $(OUTPUTDIR)mps_regular_fusion src/mps_arrayList.cpp
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -D MPS_FUSION_POST_NEW -o $(OUTPUTDIR)mps_regular_circRNA_fusion src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps3_fusion_post_new src/otherProjects/fusion_post_new/fusion_post_main.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

rm_mps:
	rm -f $(OUTPUTDIR)mps_regular
	rm -f $(OUTPUTDIR)mps_regular_circRNA
	rm -f $(OUTPUTDIR)mps_regular_fusion
	rm -f $(OUTPUTDIR)mps_regular_circRNA_fusion
	rm -f $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm -f $(OUTPUTDIR)reportJunc_embeded
	rm -f $(OUTPUTDIR)removeFalseJuncSam
	rm -f $(OUTPUTDIR)mps3_fusion_post_new
	rm -f $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify
	rm -f $(OUTPUTDIR)detection_globalMap
	rm -f $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile
	rm -f $(OUTPUTDIR)separateFusionJuncStrandedOrNot
	rm -f $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation 
	rm -f $(OUTPUTDIR)filtering_supNum_repeat

mps_regular:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps_regular src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

mps_regular_circRNA:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_regular_circRNA src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

mps_regular_fusion: sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify \
	detection_globalMap generateFusionJuncInfoFromAdjustedBreakPointFile separateFusionJuncStrandedOrNot \
	filterFusionJunc_anchorSeqSimilarity filterNonStrandedFusionJunc_anchorSeqSimilarity \
	compareFusionJuncResultsWithGTFannotation filtering_supNum_repeat
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps_regular_fusion src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps3_fusion_post_new src/otherProjects/fusion_post_new/fusion_post_main.cpp

mps_regular_circRNA_fusion: sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify \
	detection_globalMap generateFusionJuncInfoFromAdjustedBreakPointFile separateFusionJuncStrandedOrNot \
	filterFusionJunc_anchorSeqSimilarity filterNonStrandedFusionJunc_anchorSeqSimilarity \
	compareFusionJuncResultsWithGTFannotation filtering_supNum_repeat
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_regular_circRNA_fusion src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

rm_mps_regular:
	rm -f $(OUTPUTDIR)mps_regular
	rm -f $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm -f $(OUTPUTDIR)reportJunc_embeded
	rm -f $(OUTPUTDIR)removeFalseJuncSam

rm_mps_regular_circRNA:
	rm -f $(OUTPUTDIR)mps_regular_circRNA
	rm -f $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm -f $(OUTPUTDIR)reportJunc_embeded
	rm -f $(OUTPUTDIR)removeFalseJuncSam

rm_mps_regular_fusion:
	rm -f $(OUTPUTDIR)mps_regular_fusion
	rm -f $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm -f $(OUTPUTDIR)reportJunc_embeded
	rm -f $(OUTPUTDIR)removeFalseJuncSam
	rm -f $(OUTPUTDIR)mps3_fusion_post_new
	rm -f $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify
	rm -f $(OUTPUTDIR)detection_globalMap
	rm -f $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile
	rm -f $(OUTPUTDIR)separateFusionJuncStrandedOrNot
	rm -f $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation 
	rm -f $(OUTPUTDIR)filtering_supNum_repeat

rm_mps_regular_circRNA_fusion:
	rm -f $(OUTPUTDIR)mps_regular_circRNA_fusion
	rm -f $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm -f $(OUTPUTDIR)reportJunc_embeded
	rm -f $(OUTPUTDIR)removeFalseJuncSam
	rm -f $(OUTPUTDIR)mps3_fusion_post_new
	rm -f $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify
	rm -f $(OUTPUTDIR)detection_globalMap
	rm -f $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile
	rm -f $(OUTPUTDIR)separateFusionJuncStrandedOrNot
	rm -f $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity
	rm -f $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation 
	rm -f $(OUTPUTDIR)filtering_supNum_repeat

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#mps:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)mps src/mps_arrayList.cpp
#	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
#	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
#	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

mps_circRNA:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_circRNA src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

mps_arrayList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps_arrayList src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

mps_arrayList_backSpliceIncluded:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_arrayList_backSpliceIncluded src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp

mps_backSpliceIncluded:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_backSpliceIncluded src/mps_arrayList.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp

mps_batch:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps_batch src/mps_batch.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_batch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp

mps_batch_backSpliceIncluded:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)mps_batch_backSpliceIncluded src/mps_batch.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_batch src/otherProjects/indexing/utils/build_snpMerIndex.cpp
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp

reportJunc_embeded:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reportJunc_embeded src/otherProjects/utilities/reportJunc_embeded.cpp

removeFalseJuncSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeFalseJuncSam src/otherProjects/utilities/removeFalseJuncSam.cpp

#mps_fusion_post

mps3_fusion_post_new: sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify \
	detection_globalMap generateFusionJuncInfoFromAdjustedBreakPointFile separateFusionJuncStrandedOrNot \
	filterFusionJunc_anchorSeqSimilarity filterNonStrandedFusionJunc_anchorSeqSimilarity \
	compareFusionJuncResultsWithGTFannotation filtering_supNum_repeat
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps3_fusion_post_new src/otherProjects/fusion_post_new/fusion_post_main.cpp

detection_globalMap:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detection_globalMap src/otherProjects/fusion_post_new/utils/detection_globalMap.cpp

generateFusionJuncInfoFromAdjustedBreakPointFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile src/otherProjects/fusion_post_new/utils/generateFusionJuncInfoFromAdjustedBreakPointFile.cpp

separateFusionJuncStrandedOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateFusionJuncStrandedOrNot src/otherProjects/fusion_post_new/utils/separateFusionJuncStrandedOrNot.cpp

filterFusionJunc_anchorSeqSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity src/otherProjects/fusion_post_new/utils/fusionJuncFiltering_anchorSeqSimilarity.cpp

filterNonStrandedFusionJunc_anchorSeqSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity src/otherProjects/fusion_post_new/utils/filterNonStrandedFusionJunc_anchorSeqSimilarity.cpp

compareFusionJuncResultsWithGTFannotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation src/otherProjects/fusion_post_new/utils/compareFusionJuncResultsWithGTFannotation.cpp

filtering_supNum_repeat:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filtering_supNum_repeat src/otherProjects/fusion_post_new/utils/filtering_supNum_repeat.cpp

convertMPS3fusionResults2bedpe:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertMPS3fusionResults2bedpe src/otherProjects/fusion_post_new/utils/convertMPS3fusionResults2bedpe.cpp

cmpFusionRes2gt:
	g++ $(CFLAGS) -o $(OUTPUTDIR)cmpFusionRes2gt src/otherProjects/MPS3_fusion/fusionReviewNAR/cmpFusionRes2gt.cpp

# seqQueryTree
#parseSeqOthelloResults_uniformSampleTranscript:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSeqOthelloResults_uniformSampleTranscript src/otherProjects/seqQueryTree/parseSeqOthelloResults_uniformSampleTranscript.cpp

#parseSbtResults_uniformSampleTranscript:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSbtResults_uniformSampleTranscript src/otherProjects/seqQueryTree/parseSbtResults_uniformSampleTranscript.cpp

#parseSbtAsResults_uniformSampleTranscript:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSbtAsResults_uniformSampleTranscript src/otherProjects/seqQueryTree/parseSbtAsResults_uniformSampleTranscript.cpp

#parseSplitSbtResults_uniformSampleTranscript:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSplitSbtResults_uniformSampleTranscript src/otherProjects/seqQueryTree/parseSplitSbtResults_uniformSampleTranscript.cpp

#mergeMultiSampleSailfishResults_uniformSampleTranscript:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeMultiSampleSailfishResults_uniformSampleTranscript src/otherProjects/seqQueryTree/mergeMultiSampleSailfishResults_uniformSampleTranscript.cpp

checkQueryResultsDiffSizeCategoryDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkQueryResultsDiffSizeCategoryDistribution src/otherProjects/seqQueryTree/checkQueryResultsDiffSizeCategoryDistribution.cpp

diffTwoUnsortedList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)diffTwoUnsortedList src/otherProjects/seqQueryTree/diffTwoUnsortedList.cpp

bam2fq2jf2kmer:
	g++ $(CFLAGS) -o $(OUTPUTDIR)bam2fq2jf2kmer src/otherProjects/seqQueryTree/bam2fq2jf2kmer.cpp

checkTwoSeqQueryToolDiff:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkTwoSeqQueryToolDiff src/otherProjects/seqQueryTree/checkTwoSeqQueryToolDiff.cpp

examineSeqOthelloKmerQueryCorrectness:
	g++ $(CFLAGS) -o $(OUTPUTDIR)examineSeqOthelloKmerQueryCorrectness src/otherProjects/seqQueryTree/examineSeqOthelloKmerQueryCorrectness.cpp

summarize_seqQueryTools_sailfish_results:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarize_seqQueryTools_sailfish_results src/otherProjects/seqQueryTree/summarize_seqQueryTools_sailfish_results.cpp

parseSeqOthelloResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSeqOthelloResults src/otherProjects/seqQueryTree/parseSeqOthelloResults.cpp

parseSbtResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSbtResults src/otherProjects/seqQueryTree/parseSbtResults.cpp

parseSbtAsResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSbtAsResults src/otherProjects/seqQueryTree/parseSbtAsResults.cpp

parseSplitSbtResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSplitSbtResults src/otherProjects/seqQueryTree/parseSplitSbtResults.cpp

mergeMultiSampleSailfishResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeMultiSampleSailfishResults src/otherProjects/seqQueryTree/mergeMultiSampleSailfishResults.cpp

querySpecificSampleTranscriptExpression:
	g++ $(CFLAGS) -o $(OUTPUTDIR)querySpecificSampleTranscriptExpression src/otherProjects/seqQueryTree/querySpecificSampleTranscriptExpression.cpp

bam2fa2kmer:
	g++ $(CFLAGS) -o $(OUTPUTDIR)bam2fa2kmer src/otherProjects/seqQueryTree/bam2fa2kmer.cpp

faVec2KmerVecDir_minCovVarOnFileSize:
	g++ $(CFLAGS) -o $(OUTPUTDIR)faVec2KmerVecDir_minCovVarOnFileSize src/otherProjects/lothello/utils/faVec2KmerVecDir_minCovVarOnFileSize.cpp

generateNewManifest_removeExistData:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateNewManifest_removeExistData src/otherProjects/seqQueryTree/generateNewManifest_removeExistData.cpp

extractReadSeqFromSamWithinRegion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractReadSeqFromSamWithinRegion src/otherProjects/seqQueryTree/extractReadSeqFromSamWithinRegion.cpp

extractExonReferenceSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractExonReferenceSeq src/otherProjects/seqQueryTree/extractExonReferenceSeq.cpp

#imps on 1000 genome
getExonicSNPhomHetDistributionStats:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getExonicSNPhomHetDistributionStats src/otherProjects/incorporateGenomicVariants/utils/getExonicSNPhomHetDistributionStats.cpp

maskSNPposInSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)maskSNPposInSeq src/otherProjects/impsOn1000genome/utils/maskSNPposInSeq.cpp

summarizePsSJstatsList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizePsSJstatsList src/otherProjects/impsOn1000genome/utils/summarizePsSJstatsList.cpp

summarizePsSJcountInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizePsSJcountInPopulation src/otherProjects/impsOn1000genome/utils/summarizePsSJcountInPopulation.cpp

clean_summarizePsSJcountInPopulation:
	rm summarizePsSJcountInPopulation

summarizePsSJexpressionInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizePsSJexpressionInPopulation src/otherProjects/impsOn1000genome/utils/summarizePsSJexpressionInPopulation.cpp	

getData2plotInRforPsSJexpression:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getData2plotInRforPsSJexpression src/otherProjects/impsOn1000genome/utils/getData2plotInRforPsSJexpression.cpp

getData2putInNon2canJuncSupCorTable:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getData2putInNon2canJuncSupCorTable src/otherProjects/impsOn1000genome/utils/getData2putInNon2canJuncSupCorTable.cpp

checkSJflankStringChange:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringChange src/otherProjects/impsOn1000genome/utils/checkSJflankStringChange.cpp

generateFastqDownloadScript:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFastqDownloadScript src/otherProjects/impsOn1000genome/utils/generateFastqDownloadScript.cpp

getHaplotypedSNPforHisat2indexBuilding:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getHaplotypedSNPforHisat2indexBuilding src/otherProjects/impsOn1000genome/utils/getHaplotypedSNPforHisat2indexBuilding.cpp

getRefBaseFreq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getRefBaseFreq src/otherProjects/impsOn1000genome/utils/getRefBaseFreq.cpp

getRefBaseFreqDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getRefBaseFreqDistribution src/otherProjects/impsOn1000genome/utils/getRefBaseFreqDistribution.cpp

getRefBaseFreqDistribution_oddBinNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getRefBaseFreqDistribution_oddBinNum src/otherProjects/impsOn1000genome/utils/getRefBaseFreqDistribution_oddBinNum.cpp

removeNonMainChrMapSamFromHISAT2results:
	g++ $(CFLAGS) -o $(OUTPUTDIR)removeNonMainChrMapSamFromHISAT2results src/otherProjects/impsOn1000genome/utils/removeNonMainChrMapSamFromHISAT2results.cpp

convertBeersSNP2hisat2format:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertBeersSNP2hisat2format src/otherProjects/impsOn1000genome/utils/convertBeersSNP2hisat2format.cpp						

generateRefGstarScript:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateRefGstarScript src/otherProjects/impsOn1000genome/utils/generateRefGstarScript.cpp

generateRefGmps3Script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateRefGmps3Script src/otherProjects/impsOn1000genome/utils/generateRefGmps3Script.cpp	

examineDataValidOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)examineDataValidOrNot src/otherProjects/impsOn1000genome/utils/examineDataValidOrNot.cpp

transferAndUnzipData:
	g++ $(CFLAGS) -o $(OUTPUTDIR)transferAndUnzipData src/otherProjects/impsOn1000genome/utils/transferAndUnzipData.cpp

#map to equine genome cab2 cab3 diff:
getMaxSegLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getMaxSegLength src/otherProjects/equineCab2Cab3/getMaxSegLength.cpp

clean_getMaxSegLength:
	rm getMaxSegLength

extractUnmapPeReadSeqFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractUnmapPeReadSeqFromSam src/otherProjects/equineCab2Cab3/extractUnmapPeReadSeqFromSam.cpp

extractMap2chrUnPeReadSeqFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractMap2chrUnPeReadSeqFromSam src/otherProjects/equineCab2Cab3/extractMap2chrUnPeReadSeqFromSam.cpp

#BEGIN Lothello 

Kmer2sortedKmer:
	g++ $(CFLAGS) -o $(OUTPUTDIR)Kmer2sortedKmer src/otherProjects/lothello/utils/Kmer2sortedKmer.cpp

#ortKmerVec2KmerOccurrenceProfile:

get_geneSeq_transcriptSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)get_geneSeq_transcriptSeq src/otherProjects/lothello/geneCount/buildIndex/get_geneSeq_transcriptSeq.cpp

KmerVec2setSpecificKmerProfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)KmerVec2setSpecificKmerProfile src/otherProjects/lothello/utils/KmerVec2setSpecificKmerProfile.cpp

faVec2KmerVecDir:
	g++ $(CFLAGS) -o $(OUTPUTDIR)faVec2KmerVecDir src/otherProjects/lothello/utils/faVec2KmerVecDir.cpp

diffSRRlist:
	g++ $(CFLAGS) -o $(OUTPUTDIR)diffSRRlist src/otherProjects/lothello/utils/diffSRRlist.cpp

clean_diffSRRlist:
	rm diffSRRlist

buildGenomeSubRegionFaDir:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildGenomeSubRegionFaDir src/otherProjects/lothello/readAssignment2genomeLocalRegion/buildIndex/buildGenomeSubRegionFaDir.cpp

separateRegionAssignedRead:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateRegionAssignedRead src/otherProjects/lothello/readAssignment2genomeLocalRegion/utils/separateRegionAssignedRead.cpp

summaryJfStats:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summaryJfStats src/otherProjects/lothello/metagenomics/others/summaryJfStats.cpp

#END Lothello

#BEGIN querySeq
filterOutSNPmerWithN:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterOutSNPmerWithN src/otherProjects/bloomFilter/querySeq/filterOutSNPmerWithN.cpp

#generateCmdForKmersetInEachNode:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCmdForKmersetInEachNode src/otherProjects/bloomFilter/othello/generateCmdForKmersetInEachNode.cpp
generateCmdForKmersetInEachNode_8set:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCmdForKmersetInEachNode_8set src/otherProjects/bloomFilter/othello/generateCmdForKmersetInEachNode_8set.cpp

assignOthelloNodeKmer2disjointGroup:
	g++ $(CFLAGS) -o $(OUTPUTDIR)assignOthelloNodeKmer2disjointGroup src/otherProjects/bloomFilter/othello/assignOthelloNodeKmer2disjointGroup.cpp

assignOthelloNodeKmer2disjointGroup_varyKmerSetNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)assignOthelloNodeKmer2disjointGroup_varyKmerSetNum src/otherProjects/bloomFilter/othello/assignOthelloNodeKmer2disjointGroup_varyKmerSetNum.cpp

buildLothelloTreeFromScratch:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildLothelloTreeFromScratch src/otherProjects/bloomFilter/othello/buildLothelloTreeFromScratch.cpp

buildLothelloNode:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildLothelloNode src/otherProjects/bloomFilter/othello/buildLothelloNode.cpp

HCtransformKmerCountFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)HCtransformKmerCountFile src/otherProjects/bloomFilter/kmer_utilities/HCtransformKmerCountFile.cpp

reformReadableKmerCountFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformReadableKmerCountFile src/otherProjects/bloomFilter/kmer_utilities/reformReadableKmerCountFile.cpp
#END querySeq

transcriptCount:
	g++	$(CFLAGS) -o $(OUTPUTDIR)transcriptCount src/alignAll_transcriptCount.cpp

geneCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)geneCount src/otherProjects/geneCount/geneCount.cpp

map2annotationOnly:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2annotationOnly src/otherProjects/map2annotationOnly/map2annotationOnly.cpp

compareIntraChromSJwithGTF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareIntraChromSJwithGTF src/MPS3SeqTools/compareIntraChromSJwithGTF.cpp

#xenograft_begin
renameXenograftChrName:
	g++  $(CFLAGS) -o $(OUTPUTDIR)renameXenograftChrName  src/otherProjects/xenograft/renameXenograftChrName.cpp
#xenograft_end

#assign read 2 genomic region BEGIN
checkReadAssignmentAccuracy:
	g++  $(CFLAGS) -o $(OUTPUTDIR)checkReadAssignmentAccuracy src/otherProjects/assignRead2genomeRegion/checkReadAssignmentAccuracy.cpp

beersSam2fa_assignedRegionId:
	g++  $(CFLAGS) -o $(OUTPUTDIR)beersSam2fa_assignedRegionId src/otherProjects/assignRead2genomeRegion/beersSam2fa_assignedRegionId.cpp

#assign read 2 genomic region END

#CircRNA_begin

cirRNA:
	g++  $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)alignAll_cirRNA src/alignAll.cpp

circularRNA_post:
	g++ $(CFLAGS) -o $(OUTPUTDIR)circularRNA_post src/otherProjects/cirRNA/circularRNA_post.cpp

compareBackSpliceResultsWithValidatedList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSpliceResultsWithValidatedList src/otherProjects/cirRNA/utils/compareBackSpliceResultsWithValidatedList.cpp

simulateCirRNAalignment2IGV:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateCirRNAalignment2IGV src/otherProjects/cirRNA/utils/simulateCirRNAalignments2IGV.cpp

extractBackSpliceFromAllJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractBackSpliceFromAllJunc src/otherProjects/cirRNA/extractBackSpliceFromAllJunc.cpp

extractLariatFromBackSpliceSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLariatFromBackSpliceSJ src/otherProjects/cirRNA/utils/extractLariatFromBackSpliceSJ.cpp

extractLariatFromBackSpliceSJ_syntheticSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLariatFromBackSpliceSJ_syntheticSAM src/otherProjects/cirRNA/utils/extractLariatFromBackSpliceSJ_syntheticSAM.cpp

generateCirRNAtranscript:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCirRNAtranscript src/otherProjects/cirRNA/generateCirRNAtranscript_withBackSpliceAndNormalJunction.cpp

extractBackSpliceFromSTARfusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractBackSpliceFromSTARfusionJunc src/otherProjects/cirRNA/utils/extractBackSpliceFromSTARfusionJunc.cpp

reformatMPS2circRNA:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatMPS2circRNA src/otherProjects/cirRNA/utils/reformatMPS2circRNA.cpp

mergeCirRNAreadCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeCirRNAreadCount src/otherProjects/cirRNA/mergeCirRNAreadCount.cpp

#mergeBackSpliceJunctionWithReadCount:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount src/otherProjects/cirRNA/utils/mergeBackSpliceJunctionWithReadCount.cpp

mergeBackSpliceJunctionWithReadCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount src/otherProjects/cirRNA/utils/mergeBackSpliceJunctionWithReadCount_normalize.cpp

reformBackSpliceSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformBackSpliceSAM src/otherProjects/cirRNA/utils/reformBackSpliceSAM.cpp

compareBackSplice2gtfAnnotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSplice2gtfAnnotation src/otherProjects/cirRNA/compareBackSplice2gtfAnnotation.cpp

compareBackSplice2gtfAnnotation_geneName:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareBackSplice2gtfAnnotation_geneName src/otherProjects/cirRNA/compareBackSplice2gtfAnnotation_geneName.cpp

compareCirRNA2gtfAnnotation_withReadCountFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareCirRNA2gtfAnnotation_withReadCountFile src/otherProjects/cirRNA/compareCirRNA2gtfAnnotation_withReadCountFile.cpp

filterBackSpliceJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterBackSpliceJunc src/otherProjects/cirRNA/backSpliceFiltering_anchorSeqSimilarity.cpp

fixHeadTail_cirRNA:
	g++ $(CFLAGS) -D DETECT_CIRCULAR_RNA -o $(OUTPUTDIR)fixHeadTail_cirRNA_withCirRNAtranscriptReadCountFile src/otherMainFiles/alignAll_fixHeadTail.cpp

normalize_foldChange_comparedTo1stSample:
	g++ $(CFLAGS) -o $(OUTPUTDIR)normalize_foldChange_comparedTo1stSample src/otherProjects/cirRNA/normalize_foldChangeCompared1stSample.cpp

countBackSpliceEncompassingReads:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countBackSpliceEncompassingReads src/otherProjects/cirRNA/countBackSpliceEncompassingReads.cpp

#CircRNA_end

#Mark_begin
annotateBackSpliceJuncWithGeneNameAndReadCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)annotateBackSpliceJuncWithGeneNameAndReadCount src/otherProjects/cirRNA/forMarkData/annotateBackSpliceJuncWithGeneNameAndReadCount.cpp

#Mark_end

#Salma_begin
generateCircRNAseq_exon:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCircRNAseq_exon src/otherProjects/cirRNA/forSalmaData/generateCircRNAseq_exon.cpp

generateCircRNAseq_isoform:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCircRNAseq_isoform src/otherProjects/cirRNA/forSalmaData/generateCircRNAseq_isoform.cpp
#Salma_end

initialMappingForFusionDetection:
	g++ $(CFLAGS) -D OUTPUTUNIQUEUNPAIREDBOTHENDSMAPPEDSAM -o $(OUTPUTDIR)alignAll_initialMappingForFusionDetection src/alignAll.cpp

detectFusion_post:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detectFusion_post src/otherProjects/MPS3_fusion/detectFusion_post.cpp

mps3_fusion_post: sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion generateFusionJunc_filter_compare2gtf separateFusionJuncStrandedOrNot remapAndCountEncompRead
	g++ $(CFLAGS) -o $(OUTPUTDIR)mps3_fusion_post src/otherProjects/MPS3_fusion/mps3_fusion_post.cpp

globalMapIncompletePairSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)globalMapIncompletePairSam src/otherProjects/MPS3_fusion/globalMapIncompletePairSam.cpp

globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion src/otherProjects/MPS3_fusion/globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion.cpp

remapAndCountEncompRead:
	g++ $(CFLAGS) -o $(OUTPUTDIR)remapAndCountEncompRead src/otherProjects/MPS3_fusion/remapAndCountEncompRead.cpp

generateFusionJunc_filter_compare2gtf:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJunc_filter_compare2gtf src/otherProjects/MPS3_fusion/generateFusionJunc_filter_compare2gtf.cpp

remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint:
	g++ $(CFLAGS) -o $(OUTPUTDIR)remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint src/otherProjects/MPS3_fusion/remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint.cpp

countUniqueUnpairedReadsEncompassingFusionBreakPoint:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countUniqueUnpairedReadsEncompassingFusionBreakPoint src/otherProjects/MPS3_fusion/countUniqueUnpairedReadsEncompassingFusionBreakPoint.cpp

analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair:
	g++ $(CFLAGS) -o $(OUTPUTDIR)analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair src/otherProjects/MPS3_fusion/analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair.cpp

filterFusionJunc_localIndexMap:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJunc_localIndexMap src/otherProjects/MPS3_fusion/fusionJuncFiltering_localIndexMap.cpp

filterOutBeersGroundTruthSJfromDetectedFusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterOutBeersGroundTruthSJfromDetectedFusionJunc src/otherProjects/MPS3_fusion/utils/filterOutBeersGroundTruthSJfromDetectedFusionJunc.cpp

reformatMPS2AnnotatedFusionResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatMPS2AnnotatedFusionResults src/otherProjects/MPS3_fusion/utils/reformatMPS2AnnotatedFusionResults.cpp

clusterAndGenerateAlterFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)clusterAndGenerateAlterFusion src/otherProjects/MPS3_fusion/clusterAndGenerateAlterFusion.cpp

generateFusionSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionSeq src/otherProjects/MPS3_fusion/utils/generateFusionSeq.cpp

mergeMultiSampleMPS2results:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeMultiSampleMPS2results src/otherProjects/MPS3_fusion/utils/MPS2resultsAnalysis/mergeMultiSampleMPS2results.cpp

fixHeadTail:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fixHeadTail src/otherMainFiles/alignAll_fixHeadTail.cpp

alignAll_PE:
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll src/alignAll_PE.cpp

multi:
	g++  $(CFLAGS) -D CHECK_MULTI -o $(OUTPUTDIR)alignAll_multi src/alignAll.cpp

debug:
	g++  $(CFLAGS) -D DEBUG -o $(OUTPUTDIR)alignAll_debug src/alignAll.cpp	
	
calTime:
	g++  $(CFLAGS) -D CAL_TIME -o $(OUTPUTDIR)alignAll_calTime src/alignAll.cpp	

#annotation:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)alignAll_annotation src/alignAll_transcript.cpp

#annotation_mapInfo:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)alignAll_annotation_mapInfo src/alignAll_transcript.cpp -D MAP_INFO

mapInfo:
	g++  $(CFLAGS) -D MAP_INFO -o $(OUTPUTDIR)alignAll_mapInfo src/alignAll.cpp

noFiltering:
	g++  $(CFLAGS) -D NO_FILTERING -o $(OUTPUTDIR)alignAll_noFiltering src/alignAll.cpp

alignment2Junc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)alignment2Junc src/MPS3SeqTools/alignmentToJunc.cpp

compareJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareJunc src/MPS3SeqTools/compareJunc.cpp

mergeSJfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeSJfile src/MPS3SeqTools/mergeMultiSJfiles2One.cpp

getAlignInferJuncHashFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile src/MPS3SeqTools/getAlignInferJuncHashFromJuncFile.cpp

getAlignInferJuncHashFromSAMfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile.cpp

getAlignInferJuncHashFromJuncFile_doAlterSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile_doAlterSpliceSite src/MPS3SeqTools/getAlignInferJuncHashFromJuncFile_doAlterSpliceSite.cpp

getAlignInferJuncHashFromSAMfile_doAlterSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile_doAlterSpliceSite.cpp

#getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ src/MPS3SeqTools/getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifyJunc.cpp

getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity src/MPS3SeqTools/getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity.cpp

getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ src/MPS3SeqTools/getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ.cpp

getSmallExon:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSmallExon src/MPS3SeqTools/geSmallExonFromSAMfile.cpp

getSoftClipSeqSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSoftClipSeqSAM src/MPS3SeqTools/getSoftClipSeqSAM.cpp
nwdp:
	g++ $(CFLAGS) -o $(OUTPUTDIR)nwdp src/others/nw_DP_new.cpp

nwdp_oneEndOpen:
	g++ $(CFLAGS) -o $(OUTPUTDIR)nwdp_oneEndOpen src/others/nwdp_oneEndOpen.cpp

getShortAnchorSJalignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getShortAnchorSJalignment src/MPS3SeqTools/seperateTrueOrFalseShortAnchorSJalignment.cpp

refineSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM src/MPS3SeqTools/refineAlignment.cpp

refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM src/MPS3SeqTools/refineAlignment_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity.cpp

extendSJspliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extendSJspliceSite src/MPS3SeqTools/extendFromSpliceSiteInsteadOfJunction.cpp

splitSpliceSiteExtensionSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)splitSpliceSiteExtensionSJ src/MPS3SeqTools/spliceSiteExtension2SJ.cpp

compareSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSpliceSite src/MPS3SeqTools/compareSpliceSite.cpp

spliceSiteAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)spliceSiteAnchorSimilarity src/MPS3SeqTools/compare_anchorSimilarity_extension_multiSpliceSite.cpp

detectIncompatibleSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detectIncompatibleSAM src/MPS3SeqTools/extractIncompatiblePairedAlignments.cpp

filterJunc_anchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterJunc_anchorSimilarity src/MPS3SeqTools/classifyJunc_getAlignInferJuncHashFromJuncFile_doAlterSpliceSite.cpp

classifySJ_withAnchorSimilarity:
	g++ $(CFLAGS) -o $(OUTPUTDIR)classifySJ_withAnchorSimilarity src/MPS3SeqTools/classifyJunc_withAlterSpliceSiteAnchorSimilarity.cpp

getPositiveStrandFa_MPS3sam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getPositiveStrandFa_MPS3sam src/MPS3SeqTools/getPositiveStrandFastaFromMPS3sam.cpp

outputSJtypeNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)outputSJtypeNum src/MPS3SeqTools/countSJtypeNum.cpp

countSJcanonicalOrNotNumFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSJcanonicalOrNotNumFromJuncFile src/MPS3SeqTools/countSJcanonicalOrNotNumFromJuncFile.cpp

parseSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseSAM src/MPS3SeqTools/parseSAM.cpp

refineSAM_SE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_SE src/MPS3SeqTools/refineSAM_SE.cpp

refineSAM_PE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_PE src/MPS3SeqTools/refineSAM_PE.cpp

softClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)softClip src/MPS3SeqTools/getSoftClippedStats.cpp

softClip_varReadLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)softClip_varReadLength src/MPS3SeqTools/getSoftClippedStats_varReadLength.cpp

insertColumn2eachLine:
	g++ $(CFLAGS) -o $(OUTPUTDIR)insertColumn2eachLine src/MPS3SeqTools/insertColumn2eachLine.cpp

SJsupportNumDiff:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJsupportNumDiff src/MPS3SeqTools/generateSJsupportNumDifference.cpp

SJsupportNumDiff_targetSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJsupportNumDiff src/MPS3SeqTools/generateSJsupportNumDifference_targetSJ.cpp

getSJsupportNum_targetSJ_differentAligner:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSJsupportNum_targetSJ_differentAligner src/MPS3SeqTools/compareJunctionSupportNumDifference.cpp

SJdiscoveryStats_supportNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)SJdiscoveryStats_supportNum src/MPS3SeqTools/generateSJdiscoveryStats_supportNum.cpp

countAlignmentType:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countAlignmentType src/MPS3SeqTools/countAlignmentTypeNum.cpp

checkSJanchorSimilarity_fromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJanchorSimilarity_fromJuncFile src/MPS3SeqTools/checkSJanchorSimilarity_fromJuncFile.cpp

getBothEndsReadWholeMatchFasta_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_groundTruth src/MPS3SeqTools/getBothEndsReadWholeMatchFasta_groundTruth.cpp

getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth src/MPS3SeqTools/getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth.cpp

getBothEndsReadWithoutSJfasta_groundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBothEndsReadWithoutSJfasta_groundTruth src/MPS3SeqTools/getBothEndsReadWithoutSJfasta_groundTruth.cpp

convertBedtoolsResults2plotInR:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertBedtoolsResults2plotInR src/MPS3SeqTools/convertBedtoolsResults2plotInR.cpp

#sam2alignInferJuncHash_supportNum:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum src/MPS3SeqTools/sam2alignInferJuncHash_supportNum.cpp

#sam2alignInferJuncHash_supportNum_parallel:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_parallel.cpp

#sam2alignInferJuncHash_supportNum_anchorSize:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize.cpp

#sam2alignInferJuncHash_supportNum_anchorSize_parallel:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_parallel.cpp

#sam2alignInferJuncHash_supportNum_anchorSize_XM:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite.cpp

sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile src/MPS3SeqTools/sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile.cpp

compareSam2alignInferJuncHash:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSam2alignInferJuncHash src/MPS3SeqTools/compareSam2alignInferJuncHash.cpp

countSplicedAlignmentNum:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSplicedAlignmentNum src/MPS3SeqTools/countSplicedAlignmentNum.cpp

separateForRevSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateForRevSam src/otherProjects/separateForRevSam.cpp

extractPeAlignmentWithInnerSoftClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractPeAlignmentWithInnerSoftClip src/otherProjects/refineInnerParts_PEreads/extractPeAlignmentWithInnerSoftClip.cpp

refineInnerSoftClip:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineInnerSoftClip src/otherProjects/refineInnerParts_PEreads/refineInnerSoftClip.cpp

exonSize:
	g++ $(CFLAGS) -o $(OUTPUTDIR)exonSize src/others/exonSize.cpp

compareFusionRegion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionRegion src/otherProjects/MPS3_fusion/others/compareFusionRegion.cpp

refineSAM_insertAnnotatedSJ2matchThroughEnd:
	g++ $(CFLAGS) -o $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd src/MPS3SeqTools/refineSAM_insertAnnotatedSJ2matchThroughEnd.cpp

refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly:
	g++ $(CFLAGS) -D COMPLETEONLY_WHEN_INSERT_ANNOTATEDSJ -o $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly src/MPS3SeqTools/refineSAM_insertAnnotatedSJ2matchThroughEnd.cpp

countReadsStartEndMappingPosWithinOverlapRegionList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadsStartEndMappingPosWithinOverlapRegionList src/MPS3SeqTools/countReadsStartEndMappingPosWithinOverlapRegionList.cpp

compareSJfromMultiSamples:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSJfromMultiSamples src/MPS3SeqTools/compareSJfromMultiSamples.cpp

generateUnannotatedSJdistributionAlongRefGenome:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateUnannotatedSJdistributionAlongRefGenome src/MPS3SeqTools/generateUnannotatedSJdistributionAlongRefGenome.cpp

reformJuncFastaFromFusionCancer_oneLineSeqPerFusion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformJuncFastaFromFusionCancer_oneLineSeqPerFusion src/otherProjects/MPS3_fusion/simulate_FusionCancer/reformJuncFastaFromFusionCancer_oneLineSeqPerFusion.cpp

extractUniqueBreakPointFusionFromReformedJuncFastaFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractUniqueBreakPointFusionFromReformedJuncFastaFile src/otherProjects/MPS3_fusion/simulate_FusionCancer/extractUniqueBreakPointFusionFromReformedJuncFastaFile.cpp

simulateFusionReadsFromUniqueBreakPointFusionFasta:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateFusionReadsFromUniqueBreakPointFusionFasta src/otherProjects/MPS3_fusion/simulate_FusionCancer/simulateFusionReadsFromUniqueBreakPointFusionFasta.cpp

analysisFusionBreakPointResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)analysisFusionBreakPointResults src/otherProjects/MPS3_fusion/simulate_FusionCancer/analysisFusionBreakPointResults.cpp

generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile src/otherProjects/MPS3_fusion/utils/generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile.cpp

filterFusionJuncInfoWithSupNumThres:
	g++ $(CFLAGS) -o $(OUTPUTDIR)filterFusionJuncInfoWithSupNumThres src/otherProjects/MPS3_fusion/utils/filterFusionJuncInfoWithSupNumThres.cpp

extractFusionGeneNameFromFusionInfoFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFusionGeneNameFromFusionInfoFile src/otherProjects/MPS3_fusion/utils/extractFusionGeneNameFromFusionInfoFile.cpp

compareTwoGenePairList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareTwoGenePairList src/otherProjects/MPS3_fusion/fusionReviewNAR/compareTwoGenePairList.cpp

compareFusionGeneInfoWithGenePairList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionGeneInfoWithGenePairList src/otherProjects/MPS3_fusion/fusionReviewNAR/compareFusionGeneInfoWithGenePairList.cpp

extractGeneNamePairFromFusionGroundTruthList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractGeneNamePairFromFusionGroundTruthList src/otherProjects/MPS3_fusion/fusionReviewNAR/extractGeneNamePairFromFusionGroundTruthList.cpp

compareStrandedFusionJuncResultsWithGTFannotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareStrandedFusionJuncResultsWithGTFannotation src/otherProjects/MPS3_fusion/utils/compareStrandedFusionJuncResultsWithGTFannotation.cpp

compareFusionDetectionResultsWithFusionMapGroundTruth:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareFusionDetectionResultsWithFusionMapGroundTruth src/otherProjects/MPS3_fusion/utils/compareFusionDetectionResultsWithFusionMapGroundTruth.cpp

countMultiMapSegMaxLength_fusionJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countMultiMapSegMaxLength_fusionJunc src/MPS3SeqTools/countMultiMapSegMaxLength_fusionJunc.cpp

countMultiMapSegMaxLength_regularJunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countMultiMapSegMaxLength_regularJunc src/MPS3SeqTools/countMultiMapSegMaxLength_regularJunc.cpp

compareJunc2groundTruth_outputTotalInfo:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareJunc2groundTruth_outputTotalInfo src/MPS3SeqTools/compareJunc2groundTruth_outputTotalInfo.cpp

extractSplicedAlignmentsFromSAM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSplicedAlignmentsFromSAM src/MPS3SeqTools/extractSplicedAlignmentsFromSAM.cpp

modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount src/MPS3SeqTools/modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount.cpp

separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount src/MPS3SeqTools/separateUnmapPeSamWithinMappedGroup_causedByInvalidJumpCode_notAllowedByHTseqCount.cpp

#BEGIN splice data base
mergeJuncFromMultiFile_readCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeJuncFromMultiFile_readCount src/otherProjects/spliceDataBase/utils/mergeJuncFromMultiFile_readCount.cpp

generateSpliceSiteBranchingRatio:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSpliceSiteBranchingRatio src/otherProjects/spliceDataBase/utils/generateSpliceSiteBranchingRatio.cpp

generateSpliceSiteBranchingRatio_includeContiguousBaseCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSpliceSiteBranchingRatio_includeContiguousBaseCount src/otherProjects/spliceDataBase/utils/generateSpliceSiteBranchingRatio_includeContiguousBaseCount.cpp	

generateSpliceSiteBranchingRatio_includeContiguousBaseCount_includeAnchorSize:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSpliceSiteBranchingRatio_includeContiguousBaseCount_includeAnchorSize src/otherProjects/spliceDataBase/utils/generateSpliceSiteBranchingRatio_includeContiguousBaseCount_includeAnchorSize.cpp	

cmpJuncCovBetweenGroup:
	g++ $(CFLAGS) -o $(OUTPUTDIR)cmpJuncCovBetweenGroup src/otherProjects/spliceDataBase/utils/cmpJuncCovBetweenGroup.cpp

check_altSplice_spliceSiteMutation_inTCGA:
	g++ $(CFLAGS) -o $(OUTPUTDIR)check_altSplice_spliceSiteMutation_inTCGA src/otherProjects/spliceDataBase/utils/check_altSplice_spliceSiteMutation_inTCGA.cpp

parseOutExonFromGeneAnn:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseOutExonFromGeneAnn src/otherProjects/spliceDataBase/utils/parseOutExonFromGeneAnn.cpp

parseOutTranscriptFromGeneAnn:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseOutTranscriptFromGeneAnn src/otherProjects/spliceDataBase/utils/parseOutTranscriptFromGeneAnn.cpp
#END splice data base
#BEGIN splice examiner
intropolis_junc_supSampleReadNumDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)intropolis_junc_supSampleReadNumDistribution src/otherProjects/spliceExaminer/utils/intropolis_junc_supSampleReadNumDistribution.cpp

intropolis_cmp2ann_nonOffset:
	g++ $(CFLAGS) -o $(OUTPUTDIR)intropolis_cmp2ann_nonOffset src/otherProjects/spliceExaminer/utils/intropolis_cmp2ann_nonOffset.cpp

intropolis_cmp2ann_nonOffset_atLeastOneEnd:
	g++ $(CFLAGS) -o $(OUTPUTDIR)intropolis_cmp2ann_nonOffset_atLeastOneEnd src/otherProjects/spliceExaminer/utils/intropolis_cmp2ann_nonOffset_atLeastOneEnd.cpp	

intropolis_extensionSeqSimi:
	g++ $(CFLAGS) -o $(OUTPUTDIR)intropolis_extensionSeqSimi src/otherProjects/spliceExaminer/utils/intropolis_extensionSeqSimi.cpp

spliceExaminer_withSpliceReference:
	g++ $(CFLAGS) -o $(OUTPUTDIR)spliceExaminer_withSpliceReference src/otherProjects/spliceExaminer/spliceExaminer_withSpliceReference.cpp	
#END splice examiner

#BEGIN personal genome alignment
map2syntheticSNPtransSeq_includeGlobalAlignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2syntheticSNPtransSeq_includeGlobalAlignment src/otherProjects/incorporateGenomicVariants/map2syntheticSNPtransSeq_includeGlobalAlignment.cpp

map2syntheticSNPtransSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)map2syntheticSNPtransSeq src/otherProjects/incorporateGenomicVariants/map2syntheticSNPtransSeq.cpp

convertBeersGeneInfo2standardGAF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertBeersGeneInfo2standardGAF src/otherProjects/incorporateGenomicVariants/utils/convertBeersGeneInfo2standardGAF.cpp

extractSJfromGAF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSJfromGAF src/otherProjects/incorporateGenomicVariants/utils/extractSJfromGAF.cpp

reformatLuanSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatLuanSNP src/otherProjects/incorporateGenomicVariants/utils/reformatLuanSNP.cpp

generateSNPhash_fromBeersSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSNPhash_fromBeersSNPfile src/otherProjects/incorporateGenomicVariants/utils/generateSNPhash_fromBeersSNPfile.cpp

simulateGenomReads_SNPpos:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateGenomReads_SNPpos src/otherProjects/incorporateGenomicVariants/simulateGenomReads_SNPpos.cpp

simulateTransReads_SNPpos:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateTransReads_SNPpos src/otherProjects/incorporateGenomicVariants/simulateTransReads_SNPpos.cpp

generateSyntheticSNPtransSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSyntheticSNPtransSeq src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPtransSeq.cpp

generateSyntheticSNPgenomSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSyntheticSNPgenomSeq src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPgenomSeq.cpp

generateDecoyAlignHash_fromSimulatedTransReadAlignResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateDecoyAlignHash_fromSimulatedTransReadAlignResults src/otherProjects/incorporateGenomicVariants/utils/generateDecoyAlignHash_fromSimulatedTransReadAlignResults.cpp

comparePSsj2SDsj:
	g++ $(CFLAGS) -o $(OUTPUTDIR)comparePSsj2SDsj src/otherProjects/incorporateGenomicVariants/utils/comparePSsj2SDsj.cpp

checkSJflankStringChanges:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringChanges src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringChanges.cpp

checkSJflankStringChanges_withFormattedSNPfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringChanges_withFormattedSNPfile src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringChanges_withFormattedSNPfile.cpp

extractSamWithSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithSJ src/otherProjects/incorporateGenomicVariants/utils/extractSamWithSJ.cpp

categorizeSJwithSiteAnnotatedOrNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)categorizeSJwithSiteAnnotatedOrNot src/otherProjects/incorporateGenomicVariants/utils/categorizeSJwithSiteAnnotatedOrNot.cpp

generateCandiPsSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCandiPsSJ src/otherProjects/incorporateGenomicVariants/utils/generateCandiPsSJ.cpp

getBlackListSamFromSimulatedReadSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getBlackListSamFromSimulatedReadSam src/otherProjects/incorporateGenomicVariants/getBlackListSamFromSimulatedReadSam.cpp

getClosestSNPbesideSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getClosestSNPbesideSpliceSite src/otherProjects/incorporateGenomicVariants/utils/getClosestSNPbesideSpliceSite.cpp

checkSJexistencInMultiJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJexistencInMultiJuncFile src/otherProjects/incorporateGenomicVariants/utils/checkSJexistencInMultiJuncFile.cpp

compareTwoAlignInferJuncHash:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareTwoAlignInferJuncHash src/otherProjects/personalGenome/utils/compareTwoAlignInferJuncHash.cpp

compareAlignInferJuncHash2annotatedSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSJ src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSJ.cpp

compareAlignInferJuncHash2annotatedSJ_chrNamePosOnly:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSJ_chrNamePosOnly src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSJ_chrNamePosOnly.cpp	

compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset src/otherProjects/incorporateGenomicVariants/utils/compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset.cpp

compareAlignInferJuncHash2annotatedSpliceSite:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSpliceSite.cpp

compareAlignInferJuncHash2annotatedSpliceSite_checkSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite_checkSNP src/otherProjects/personalGenome/utils/compareAlignInferJuncHash2annotatedSpliceSite_checkSNP.cpp

compareGeneCount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareGeneCount src/otherProjects/incorporateGenomicVariants/utils/compareGeneCount.cpp

checkSNPallelExpression:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSNPallelExpression src/otherProjects/incorporateGenomicVariants/utils/checkSNPallelExpression.cpp;

compareSnpAse_refBase2MPS3sd_altBase2MPS3ps:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSnpAse_refBase2MPS3sd_altBase2MPS3ps src/otherProjects/incorporateGenomicVariants/utils/compareSnpAse_refBase2MPS3sd_altBase2MPS3ps.cpp

checkSNPbesideSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSNPbesideSJ src/otherProjects/incorporateGenomicVariants/utils/checkSNPbesideSJ.cpp

checkPsSJfrequencyInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkPsSJfrequencyInPopulation src/otherProjects/incorporateGenomicVariants/utils/checkPsSJfrequencyInPopulation.cpp

checkTotalSJfrequencyInPopulation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkTotalSJfrequencyInPopulation src/otherProjects/incorporateGenomicVariants/utils/checkTotalSJfrequencyInPopulation.cpp

parseYiPaperPsSJresults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseYiPaperPsSJresults src/otherProjects/incorporateGenomicVariants/yiPaper/parseYiPaperPsSJresults.cpp

parseLuanAlignmentImprovementResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseLuanAlignmentImprovementResults src/otherProjects/incorporateGenomicVariants/utils/parseLuanAlignmentImprovementResults.cpp

labelSJexistenceAgainstAnotherSJlist:
	g++ $(CFLAGS) -o $(OUTPUTDIR)labelSJexistenceAgainstAnotherSJlist src/otherProjects/incorporateGenomicVariants/utils/labelSJexistenceAgainstAnotherSJlist.cpp

checkSJflankStringFrequency:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJflankStringFrequency src/otherProjects/incorporateGenomicVariants/utils/checkSJflankStringFrequency.cpp

checkSJcoverageInSampleWithSNPorNot:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkSJcoverageInSampleWithSNPorNot src/otherProjects/incorporateGenomicVariants/utils/checkSJcoverageInSampleWithSNPorNot.cpp

extractMismatchFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractMismatchFromSam src/otherProjects/incorporateGenomicVariants/integrated2MPS3/extractMismatchFromSam.cpp

generateCandiSNPfromMismatchInSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateCandiSNPfromMismatchInSam src/otherProjects/incorporateGenomicVariants/integrated2MPS3/generateCandiSNPfromMismatchInSam.cpp

checkAlignmentMismatchConsistencyWithMdMp:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkAlignmentMismatchConsistencyWithMdMp src/otherProjects/incorporateGenomicVariants/integrated2MPS3/checkAlignmentMismatchConsistencyWithMdMp.cpp

cmpTwoHaplotypeSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)cmpTwoHaplotypeSNP src/otherProjects/incorporateGenomicVariants/utils/cmpTwoHaplotypeSNP.cpp

extractSamWithSNPinList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithSNPinList src/otherProjects/incorporateGenomicVariants/integrated2MPS3/extractSamWithSNPinList.cpp

parseTCGAmafFIle2individualSomaticSNP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseTCGAmafFIle2individualSomaticSNP src/otherProjects/incorporateGenomicVariants/TCGA/parseTCGAmafFIle2individualSomaticSNP.cpp

detectSomaticMutatedSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)detectSomaticMutatedSJ src/otherProjects/incorporateGenomicVariants/TCGA/detectSomaticMutatedSJ.cpp

extractSamplePairWithSomaticMutation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamplePairWithSomaticMutation src/otherProjects/incorporateGenomicVariants/TCGA/extractSamplePairWithSomaticMutation.cpp
#END personal genome alignment

extractLowSupJuncFromJuncList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractLowSupJuncFromJuncList src/MPS3SeqTools/extractLowSupJuncFromJuncList.cpp

extractSamWithJuncInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithJuncInList src/MPS3SeqTools/extractSamWithJuncInList.cpp

extractSamWithReadNameInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamWithReadNameInList src/otherProjects/incorporateGenomicVariants/utils/extractSamWithReadNameInList.cpp

extractFaWithReadNameInList:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFaWithReadNameInList src/otherProjects/incorporateGenomicVariants/utils/extractFaWithReadNameInList.cpp

countSpliceSiteNumFromJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countSpliceSiteNumFromJuncFile src/MPS3SeqTools/countSpliceSiteNumFromJuncFile.cpp

countReadOverlapWithASM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countReadOverlapWithASM src/otherProjects/manohar/countReadOverlapWithASM.cpp

assignRedCloverGeneName2ASM:
	g++ $(CFLAGS) -o $(OUTPUTDIR)assignRedCloverGeneName2ASM src/otherProjects/manohar/assignRedCloverGeneName2ASM.cpp

#BEGIN MPS3_paper

reformatOtherAlignerReportJunc2mps3format:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatOtherAlignerReportJunc2mps3format src/otherProjects/MPS3_paper/utils/reformatOtherAlignerReportJunc2mps3format.cpp

count_alignmentTypeNum_softClipBaseNum_fixedLength:
	g++ $(CFLAGS) -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_fixedLength.cpp

#count_alignmentTypeNum_softClipBaseNum_varyLength:
#	g++ $(CFLAGS) -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_varyLength.cpp

#count_alignmentTypeNum_softClipBaseNum_fixedLength_includeUnpaired:
#	g++ $(CFLAGS) -D INCLUDE_UNPAIRED -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength_includeUnpaired src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_fixedLength.cpp

#count_alignmentTypeNum_softClipBaseNum_varyLength_includeUnpaired:
#	g++ $(CFLAGS) -D INCLUDE_UNPAIRED -o $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength_includeUnpaired src/otherProjects/MPS3_paper/count_alignmentTypeNum_softClipBaseNum_varyLength.cpp

compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres src/otherProjects/MPS3_paper/compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres.cpp

compareSam2alignInferJuncHash2ann_withVarySupNumThres:
	g++ $(CFLAGS) -o $(OUTPUTDIR)compareSam2alignInferJuncHash2ann_withVarySupNumThres src/otherProjects/MPS3_paper/compareSam2alignInferJuncHash2ann_withVarySupNumThres.cpp

generateSJsupNumCorrelation_sharedByAll:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateSJsupNumCorrelation_sharedByAll src/otherProjects/MPS3_paper/generateSJsupNumCorrelation_sharedByAll.cpp

convertUnmappedBasePercAfterTrimming2beforeTrimming:
	g++ $(CFLAGS) -o $(OUTPUTDIR)convertUnmappedBasePercAfterTrimming2beforeTrimming src/otherProjects/MPS3_paper/utils/convertUnmappedBasePercAfterTrimming2beforeTrimming.cpp

count_mismatchDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)count_mismatchDistribution src/otherProjects/MPS3_paper/count_mismatchDistribution.cpp

produceMappedReadSeqFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceMappedReadSeqFromSam src/otherProjects/MPS3_paper/utils/produceMappedReadSeqFromSam.cpp

produceMappedReadSeqFromSam_PE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceMappedReadSeqFromSam_PE src/otherProjects/MPS3_paper/utils/produceMappedReadSeqFromSam_PE.cpp	

checkInListJuncCoverageInMultiJuncFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkInListJuncCoverageInMultiJuncFile src/otherProjects/MPS3_paper/utils/checkInListJuncCoverageInMultiJuncFile.cpp

reformatSTARjunc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatSTARjunc src/otherProjects/MPS3_paper/utils/reformatSTARjunc.cpp

reformatTophat2junc:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatTophat2junc src/otherProjects/MPS3_paper/utils/reformatTophat2junc.cpp	

checkAlignmentAccuracy:
	g++ $(CFLAGS) -o $(OUTPUTDIR)checkAlignmentAccuracy src/otherProjects/MPS3_paper/checkAlignmentAccuracy.cpp

checkAlignmentAccuracy_perfect:
	g++ $(CFLAGS) -D PERFECT -o $(OUTPUTDIR)checkAlignmentAccuracy_perfect src/otherProjects/MPS3_paper/checkAlignmentAccuracy.cpp

checkAlignmentAccuracy_within:
	g++ $(CFLAGS) -D WITHIN -o $(OUTPUTDIR)checkAlignmentAccuracy_within src/otherProjects/MPS3_paper/checkAlignmentAccuracy.cpp	

separateUniqueMultiSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateUniqueMultiSam src/otherProjects/MPS3_paper/separateUniqueMultiSam.cpp

separateUniqueMultiSam_numericReadName:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateUniqueMultiSam_numericReadName src/otherProjects/MPS3_paper/separateUniqueMultiSam_numericReadName.cpp

countTotalBaseNumInFqFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)countTotalBaseNumInFqFile src/otherProjects/MPS3_paper/utils/countTotalBaseNumInFqFile.cpp
#END MPS3_paper

#START largeNumDataSetProcessCmdLine
mk_sbt_count:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_sbt_count src/otherProjects/largeNumDataSetProcessCmdLine/mk_sbt_count.cpp

mk_runSailfish:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_runSailfish src/otherProjects/largeNumDataSetProcessCmdLine/mk_runSailfish.cpp

mk_mergeTwoEndsFastq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_mergeTwoEndsFastq src/otherProjects/largeNumDataSetProcessCmdLine/mk_mergeTwoEndsFastq.cpp
clean_mk_mergeTwoEndsFastq:
	rm mk_mergeTwoEndsFastq
mk_mergePeFq_keepSRRid:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_mergePeFq_keepSRRid src/otherProjects/largeNumDataSetProcessCmdLine/mk_mergePeFq_keepSRRid.cpp
clean_mk_mergePeFq_keepSRRid:
	rm mk_mergePeFq_keepSRRid
mk_sra2fastq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_sra2fastq src/otherProjects/largeNumDataSetProcessCmdLine/mk_sra2fastq.cpp
clean_mk_sra2fastq:
	rm mk_sra2fastq
mkScript_buildHisat2iSNPindex:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mkScript_buildHisat2iSNPindex src/otherProjects/largeNumDataSetProcessCmdLine/mkScript_buildHisat2iSNPindex.cpp
mkScript_hisat2map_hisat2iSNPmap:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mkScript_hisat2map_hisat2iSNPmap src/otherProjects/largeNumDataSetProcessCmdLine/mkScript_hisat2map_hisat2iSNPmap.cpp	
mk_mps_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_mps_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_mps_cmd_script.cpp
mk_mpsDir_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_mpsDir_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_mpsDir_cmd_script.cpp
mk_checkMPS3completeOrNot_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_checkMPS3completeOrNot_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_checkMPS3completeOrNot_cmd_script.cpp
mk_checkSam2juncCompleteOrNot_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_checkSam2juncCompleteOrNot_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_checkSam2juncCompleteOrNot_cmd_script.cpp	
mk_sam2junc_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_sam2junc_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_sam2junc_cmd_script.cpp
mk_star_cmd_script_SE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_star_cmd_script_SE src/otherProjects/largeNumDataSetProcessCmdLine/mk_star_cmd_script_SE.cpp
mk_star_cmd_script_PE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_star_cmd_script_PE src/otherProjects/largeNumDataSetProcessCmdLine/mk_star_cmd_script_PE.cpp
mk_mkdir_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_mkdir_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_mkdir_cmd_script.cpp
mk_grepMM9sam_cmd_script:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_grepMM9sam_cmd_script src/otherProjects/largeNumDataSetProcessCmdLine/mk_grepMM9sam_cmd_script.cpp
mk_SBTcount:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_SBTcount src/otherProjects/largeNumDataSetProcessCmdLine/mk_SBTcount.cpp
mk_cmpUnmapReadSet:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_cmpUnmapReadSet src/otherProjects/largeNumDataSetProcessCmdLine/mk_cmpUnmapReadSet.cpp
mk_cat2outputSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_cat2outputSam src/otherProjects/largeNumDataSetProcessCmdLine/mk_cat2outputSam.cpp
mk_catMPS3sam2target:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_catMPS3sam2target src/otherProjects/largeNumDataSetProcessCmdLine/mk_catMPS3sam2target.cpp
rm_mk_catMPS3sam2target:
	rm mk_catMPS3sam2target
mk_catFq_fq2fa_rmFq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_catFq_fq2fa_rmFq src/otherProjects/largeNumDataSetProcessCmdLine/mk_catFq_fq2fa_rmFq.cpp
mk_rmUncompressedIndex:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mk_rmUncompressedIndex src/otherProjects/largeNumDataSetProcessCmdLine/mk_rmUncompressedIndex.cpp
#END largeNumDataSetProcessCmdLine

separateMergedChromFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateMergedChromFile src/otherProjects/manohar/separateMergedChromFile.cpp

#START ALU
separateAluFamilyRegion:
	g++ $(CFLAGS) -o $(OUTPUTDIR)separateAluFamilyRegion src/otherProjects/Alu/separateAluFamilyRegion.cpp

generateAluSeqFasta:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateAluSeqFasta src/otherProjects/Alu/generateAluSeqFasta.cpp

sumUpAluReadCountForEachFamily_sailfish:
	g++ $(CFLAGS) -o $(OUTPUTDIR)sumUpAluReadCountForEachFamily_sailfish src/otherProjects/Alu/sumUpAluReadCountForEachFamily_sailfish.cpp

segmentMap2singleAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)segmentMap2singleAluFamily src/otherProjects/Alu/segmentMap2singleAluFamily.cpp

segmentMap2multiAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)segmentMap2multiAluFamily src/otherProjects/Alu/segmentMap2multiAluFamily.cpp

mergeSegMapResultsAndAssignAluFamily:
	g++ $(CFLAGS) -o $(OUTPUTDIR)mergeSegMapResultsAndAssignAluFamily src/otherProjects/Alu/mergeSegMapResultsAndAssignAluFamily.cpp

quantifyAluElement_aluIndexMapSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)quantifyAluElement_aluIndexMapSam src/otherProjects/Alu/quantification/quantifyAluElement_aluIndexMapSam.cpp

quantifyAluElement_genomeMapSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)quantifyAluElement_genomeMapSam src/otherProjects/Alu/quantification/quantifyAluElement_genomeMapSam.cpp	

simulateCompleteReadWithFaFile_SE:
	g++ $(CFLAGS) -o $(OUTPUTDIR)simulateCompleteReadWithFaFile_SE src/otherProjects/simulator/simulateCompleteReadWithFaFile_SE.cpp

queryAluAnn:
	g++ $(CFLAGS) -o $(OUTPUTDIR)queryAluAnn src/otherProjects/Alu/queryAluAnn.cpp

annotateJuncWithAluAnn:
	g++ $(CFLAGS) -o $(OUTPUTDIR)annotateJuncWithAluAnn src/otherProjects/Alu/annotateJuncWithAluAnn.cpp
#END ALU

#BEGIN 16S rRNA

reform16SrRNA_RDP:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reform16SrRNA_RDP src/otherProjects/16SrRNA/utils/reform16SrRNA_RDP.cpp

#END

parseEquineNCBI102AnnToGenerateSJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSJ src/otherProjects/equineAnn/parseEquineNCBI102AnnToGenerateSJ.cpp

addGeneSymbolAccession2SJ:
	g++ $(CFLAGS) -o $(OUTPUTDIR)addGeneSymbolAccession2SJ src/otherProjects/equineAnn/addGeneSymbolAccession2SJ.cpp

getSpliceSiteFromSJ_withGeneSymbolIdAndExonID:
	g++ $(CFLAGS) -o $(OUTPUTDIR)getSpliceSiteFromSJ_withGeneSymbolIdAndExonID src/otherProjects/equineAnn/getSpliceSiteFromSJ_withGeneSymbolIdAndExonID.cpp

#BEGIN extract SJ from gtf
extractSJfromGTF:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSJfromGTF src/otherProjects/extractSJfromGTF/extractSJfromGTF.cpp
#END extract SJ from gtf

reformFeatureCount2deseq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformFeatureCount2deseq src/otherProjects/featureCounts/reformFeatureCount2deseq.cpp

reformSailfish2tpm:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformSailfish2tpm src/otherProjects/sailfish/reformSailfish2tpm.cpp

reformSailfishReadCount2deseq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformSailfishReadCount2deseq src/otherProjects/sailfish/reformSailfishReadCount2deseq.cpp

clean_reformSailfishReadCount2deseq:
	rm reformSailfishReadCount2deseq

buildSNPmerIndex:
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildSNPmerIndex src/otherProjects/indexing/utils/build_snpMerIndex.cpp

parseGeneAnn_generateAnnSJ_generateNovelSJwithAnnEnds:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseGeneAnn_generateAnnSJ_generateNovelSJwithAnnEnds src/otherProjects/parseGeneAnn/parseGeneAnn_generateAnnSJ_generateNovelSJwithAnnEnds.cpp

reformatGeneAnn_addChr:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatGeneAnn_addChr src/otherProjects/reAnn2doFeatureCounts/reformatGeneAnn/reformatGeneAnn_addChr.cpp

reformatGeneAnn_selectChr_1-31_X_Y:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatGeneAnn_selectChr_1-31_X_Y src/otherProjects/reAnn2doFeatureCounts/reformatGeneAnn/reformatGeneAnn_selectChr_1-31_X_Y.cpp

reformatGeneAnn_replaceChrNameWithAccessionIdMap:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatGeneAnn_replaceChrNameWithAccessionIdMap src/otherProjects/reAnn2doFeatureCounts/reformatGeneAnn/reformatGeneAnn_replaceChrNameWithAccessionIdMap.cpp

fastqReadLengthDistribution:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fastqReadLengthDistribution src/otherProjects/10x/fastqReadLengthDistribution.cpp

extractUnmapReadFromSam:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractUnmapReadFromSam src/otherProjects/unmapReadSetCmp/extractUnmapReadFromSam.cpp

fq2fa:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fq2fa src/otherProjects/utilities/fq2fa.cpp

extractChromSubSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractChromSubSeq src/otherProjects/utilities/extractChromSubSeq.cpp

extractFusionSeq_noIns:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFusionSeq_noIns src/otherProjects/utilities/extractFusionSeq_noIns.cpp

extractAnnGeneFusionAnchorSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractAnnGeneFusionAnchorSeq src/otherProjects/utilities/extractAnnGeneFusionAnchorSeq.cpp		

extractSeqOnlyFromFaOrFq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fa2seq src/otherProjects/utilities/extractSeqOnlyFromFaOrFq.cpp

HCtransformation_fa:
	g++ $(CFLAGS) -o $(OUTPUTDIR)HCtransformation_fa src/otherProjects/utilities/HCtransformation_fa.cpp

HCtransformation_genome:
	g++ $(CFLAGS) -o $(OUTPUTDIR)HCtransformation_genome src/otherProjects/utilities/HCtransformation_genome.cpp

produceGenomeSubSeqRead:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceGenomeSubSeqRead src/otherProjects/xenograft/SBTbased/utils/produceGenomeSubSeqRead.cpp

produceGenomeLocalRegionSpecificKmerProfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceGenomeLocalRegionSpecificKmerProfile src/otherProjects/xenograft/SBTbased/utils/produceGenomeLocalRegionSpecificKmerProfile.cpp

produceSetSpecificKmerProfile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)produceSetSpecificKmerProfile src/otherProjects/xenograft/SBTbased/utils/produceSetSpecificKmerProfile.cpp		

summarizeBloomFilterCheckResults:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizeBloomFilterCheckResults src/otherProjects/unmapReadSetCmp/summarizeBloomFilterCheckResults.cpp

#BEGIN Dream Challenge
fusionDetection_dreamChallengeOnly:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fusionDetection_dreamChallengeOnly src/otherProjects/dreamChallenge/fusionDetection_dreamChallengeOnly.cpp

clean_fusionDetection_dreamChallengeOnly:
	rm fusionDetection_dreamChallengeOnly

breakPointSummary_2bedpeFormat:
	g++ $(CFLAGS) -o $(OUTPUTDIR)breakPointSummary_2bedpeFormat src/otherProjects/dreamChallenge/breakPointSummary_2bedpeFormat.cpp

clean_breakPointSummary_2bedpeFormat:
	rm breakPointSummary_2bedpeFormat

parseOneLinePerReadFormatFile2peFqFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseOneLinePerReadFormatFile2peFqFile src/otherProjects/dreamChallenge/utils/parseOneLinePerReadFormatFile2peFqFile.cpp

clean_parseOneLinePerReadFormatFile2peFqFile:
	rm parseOneLinePerReadFormatFile2peFqFile

extractFusionBreakPointAnchorSeq:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractFusionBreakPointAnchorSeq src/otherProjects/dreamChallenge/utils/extractFusionBreakPointAnchorSeq.cpp

clean_extractFusionBreakPointAnchorSeq:
	rm extractFusionBreakPointAnchorSeq

parseOutFeatureTypePosGeneTransIdNameFromGeneAnn:
	g++ $(CFLAGS) -o $(OUTPUTDIR)parseOutFeatureTypePosGeneTransIdNameFromGeneAnn src/otherProjects/dreamChallenge/geneAnn_bedPe/parseOutFeatureTypePosGeneTransIdNameFromGeneAnn.cpp

queryGeneAnnotation:
	g++ $(CFLAGS) -o $(OUTPUTDIR)queryGeneAnnotation src/otherProjects/dreamChallenge/geneAnn_bedPe/queryGeneAnnotation.cpp

queryGeneAnnotation_withIndexFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)queryGeneAnnotation_withIndexFile src/otherProjects/dreamChallenge/geneAnn_bedPe/queryGeneAnnotation_withIndexFile.cpp	

generateGeneAnnotationIndex:
	g++ $(CFLAGS) -o $(OUTPUTDIR)generateGeneAnnotationIndex src/otherProjects/dreamChallenge/geneAnn_bedPe/generateGeneAnnotationIndex.cpp

reformatDreamChallengeBedPeFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)reformatDreamChallengeBedPeFile src/otherProjects/dreamChallenge/geneAnn_bedPe/reformatDreamChallengeBedPeFile.cpp		

cmpMPS3fusionResults2gt:
	g++ $(CFLAGS) -o $(OUTPUTDIR)cmpMPS3fusionResults2gt src/otherProjects/dreamChallenge/utils/cmpMPS3fusionResults2gt.cpp

cmpMPS3genePairResults2gt:
	g++ $(CFLAGS) -o $(OUTPUTDIR)cmpMPS3genePairResults2gt src/otherProjects/dreamChallenge/utils/cmpMPS3genePairResults2gt.cpp	

loadParalogGeneGroup_duplicatedGeneDatabaseFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)loadParalogGeneGroup_duplicatedGeneDatabaseFile src/otherProjects/dreamChallenge/paralogGene/loadParalogGeneGroup_duplicatedGeneDatabaseFile.cpp

loadParalogGeneGroup_reformattedParalogGeneGroupFile:
	g++ $(CFLAGS) -o $(OUTPUTDIR)loadParalogGeneGroup_reformattedParalogGeneGroupFile src/otherProjects/dreamChallenge/paralogGene/loadParalogGeneGroup_reformattedParalogGeneGroupFile.cpp

fixFusion_mainBody_clippedSeg_alignment:
	g++ $(CFLAGS) -o $(OUTPUTDIR)fixFusion_mainBody_clippedSeg_alignment src/otherProjects/dreamChallenge/fixFusion_mainBody_clippedSeg_alignment.cpp

spliceFilter_intropolis_Yi:
	g++ $(CFLAGS) -o $(OUTPUTDIR)spliceFilter_intropolis_Yi src/otherProjects/spliceExaminer/spliceFilter_intropolis_Yi.cpp

summarizeJuncFromMultiSample:
	g++ $(CFLAGS) -o $(OUTPUTDIR)summarizeJuncFromMultiSample src/otherProjects/utilities/summarizeJuncFromMultiSample.cpp

annotatedBreakPoint2genePair:
	g++ $(CFLAGS) -o $(OUTPUTDIR)annotatedBreakPoint2genePair src/otherProjects/dreamChallenge/utils/annotatedBreakPoint2genePair.cpp

#END Dream Challenge

extractSamInSpecificRegions:
	g++ $(CFLAGS) -o $(OUTPUTDIR)extractSamInSpecificRegions src/otherProjects/utilities/extractSamInSpecificRegions.cpp

clean: clean_annotation
	rm $(OUTPUTDIR)alignAll

clean_alignAll_speedUp:
	rm $(OUTPUTDIR)alignAll_speedUp

clean_alignAll:
	rm $(OUTPUTDIR)alignAll

clean_alignAll_phase1only:
	rm $(OUTPUTDIR)alignAll_phase1only:

clean_alignAll_initialMappingForFusionDetection:
	rm $(OUTPUTDIR)alignAll_initialMappingForFusionDetection

clean_detectFusion_post:
	rm $(OUTPUTDIR)detectFusion_post

clean_cirRNA:
	rm $(OUTPUTDIR)alignAll_cirRNA

clean_circularRNA_post:
	rm $(OUTPUTDIR)circularRNA_post

clean_map2annotationOnly:
	rm $(OUTPUTDIR)map2annotationOnly

clean_compareBackSpliceResultsWithValidatedList:
	rm $(OUTPUTDIR)compareBackSpliceResultsWithValidatedList

clean_fixHeadTail:
	rm $(OUTPUTDIR)fixHeadTail

clean_simulateCirRNAalignment2IGV:
	rm $(OUTPUTDIR)simulateCirRNAalignment2IGV

clean_extractBackSpliceFromAllJunc:
	rm $(OUTPUTDIR)extractBackSpliceFromAllJunc

clean_extractLariatFromBackSpliceSJ:
	rm $(OUTPUTDIR)extractLariatFromBackSpliceSJ

clean_extractLariatFromBackSpliceSJ_syntheticSAM:
	rm $(OUTPUTDIR)extractLariatFromBackSpliceSJ_syntheticSAM

clean_generateCirRNAtranscript:
	rm $(OUTPUTDIR)generateCirRNAtranscript

clean_extractBackSpliceFromSTARfusionJunc:
	rm $(OUTPUTDIR)extractBackSpliceFromSTARfusionJunc

clean_reformatMPS2circRNA:
	rm $(OUTPUTDIR)reformatMPS2circRNA

clean_mergeCirRNAreadCount:
	rm $(OUTPUTDIR)mergeCirRNAreadCount

clean_mergeBackSpliceJunctionWithReadCount:
	rm $(OUTPUTDIR)mergeBackSpliceJunctionWithReadCount

clean_reformBackSpliceSAM:
	rm $(OUTPUTDIR)reformBackSpliceSAM

clean_compareBackSplice2gtfAnnotation:
	rm $(OUTPUTDIR)compareBackSplice2gtfAnnotation

clean_compareBackSplice2gtfAnnotation_geneName:
	rm $(OUTPUTDIR)compareBackSplice2gtfAnnotation_geneName

clean_compareCirRNA2gtfAnnotation_withReadCountFile:
	rm $(OUTPUTDIR)compareCirRNA2gtfAnnotation_withReadCountFile

clean_fixHeadTail_cirRNA:
	rm $(OUTPUTDIR)fixHeadTail_cirRNA

clean_multi:
	rm $(OUTPUTDIR)alignAll_multi

clean_debug:
	rm $(OUTPUTDIR)alignAll_debug

clean_calTime:
	rm $(OUTPUTDIR)alignAll_calTime

clean_annotation:
	rm $(OUTPUTDIR)alignAll_annotation

clean_mapInfo:
	rm $(OUTPUTDIR)alignAll_mapInfo

clean_annotation_mapInfo:
	rm $(OUTPUTDIR)alignAll_annotation_mapInfo

clean_arrayList:
	rm $(OUTPUTDIR)alignAll_arrayList

clean_mergeSJfile:
	rm $(OUTPUTDIR)mergeSJfile

clean_getAlignInferJuncHashFromJuncFile:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile

clean_getAlignInferJuncHashFromSAMfile:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile

clean_getAlignInferJuncHashFromJuncFile_doAlterSpliceSite:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncFile_doAlterSpliceSite

clean_getAlignInferJuncHashFromSAMfile_doAlterSpliceSite:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite

clean_getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromSAMfile_doAlterSpliceSite_classifySJ

clean_getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromJuncfile_withAlterSpliceSiteAnchorSimilarity

clean_getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ:
	rm $(OUTPUTDIR)getAlignInferJuncHashFromChrNamePosSupNumOnlyJuncFile_checkSJanchorSimilarity_doAlterSpliceSite_classifySJ

clean_getSmallExon:
	rm $(OUTPUTDIR)getSmallExon

clean_getSoftClipSeqSAM:
	rm $(OUTPUTDIR)getSoftClipSeqSAM

clean_nwdp:
	rm $(OUTPUTDIR)nwdp

clean_nwdp_oneEndOpen:
	rm $(OUTPUTDIR)nwdp_oneEndOpen

clean_getShortAnchorSJalignment:
	rm $(OUTPUTDIR)getShortAnchorSJalignment

clean_refineSAM:
	rm $(OUTPUTDIR)refineSAM

clean_refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)refineSAM_withAlignInferJuncHash_withAlterSpliceSiteAnchorSimilarity

clean_extendSJspliceSite:
	rm $(OUTPUTDIR)extendSJspliceSite

clean_splitSpliceSiteExtensionSJ:
	rm $(OUTPUTDIR)splitSpliceSiteExtensionSJ

clean_compareSpliceSite:
	rm $(OUTPUTDIR)compareSpliceSite

clean_spliceSiteAnchorSimilarity:
	rm $(OUTPUTDIR)spliceSiteAnchorSimilarity

clean_detectIncompatibleSAM:
	rm $(OUTPUTDIR)detectIncompatibleSAM

clean_filterJunc_anchorSimilarity:
	rm $(OUTPUTDIR)filterJunc_anchorSimilarity

clean_getPositiveStrandFa_MPS3sam:
	rm $(OUTPUTDIR)getPositiveStrandFa_MPS3sam

clean_classifySJ_withAnchorSimilarity:
	rm $(OUTPUTDIR)classifySJ_withAnchorSimilarity

clean_outputSJtypeNum:
	rm $(OUTPUTDIR)outputSJtypeNum

clean_parseSAM:
	rm $(OUTPUTDIR)parseSAM

clean_refineSAM_SE:
	rm $(OUTPUTDIR)refineSAM_SE

clean_refineSAM_PE:
	rm $(OUTPUTDIR)refineSAM_PE

clean_softClip:
	rm $(OUTPUTDIR)softClip

clean_softClip_varReadLength:
	rm $(OUTPUTDIR)softClip_varReadLength

clean_insertColumn2eachLine:
	rm $(OUTPUTDIR)insertColumn2eachLine

clean_SJsupportNumDiff:
	rm $(OUTPUTDIR)SJsupportNumDiff

clean_SJsupportNumDiff_targetSJ:
	rm $(OUTPUTDIR)SJsupportNumDiff_targetSJ

clean_getSJsupportNum_targetSJ_differentAligner:
	rm $(OUTPUTDIR)getSJsupportNum_targetSJ_differentAligner

clean_SJdiscoveryStats_supportNum:
	rm $(OUTPUTDIR)SJdiscoveryStats_supportNum

clean_countAlignmentType:
	rm $(OUTPUTDIR)countAlignmentType

clean_checkSJanchorSimilarity_fromJuncFile:
	rm $(OUTPUTDIR)checkSJanchorSimilarity_fromJuncFile

clean_getBothEndsReadWholeMatchFasta_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_groundTruth

clean_getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWholeMatchFasta_forwardDirection_groundTruth

clean_getBothEndsReadWithoutSJfasta_groundTruth:
	rm $(OUTPUTDIR)getBothEndsReadWithoutSJfasta_groundTruth

clean_convertBedtoolsResults2plotInR:
	rm $(OUTPUTDIR)convertBedtoolsResults2plotInR

clean_sam2alignInferJuncHash_supportNum:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum

clean_sam2alignInferJuncHash_supportNum_parallel:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_parallel

clean_sam2alignInferJuncHash_supportNum_anchorSize:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize

clean_sam2alignInferJuncHash_supportNum_anchorSize_parallel:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_parallel

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_withFormattedSNPfile

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepHighSupAlterSite_withFormattedSNPfile

clean_compareSam2alignInferJuncHash:
	rm $(OUTPUTDIR)compareSam2alignInferJuncHash

clean_countSplicedAlignmentNum:
	rm $(OUTPUTDIR)countSplicedAlignmentNum

clean_separateForRevSam:
	rm $(OUTPUTDIR)separateForRevSam

clean_extractPeAlignmentWithInnerSoftClip:
	rm $(OUTPUTDIR)extractPeAlignmentWithInnerSoftClip

clean_transcriptCount:
	rm $(OUTPUTDIR)transcriptCount

clean_refineInnerSoftClip:
	rm $(OUTPUTDIR)refineInnerSoftClip

clean_exonSize:
	rm $(OUTPUTDIR)exonSize

clean_compareFusionRegion:
	rm $(OUTPUTDIR)compareFusionRegion

clean_normalize_foldChange_comparedTo1stSample:
	rm $(OUTPUTDIR)normalize_foldChange_comparedTo1stSample

clean_refineSAM_insertAnnotatedSJ2matchThroughEnd:
	rm $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd

clean_refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly:	
	rm $(OUTPUTDIR)refineSAM_insertAnnotatedSJ2matchThroughEnd_completeOnly

clean_countReadsStartEndMappingPosWithinOverlapRegionList:
	rm $(OUTPUTDIR)countReadsStartEndMappingPosWithinOverlapRegionList

clean_compareSJfromMultiSamples:
	rm $(OUTPUTDIR)compareSJfromMultiSamples

clean_generateUnannotatedSJdistributionAlongRefGenome:
	rm $(OUTPUTDIR)generateUnannotatedSJdistributionAlongRefGenome

clean_globalMapIncompletePairSam:
	rm $(OUTPUTDIR)globalMapIncompletePairSam

clean_globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion:
	rm $(OUTPUTDIR)globalMapOuterSoftClipUniquePairedAlignmentToDetectFusion

clean_remapAndCountEncompRead:
	rm $(OUTPUTDIR)remapAndCountEncompRead

clean_generateFusionJunc_filter_compare2gtf:
	rm $(OUTPUTDIR)generateFusionJunc_filter_compare2gtf

clean_remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint:
	rm $(OUTPUTDIR)remapOuterSoftClipUniquePairedAlignmentAgainstFusionBreakPoint

clean_countUniqueUnpairedReadsEncompassingFusionBreakPoint:
	rm $(OUTPUTDIR)countUniqueUnpairedReadsEncompassingFusionBreakPoint

clean_analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair:
	rm $(OUTPUTDIR)analysisUniqueUnpairedCompleteAlignmentToDetectFusionRegionPair

clean_reformJuncFastaFromFusionCancer_oneLineSeqPerFusion:
	rm $(OUTPUTDIR)reformJuncFastaFromFusionCancer_oneLineSeqPerFusion

clean_extractUniqueBreakPointFusionFromReformedJuncFastaFile:
	rm $(OUTPUTDIR)extractUniqueBreakPointFusionFromReformedJuncFastaFile

clean_simulateFusionReadsFromUniqueBreakPointFusionFasta:
	rm $(OUTPUTDIR)simulateFusionReadsFromUniqueBreakPointFusionFasta

clean_analysisFusionBreakPointResults:
	rm $(OUTPUTDIR)analysisFusionBreakPointResults

clean_generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile:
	rm $(OUTPUTDIR)generateFusionJuncInfoFromFusionReadWithBreakPointInfoFile

clean_generateFusionJuncInfoFromAdjustedBreakPointFile:
	rm $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile

clean_compareFusionJuncResultsWithGTFannotation:
	rm $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation

clean_filterFusionJuncInfoWithSupNumThres:
	rm $(OUTPUTDIR)filterFusionJuncInfoWithSupNumThres

clean_extractFusionGeneNameFromFusionInfoFile:
	rm $(OUTPUTDIR)extractFusionGeneNameFromFusionInfoFile

clean_compareTwoGenePairList:
	rm $(OUTPUTDIR)compareTwoGenePairList

clean_compareFusionGeneInfoWithGenePairList:
	rm $(OUTPUTDIR)compareFusionGeneInfoWithGenePairList

clean_compareStrandedFusionJuncResultsWithGTFannotation:
	rm $(OUTPUTDIR)compareStrandedFusionJuncResultsWithGTFannotation

clean_clusterFusionJunc:
	rm $(OUTPUTDIR)clusterFusionJunc

clean_separateFusionJuncStrandedOrNot:
	rm $(OUTPUTDIR)separateFusionJuncStrandedOrNot

clean_filterFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity

clean_filterFusionJunc_localIndexMap:
	rm $(OUTPUTDIR)filterFusionJunc_localIndexMap

clean_filterNonStrandedFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity

clean_filterOutBeersGroundTruthSJfromDetectedFusionJunc:
	rm $(OUTPUTDIR)filterOutBeersGroundTruthSJfromDetectedFusionJunc

clean_filterBackSpliceJunc:
	rm $(OUTPUTDIR)filterBackSpliceJunc

clean_clusterAndGenerateAlterFusion:
	rm $(OUTPUTDIR)clusterAndGenerateAlterFusion

clean_compareIntraChromSJwithGTF:
	rm $(OUTPUTDIR)compareIntraChromSJwithGTF

clean_reformatMPS2AnnotatedFusionResults:
	rm $(OUTPUTDIR)reformatMPS2AnnotatedFusionResults

clean_generateFusionSeq:
	rm $(OUTPUTDIR)generateFusionSeq

clean_countMultiMapSegMaxLength_fusionJunc:
	rm $(OUTPUTDIR)countMultiMapSegMaxLength_fusionJunc

clean_countMultiMapSegMaxLength_regularJunc:
	rm $(OUTPUTDIR)countMultiMapSegMaxLength_regularJunc

clean_compareJunc2groundTruth_outputTotalInfo:
	rm $(OUTPUTDIR)compareJunc2groundTruth_outputTotalInfo

clean_geneCount:
	rm $(OUTPUTDIR)geneCount

clean_compareFusionDetectionResultsWithFusionMapGroundTruth:
	rm $(OUTPUTDIR)compareFusionDetectionResultsWithFusionMapGroundTruth

clean_extractSplicedAlignmentsFromSAM:
	rm $(OUTPUTDIR)extractSplicedAlignmentsFromSAM

clean_cluster_filterAlterFusionJunc_anchorSeqSimilarity:
	rm $(OUTPUTDIR)cluster_filterAlterFusionJunc_anchorSeqSimilarity

clean_modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount:
	rm $(OUTPUTDIR)modifyReadNameToMakePEreadsTheSame_requiredByHTseqCount

clean_reformatBeersSNP:
	rm $(OUTPUTDIR)reformatBeersSNP

clean_generateSNPhash_fromBeersSNPfile:
	rm $(OUTPUTDIR)generateSNPhash_fromBeersSNPfile

clean_simulateGenomReads_SNPpos:
	rm $(OUTPUTDIR)simulateGenomReads_SNPpos

clean_simulateTransReads_SNPpos:
	rm $(OUTPUTDIR)simulateTransReads_SNPpos

clean_generateDecoyAlignHash_fromSimulatedTransReadAlignResults:
	rm $(OUTPUTDIR)generateDecoyAlignHash_fromSimulatedTransReadAlignResults

clean_comparePSsj2SDsj:
	rm $(OUTPUTDIR)comparePSsj2SDsj

clean_checkSJflankStringChanges:
	rm $(OUTPUTDIR)checkSJflankStringChanges

clean_extractSamWithSJ:
	rm $(OUTPUTDIR)extractSamWithSJ

clean_getBlackListSamFromSimulatedReadSam:
	rm $(OUTPUTDIR)getBlackListSamFromSimulatedReadSam

clean_getClosestSNPbesideSpliceSite:
	rm $(OUTPUTDIR)getClosestSNPbesideSpliceSite

clean_extractSamWithJuncInList:
	rm $(OUTPUTDIR)extractSamWithJuncInList

clean_countSJcanonicalOrNotNumFromJuncFile:
	rm $(OUTPUTDIR)countSJcanonicalOrNotNumFromJuncFile

clean_countReadOverlapWithASM:
	rm $(OUTPUTDIR)countReadOverlapWithASM

clean_assignRedCloverGeneName2ASM:
	rm $(OUTPUTDIR)assignRedCloverGeneName2ASM

clean_count_alignmentTypeNum_softClipBaseNum_fixedLength:
	rm $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_fixedLength

clean_count_alignmentTypeNum_softClipBaseNum_varyLength:
	rm $(OUTPUTDIR)count_alignmentTypeNum_softClipBaseNum_varyLength

clean_compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres:
	rm $(OUTPUTDIR)compareSam2alignInferJuncHash2groundTruth_withVarySupNumThres

clean_generateSJsupNumCorrelation_sharedByAll:
	rm $(OUTPUTDIR)generateSJsupNumCorrelation_sharedByAll

clean_separateMergedChromFile:
	rm $(OUTPUTDIR)separateMergedChromFile

clean_separateAluFamilyRegion:
	rm $(OUTPUTDIR)separateAluFamilyRegion

clean_generateAluSeqFasta:
	rm $(OUTPUTDIR)generateAluSeqFasta

clean_sumUpAluReadCountForEachFamily_sailfish:
	rm $(OUTPUTDIR)sumUpAluReadCountForEachFamily_sailfish

clean_segmentMap2singleAluFamily:
	rm $(OUTPUTDIR)segmentMap2singleAluFamily

clean_parseEquineNCBI102AnnToGenerateSJ:
	rm $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSJ

clean_parseEquineNCBI102AnnToGenerateSpliceSite:
	rm $(OUTPUTDIR)parseEquineNCBI102AnnToGenerateSpliceSite

clean_extractSJfromGTF:
	rm $(OUTPUTDIR)extractSJfromGTF

clean_reformFeatureCount2deseq:
	rm $(OUTPUTDIR)reformFeatureCount2deseq

clean_reformSailfish2tpm:
	rm $(OUTPUTDIR)reformSailfish2tpm

clean_generateCandiPsSJ:
	rm $(OUTPUTDIR)generateCandiPsSJ

clean_addGeneSymbolAccession2SJ:
	rm $(OUTPUTDIR)addGeneSymbolAccession2SJ

clean_getSpliceSiteFromSJ_withGeneSymbolIdAndExonID:
	rm $(OUTPUTDIR)getSpliceSiteFromSJ_withGeneSymbolIdAndExonID

clean_extractSpecificAluFamilySeq:
	rm $(OUTPUTDIR)extractSpecificAluFamilySeq

clean_segmentMap2multiAluFamily:
	rm $(OUTPUTDIR)segmentMap2multiAluFamily

clean_mergeSegMapResultsAndAssignAluFamily:
	rm $(OUTPUTDIR)mergeSegMapResultsAndAssignAluFamily

clean_alignAll_outputUnpaired:
	rm $(OUTPUTDIR)alignAll_outputUnpaired

clean_convertUnmappedBasePercAfterTrimming2beforeTrimming:
	rm $(OUTPUTDIR)convertUnmappedBasePercAfterTrimming2beforeTrimming

clean_count_mismatchDistribution:
	rm $(OUTPUTDIR)count_mismatchDistribution

clean_produceMappedReadSeqFromSam:
	rm $(OUTPUTDIR)produceMappedReadSeqFromSam

clean_compareTwoAlignInferJuncHash:
	rm $(OUTPUTDIR)compareTwoAlignInferJuncHash

clean_compareAlignInferJuncHash2annotatedSJ:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSJ

clean_compareAlignInferJuncHash2annotatedSpliceSite:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite

clean_compareAlignInferJuncHash2annotatedSpliceSite_checkSNP:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2annotatedSpliceSite_checkSNP

clean_generateSyntheticSNPtransSeq:
	rm $(OUTPUTDIR)generateSyntheticSNPtransSeq

clean_alignAll_personalizedSeq_SNPseqIndex:
	rm $(OUTPUTDIR)alignAll_personalizedSeq_SNPseqIndex

clean_convertBeersGeneInfo2standardGAF:
	rm $(OUTPUTDIR)convertBeersGeneInfo2standardGAF

clean_map2syntheticSNPtransSeq:
	rm $(OUTPUTDIR)map2syntheticSNPtransSeq

clean_map2syntheticSNPtransSeq_includeGlobalAlignment:
	rm $(OUTPUTDIR)map2syntheticSNPtransSeq_includeGlobalAlignment

clean_extractSJfromGAF:
	rm $(OUTPUTDIR)extractSJfromGAF

clean_sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ:
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify_keepAnnotatedSJ

clean_checkSNPallelExpression:
	rm $(OUTPUTDIR)checkSNPallelExpression

clean_compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset:
	rm $(OUTPUTDIR)compareAlignInferJuncHash2another_chrNamePosOnly_nonOffset

clean_compareSnpAse_refBase2MPS3sd_altBase2MPS3ps:
	rm $(OUTPUTDIR)compareSnpAse_refBase2MPS3sd_altBase2MPS3ps

clean_checkSNPbesideSJ:
	rm $(OUTPUTDIR)checkSNPbesideSJ

clean_checkPsSJfrequencyInPopulation:
	rm $(OUTPUTDIR)checkPsSJfrequencyInPopulation

clean_checkTotalSJfrequencyInPopulation:
	rm $(OUTPUTDIR)checkTotalSJfrequencyInPopulation

clean_checkSJflankStringChanges_withFormattedSNPfile:
	rm $(OUTPUTDIR)checkSJflankStringChanges_withFormattedSNPfile

clean_parseYiPaperPsSJresults:
	rm $(OUTPUTDIR)parseYiPaperPsSJresults

clean_parseLuanAlignmentImprovementResults:
	rm $(OUTPUTDIR)parseLuanAlignmentImprovementResults

clean_labelSJexistenceAgainstAnotherSJlist:
	rm $(OUTPUTDIR)labelSJexistenceAgainstAnotherSJlist

clean_checkSJcoverageInSampleWithSNPorNot:
	rm $(OUTPUTDIR)checkSJcoverageInSampleWithSNPorNot

clean_generateCircRNAseq:
	rm $(OUTPUTDIR)generateCircRNAseq

clean_extractMismatchFromSam:
	rm $(OUTPUTDIR)extractMismatchFromSam

clean_generateCandiSNPfromMismatchInSam:
	rm $(OUTPUTDIR)generateCandiSNPfromMismatchInSam

clean_cmpTwoHaplotypeSNP:
	rm $(OUTPUTDIR)cmpTwoHaplotypeSNP

clean_MapSplice:
	rm $(OUTPUTDIR)MapSplice

clean_iMapSplice:
	rm $(OUTPUTDIR)iMapSplice

clean_extractSamWithSNPinList:
	rm $(OUTPUTDIR)extractSamWithSNPinList

clean_reform16SrRNA_RDP:
	rm $(OUTPUTDIR)reform16SrRNA_RDP

clean_mpsI_arrayList:
	rm $(OUTPUTDIR)mpsI_arrayList

clean_mergeJuncFromMultiFile_readCount:
	rm $(OUTPUTDIR)mergeJuncFromMultiFile_readCount

clean_buildSNPmerIndex:
	rm $(OUTPUTDIR)buildSNPmerIndex

clean_buildSNPmerIndex_embeded:
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded


clean_mps:
	rm $(OUTPUTDIR)mps
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm $(OUTPUTDIR)reportJunc_embeded
	rm $(OUTPUTDIR)removeFalseJuncSam

clean_mps_arrayList:
	rm $(OUTPUTDIR)mps_arrayList
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm $(OUTPUTDIR)reportJunc_embeded
	rm $(OUTPUTDIR)removeFalseJuncSam

clean_mps_batch:
	rm $(OUTPUTDIR)mps_batch
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_batch
	rm $(OUTPUTDIR)reportJunc_embeded

clean_mps_backSpliceIncluded:
	rm $(OUTPUTDIR)mps_backSpliceIncluded
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm $(OUTPUTDIR)reportJunc_embeded

clean_mps_arrayList_backSpliceIncluded:
	rm $(OUTPUTDIR)mps_arrayList_backSpliceIncluded
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_nonBatch
	rm $(OUTPUTDIR)reportJunc_embeded

clean_mps_batch_backSpliceIncluded:
	rm $(OUTPUTDIR)mps_batch_backSpliceIncluded
	rm $(OUTPUTDIR)buildSNPmerIndex_embeded_mps_batch
	rm $(OUTPUTDIR)reportJunc_embeded

clean_checkInListJuncCoverageInMultiJuncFile:
	rm checkInListJuncCoverageInMultiJuncFile

clean_reformatSTARjunc:
	rm reformatSTARjunc

clean_reformatTophat2junc:
	rm reformatTophat2junc	

clean_intropolis_junc_supSampleReadNumDistribution:
	rm intropolis_junc_supSampleReadNumDistribution

clean_intropolis_cmp2ann_nonOffset:
	rm intropolis_cmp2ann_nonOffset

clean_intropolis_cmp2ann_nonOffset_atLeastOneEnd:
	rm intropolis_cmp2ann_nonOffset_atLeastOneEnd

clean_intropolis_extensionSeqSimi:
	rm intropolis_extensionSeqSimi

clean_compareAlignInferJuncHash2annotatedSJ_chrNamePosOnly:
	rm compareAlignInferJuncHash2annotatedSJ_chrNamePosOnly

clean_parseGeneAnn_generateAnnSJ_generateNovelSJwithAnnEnds:
	rm parseGeneAnn_generateAnnSJ_generateNovelSJwithAnnEnds

clean_generateSpliceSiteBranchingRatio:
	rm generateSpliceSiteBranchingRatio

clean_generateSpliceSiteBranchingRatio_includeContiguousBaseCount:
	rm generateSpliceSiteBranchingRatio_includeContiguousBaseCount

clean_generateSpliceSiteBranchingRatio_includeContiguousBaseCount_includeAnchorSize:
	rm generateSpliceSiteBranchingRatio_includeContiguousBaseCount_includeAnchorSize

clean_check_altSplice_spliceSiteMutation_inTCGA:
	rm check_altSplice_spliceSiteMutation_inTCGA

clean_parseOutExonFromGeneAnn:
	rm parseOutExonFromGeneAnn

clean_reformatGeneAnn_addChr:
	rm reformatGeneAnn_addChr

clean_reformatGeneAnn_selectChr_1-31_X_Y:
	rm reformatGeneAnn_selectChr_1-31_X_Y

clean_reformatGeneAnn_replaceChrNameWithAccessionIdMap:
	rm reformatGeneAnn_replaceChrNameWithAccessionIdMap

clean_fastqReadLengthDistribution:
	rm fastqReadLengthDistribution

clean_spliceExaminer_withSpliceReference:
	rm spliceExaminer_withSpliceReference

clean_quantifyAluElement_aluIndexMapSam:
	rm quantifyAluElement_aluIndexMapSam

clean_simulateCompleteReadWithFaFile_SE:
	rm simulateCompleteReadWithFaFile_SE

clean_mergeMultiSampleMPS2results:
	rm mergeMultiSampleMPS2results

clean_mps3_fusion_post:
	rm mps3_fusion_post

clean_cmpJuncCovBetweenGroup:
	rm cmpJuncCovBetweenGroup

clean_checkAlignmentAccuracy:
	rm checkAlignmentAccuracy

clean_checkAlignmentAccuracy_perfect:
	rm checkAlignmentAccuracy_perfect

clean_checkAlignmentAccuracy_within:
	rm checkAlignmentAccuracy_within

clean_compareSam2alignInferJuncHash2ann_withVarySupNumThres:
	rm compareSam2alignInferJuncHash2ann_withVarySupNumThres

clean_renameXenograftChrName:
	rm renameXenograftChrName

clean_separateUniqueMultiSam:
	rm separateUniqueMultiSam

clean_separateUniqueMultiSam_numericReadName:
	rm separateUniqueMultiSam_numericReadName

clean_mk_mps_cmd_script:
	rm mk_mps_cmd_script

clean_mk_mpsDir_cmd_script:
	rm mk_mpsDir_cmd_script

clean_mk_checkMPS3completeOrNot_cmd_script:
	rm mk_checkMPS3completeOrNot_cmd_script

clean_countTotalBaseNumInFqFile:
	rm countTotalBaseNumInFqFile

clean_mk_sam2junc_cmd_script:
	rm mk_sam2junc_cmd_script

clean_mk_star_cmd_script_SE:
	rm mk_star_cmd_script_SE

clean_mk_star_cmd_script_PE:
	rm mk_star_cmd_script_PE

clean_mk_mkdir_cmd_script:
	rm mk_mkdir_cmd_script

clean_mk_grepMM9sam_cmd_script:
	rm mk_grepMM9sam_cmd_script

clean_mk_SBTcount:
	rm mk_SBTcount
	
clean_mk_cmpUnmapReadSet:
	rm mk_cmpUnmapReadSet

clean_mk_cat2outputSam:
	rm mk_cat2outputSam	

clean_parseOutTranscriptFromGeneAnn:
	rm parseOutTranscriptFromGeneAnn

clean_parseTCGAmafFIle2individualSomaticSNP:
	rm parseTCGAmafFIle2individualSomaticSNP

clean_detectSomaticMutatedSJ:
	rm detectSomaticMutatedSJ

clean_extractUnmapReadFromSam:
	rm extractUnmapReadFromSam

clean_fq2fa:
	rm fq2fa

clean_extractSeqOnlyFromFaOrFq:
	rm extractSeqOnlyFromFaOrFq

clean_HCtransformation_fa:
	rm HCtransformation_fa

clean_HCtransformation_genome:
	rm HCtransformation_genome

clean_produceGenomeSubSeqRead:
	rm produceGenomeSubSeqRead

clean_summarizeBloomFilterCheckResults:
	rm summarizeBloomFilterCheckResults

clean_parseOutFeatureTypePosGeneTransIdNameFromGeneAnn:
	rm parseOutFeatureTypePosGeneTransIdNameFromGeneAnn

clean_queryGeneAnnotation:
	rm queryGeneAnnotation

clean_queryGeneAnnotation_withIndexFile:
	rm queryGeneAnnotation_withIndexFile	

clean_generateGeneAnnotationIndex:
	rm generateGeneAnnotationIndex

clean_reformatDreamChallengeBedPeFile:
	rm reformatDreamChallengeBedPeFile

clean_extractChromSubSeq:
	rm extractChromSubSeq

clean_reportJunc_embeded:
	rm reportJunc_embeded

clean_removeFalseJuncSam:
	rm removeFalseJuncSam

clean_loadParalogGeneGroup_duplicatedGeneDatabaseFile:
	rm loadParalogGeneGroup_duplicatedGeneDatabaseFile

clean_loadParalogGeneGroup_reformattedParalogGeneGroupFile:
	rm loadParalogGeneGroup_reformattedParalogGeneGroupFile

clean_extractSamplePairWithSomaticMutation:
	rm extractSamplePairWithSomaticMutation

clean_mk_catFq_fq2fa_rmFq:
	rm mk_catFq_fq2fa_rmFq

clean_extractFusionSeq_noIns:
	rm extractFusionSeq_noIns

clean_queryAluAnn:
	rm queryAluAnn

clean_annotateJuncWithAluAnn:
	rm annotateJuncWithAluAnn

clean_fixFusion_mainBody_clippedSeg_alignment:
	rm fixFusion_mainBody_clippedSeg_alignment

clean_filterOutSNPmerWithN:
	rm filterOutSNPmerWithN

clean_mps3_fusion_post_dreamChallenge:
	rm $(OUTPUTDIR)mps3_fusion_post_dreamChallenge
	rm $(OUTPUTDIR)sam2alignInferJuncHash_supportNum_anchorSize_XM_parallel_classify 
	rm $(OUTPUTDIR)detection_globalMap 
	rm $(OUTPUTDIR)generateFusionJuncInfoFromAdjustedBreakPointFile 
	rm $(OUTPUTDIR)separateFusionJuncStrandedOrNot 
	rm $(OUTPUTDIR)filterFusionJunc_anchorSeqSimilarity 
	rm $(OUTPUTDIR)filterNonStrandedFusionJunc_anchorSeqSimilarity 
	rm $(OUTPUTDIR)compareFusionJuncResultsWithGTFannotation
	rm $(OUTPUTDIR)filtering_supNum_repeat

clean_spliceFilter_intropolis_Yi:
	rm spliceFilter_intropolis_Yi

clean_summarizeJuncFromMultiSample:
	rm summarizeJuncFromMultiSample

clean_extractAnnGeneFusionAnchorSeq:
	rm extractAnnGeneFusionAnchorSeq

clean_annotatedBreakPoint2genePair:
	rm annotatedBreakPoint2genePair

clean_generateCmdForKmersetInEachNode:
	rm generateCmdForKmersetInEachNode

clean_filtering_supNum_repeat:
	rm filtering_supNum_repeat

clean_cmpMPS3genePairResults2gt:
	rm cmpMPS3genePairResults2gt

clean_annotateBackSpliceJuncWithGeneNameAndReadCount:
	rm annotateBackSpliceJuncWithGeneNameAndReadCount

clean_assignOthelloNodeKmer2disjointGroup:
	rm assignOthelloNodeKmer2disjointGroup

clean_generateCmdForKmersetInEachNode_8set:
	rm generateCmdForKmersetInEachNode_8set

clean_mk_rmUncompressedIndex:
	rm mk_rmUncompressedIndex

clean_convertMPS3fusionResults2bedpe:
	rm convertMPS3fusionResults2bedpe

clean_extractSamInSpecificRegions:
	rm extractSamInSpecificRegions

clean_HCtransformKmerCountFile:
	rm HCtransformKmerCountFile

clean_reformReadableKmerCountFile:
	rm reformReadableKmerCountFile

clean_assignOthelloNodeKmer2disjointGroup_varyKmerSetNum:
	rm assignOthelloNodeKmer2disjointGroup_varyKmerSetNum

clean_buildLothelloTreeFromScratch:
	rm buildLothelloTreeFromScratch

clean_buildLothelloNode:
	rm buildLothelloNode

clean_produceGenomeLocalRegionSpecificKmerProfile:
	rm produceGenomeLocalRegionSpecificKmerProfile

clean_produceSetSpecificKmerProfile:
	rm produceSetSpecificKmerProfile

clean_Kmer2sortedKmer:
	rm Kmer2sortedKmer

clean_get_geneSeq_transcriptSeq:
	rm get_geneSeq_transcriptSeq

clean_KmerVec2setSpecificKmerProfile:
	rm KmerVec2setSpecificKmerProfile

clean_faVec2KmerVecDir:
	rm faVec2KmerVecDir

clean_buildGenomeSubRegionFaDir:
	rm buildGenomeSubRegionFaDir

clean_summaryJfStats:
	rm summaryJfStats

clean_examineDataValidOrNot:
	rm examineDataValidOrNot

clean_transferAndUnzipData:
	rm transferAndUnzipData

clean_extractReadSeqFromSamWithinRegion:
	rm extractReadSeqFromSamWithinRegion

clean_extractUnmapPeReadSeqFromSam:
	rm extractUnmapPeReadSeqFromSam

clean_extractMap2chrUnPeReadSeqFromSam:
	rm extractMap2chrUnPeReadSeqFromSam

clean_convertBeersSNP2hisat2format:
	rm convertBeersSNP2hisat2format

clean_maskSNPposInSeq:
	rm maskSNPposInSeq

clean_checkSJflankStringChange:
	rm checkSJflankStringChanges

clean_generateFastqDownloadScript:
	rm generateFastqDownloadScript

clean_getHaplotypedSNPforHisat2indexBuilding:
	rm getHaplotypedSNPforHisat2indexBuilding

clean_getRefBaseFreq:
	rm getRefBaseFreq

clean_getRefBaseFreqDistribution:
	rm getRefBaseFreqDistribution

clean_removeNonMainChrMapSamFromHISAT2results:
	rm removeNonMainChrMapSamFromHISAT2results

clean_generateRefGstarScript:
	rm generateRefGstarScript

clean_generateRefGmps3Script:
	rm generateRefGmps3Script
