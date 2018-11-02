// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef OPTION_INFO_H
#define OPTION_INFO_H

#include <string.h>
#include <string>

#include <getopt.h>

using namespace std;

class Option_Info
{
public:
	string read_file_path_1;
	//bool file_1_fasta_or_fastq_bool;
	string read_file_path_2;
	//bool file_2_fasta_or_fastq_bool;
	string read_file_path_SE;

	int threads_num;
	string global_index_file_path_prefix;
	string local_index_file_path_prefix;
	//string chromsome_file_path_prefix;
	bool fasta_or_fastq_bool;
	string outputFolder_path;
	bool Do_phase1_only_bool;

	bool SE_or_PE_bool;

	int mappedLength_perc_min;

	int mappedLength_base_min;

	string annotation_file_path;
	bool annotation_provided_bool;

	string spliceJunctionAlignInferHash_file_path;
	bool spliceJunctionAlignInferHash_provided_bool;

	string realGenome_file_path;
	bool realGenome_provided_bool;

	string transcript_file_path;
	bool transcript_provided_bool;

	string transcript_type;
	bool transcript_type_assigned_bool;
	string optionStr;

	string inputSamFilePath;
	string alreadyMappedAlignmentFile;
	bool extractUnmapAlignment2ReadFile_bool;

	string SNPfilePath;
	bool SNP_provided_bool;

	string SNP_seq_index_path;
	bool SNP_seq_index_provided_bool;

	bool segMap2SNPmer_phase_defined_bool;
	bool segMap2SNPmer_bothPhase_or_phase2Only_bool;

	bool sharedThreadForIO_bool;

	bool readFileIn_listFile_or_commandLine_bool;

	bool sharing_splicing_context_bool;
	bool sharing_SNP_context_bool;

	bool stopLearningSNP_bool;
	bool stopLearningSplicing_bool;

	bool updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool;
	bool updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool;

	bool backSplice_search_bool;
	bool fusion_search_bool;

	bool avoid_creating_folder_bool;

	bool report_sam_only_bool;

	bool keep_tmp_bool;
	
	int fusion_post_sup_min;
	
	string fusion_post_formatted_gtf_path;
	bool fusion_post_formatted_gtf_provided_bool;
	
	string fusion_post_paralog_gene_path;
	bool fusion_post_paralog_gene_provided_bool;
	
	Option_Info()
	{
		optionStr = "W:1:2:S:G:L:T:F:O:M:I:Z:P:Q:R:U:Y:D:C:7:8:9:ABEXNJKH3456";
		SE_or_PE_bool = false;
		Do_phase1_only_bool = false;
		annotation_provided_bool = false;
		spliceJunctionAlignInferHash_provided_bool = false;
		realGenome_provided_bool = false;
		transcript_provided_bool = false;
		transcript_type_assigned_bool = false;
		extractUnmapAlignment2ReadFile_bool = false;
		mappedLength_perc_min = 0;
		mappedLength_base_min = 0;
		SNP_provided_bool = false;
		string transcript_type = "NULL";
		SNP_seq_index_provided_bool = false;
		segMap2SNPmer_phase_defined_bool = false;
		sharedThreadForIO_bool = true; // no need to set -J
		readFileIn_listFile_or_commandLine_bool = false;
		sharing_splicing_context_bool = false;
		sharing_SNP_context_bool = true; // no need to set -B
		stopLearningSNP_bool = false;
		stopLearningSplicing_bool = false;
		updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool = false;

		backSplice_search_bool = false;
		#ifdef DETECT_CIRCULAR_RNA
		backSplice_search_bool = true;
		#endif

		fusion_search_bool = false;
		#ifdef MPS_FUSION_POST_NEW
		fusion_search_bool = true;
		#endif

		avoid_creating_folder_bool = false;
		report_sam_only_bool = false;

		keep_tmp_bool = false;
		
		fusion_post_sup_min = 5;

		fusion_post_formatted_gtf_provided_bool = false;
		fusion_post_paralog_gene_provided_bool = false;
	}

	int return_fusion_post_sup_min()
	{
		return fusion_post_sup_min;
	}
	
	string return_fusion_post_formatted_gtf_path()
	{
		return fusion_post_formatted_gtf_path;
	}
	
	string return_fusion_post_paralog_gene_path()
	{
		return fusion_post_paralog_gene_path;
	}

	bool return_keep_tmp_bool()
	{
		return keep_tmp_bool;
	}

	bool return_report_sam_only_bool()
	{
		return report_sam_only_bool;
	}

	bool return_report_junc_bool()
	{
		return (!report_sam_only_bool);
	}

	bool return_avoid_creating_folder_bool()
	{
		return avoid_creating_folder_bool;
	}

	bool return_backSplice_search_bool()
	{
		return backSplice_search_bool;
	}

	bool return_fusion_search_bool()
	{
		return fusion_search_bool;
	}

	bool return_updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool()
	{
		return updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool;
	}

	bool return_updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool()
	{
		return updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool;
	}

	bool return_stopLearningSplicing_bool()
	{
		return stopLearningSplicing_bool;
	}

	bool return_stopLearningSNP_bool()
	{
		return stopLearningSNP_bool;
	}

	bool return_sharing_splicing_context_bool()
	{
		return sharing_splicing_context_bool;
	}

	bool return_sharing_SNP_context_bool()
	{
		return sharing_SNP_context_bool;
	}

	bool return_readFileIn_listFile_or_commandLine_bool()
	{
		return readFileIn_listFile_or_commandLine_bool;
	}

	bool return_sharedThreadForIO_bool()
	{
		return sharedThreadForIO_bool;
	}

	string returnInputSamFilePath()
	{
		return inputSamFilePath;
	}

	bool return_SNP_provided_bool()
	{
		return SNP_provided_bool;
	}

	bool return_SNP_seq_index_provided_bool()
	{
		return SNP_seq_index_provided_bool;
	}

	bool return_extractUnmapAlignment2ReadFile_bool()
	{
		return extractUnmapAlignment2ReadFile_bool;
	}

	bool returnSEorPE_bool()
	{
		// fix me
		return SE_or_PE_bool;
	}

	void getOpt_long(int argc, char**argv)
	{
		struct option long_options[] {
			{"Do_phase1_only", 0, NULL, 'A'},
			//{"Avoid_creating_folder", 0, NULL, '7'},
			{"sharedThreadForIO", 0, NULL, 'B'},
			{"updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bothPhase_or_phase2only", 1, NULL, 'C'},// C
			{"SNPmerMap", 1, NULL, 'D'},
			{"readFileInListFileOrCommandLine", 0, NULL, 'E'},
			{"read_format", 1, NULL, 'F'},
			{"global_index_path", 1, NULL, 'G'},
			{"help", 0, NULL, 'H'},
			{"annotation_path", 1, NULL, 'I'},
			{"stopLearningSNP", 0, NULL, 'J'},
			{"stopLearningSplicing", 0, NULL, 'K'},
			{"local_index_path", 1, NULL, 'L'},
			{"min_mapped_perc", 1, NULL, 'M'},
			{"sharing_SNP_context", 0, NULL, 'N'},
			{"output_path", 1, NULL, 'O'},
			{"SNPfile", 1, NULL, 'P'},
			{"SNP_seq_index_path", 1, NULL, 'Q'},
			{"real_genome_path", 1, NULL, 'R'},		
			{"read_singleEnd_path", 1, NULL, 'S'},
			{"threads", 1, NULL, 'T'},
			{"transcript_path", 1, NULL, 'U'},
			// V
			{"input_SAM_file_path", 1, NULL, 'W'},
			{"sharing_splicing_context", 0, NULL, 'X'},
			{"transcript_type", 1, NULL, 'Y'},
			{"spliceJunctionAlignInferHash_path", 1, NULL, 'Z'},
			{"read_end1_path", 1, NULL, '1'},
			{"read_end2_path", 1, NULL, '2'},
			{"backSplice-search", 0, NULL, '3'},
			{"fusion-search", 0, NULL, '4'},
			{"report-sam-only", 0, NULL, '5'},
			{"keep-tmp", 0, NULL, '6'},
			{"fusion-post_sup_min", 1, NULL, '7'},
			{"fusion-post-gene-ann", 1, NULL, '8'},
			{"fusion-post-paralog-gene", 1, NULL, '9'},
			{NULL, 0, NULL, 0},
		};		
		char ch;
		string threadsStr, formatStr, min_mapped_perc_str, optargStr, SNPmerMap_phase_str, 
			updateChrSeqWithSNP_phase_str, fusion_post_sup_min_str;

		bool option_S_set = false;
		bool option_1_set = false;
		bool option_2_set = false;
		bool option_O_set = false;
		bool option_T_set = false;
		bool option_F_set = false;
		bool option_G_set = false;
		bool option_L_set = false;
		//bool option_C_set = false;
		bool option_D_set = false;

		//cout << "command line: " << endl;
		while( (ch = getopt_long(argc, argv, optionStr.c_str(), long_options, NULL)) != -1 )
		{
		    //printf("optind:%d\n",optind);
		    //printf("optarg:%s\n",optarg);
		    //printf("ch:%c\n",ch);
		    switch(ch)
		    {
		    	case 'Y':
		    		transcript_type = optarg;
		    		if((transcript_type == "GAF")||(transcript_type == "BEER"))
		    		{
		    			transcript_type_assigned_bool = true;
		    		}
		    		else
		    		{
		    			cout << "transcript_type invalid! Must be GAF or BEER" << endl;
		    			exit(1);
		    		}
		    		break;
		    	case 'D':
		    		segMap2SNPmer_phase_defined_bool = true;
		    		SNPmerMap_phase_str = optarg;
		    		if((SNPmerMap_phase_str == "TRUE")||(SNPmerMap_phase_str == "True")||(SNPmerMap_phase_str == "true"))
		    			segMap2SNPmer_bothPhase_or_phase2Only_bool = true;
		    		else if((SNPmerMap_phase_str == "FALSE")||(SNPmerMap_phase_str == "False")||(SNPmerMap_phase_str == "false"))
		    			segMap2SNPmer_bothPhase_or_phase2Only_bool = false;
		    		else
		    		{
		    			cout << "option D for segMap2SNPmer_bothPhase_or_phase2Only_bool must be True or False" << endl;
		    			exit(1);
		    		}
		    		break;
		    	case 'E':
		    		readFileIn_listFile_or_commandLine_bool = true;
		    		break;
		    	//case '7':
		    	//	avoid_creating_folder_bool = true;
		    	//	break;
		    	case 'C':
		    		updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool = true;
		    		updateChrSeqWithSNP_phase_str = optarg;
		    		if((updateChrSeqWithSNP_phase_str == "TRUE")||(updateChrSeqWithSNP_phase_str == "True")||(updateChrSeqWithSNP_phase_str == "true"))
		    			updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool = true;
		    		else if((updateChrSeqWithSNP_phase_str == "FALSE")||(updateChrSeqWithSNP_phase_str == "False")||(updateChrSeqWithSNP_phase_str == "false"))
		    			updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool = false;
		    		else
		    		{
		    			cout << "option C for updateChrSeqWithSNPonly_bothPhase_or_phase2only_bool must be True or False" << endl;
		    			exit(1);
		    		}
		    		break;
		    	case 'X':
		    		sharing_splicing_context_bool = true;
		    		break;
		    	case 'N':
		    		sharing_SNP_context_bool = true;
		    		break;
		    	case 'J':
		    		stopLearningSNP_bool = true;
		    		break;
		    	case 'K':
		    		stopLearningSplicing_bool = true;
		    		break;
		    	case 'Q':
		    		SNP_seq_index_path = optarg;
		    		SNP_seq_index_provided_bool = true;
		    		break;		    		
		    	case 'P':
		    		SNPfilePath = optarg;
		    		SNP_provided_bool = true;
		    		break;
		    	case 'U':
		    		transcript_file_path = optarg;
		    		transcript_provided_bool = true;
		    		break;
		    	case 'R':
		    		realGenome_file_path = optarg;
		    		realGenome_provided_bool = true;
		    		break;
		    	case 'I':
		    		annotation_file_path = optarg;
		    		annotation_provided_bool = true;
		    		break;
		    	case 'Z':
		    		spliceJunctionAlignInferHash_file_path = optarg;
		    		spliceJunctionAlignInferHash_provided_bool = true;
		    		break;
		    	case 'M':
			        min_mapped_perc_str = optarg;
			        mappedLength_perc_min = atoi(min_mapped_perc_str.c_str());		    		
		    		break;
		    	case 'H':
		    		cout << this->optionInfoHelpStr() << endl;
		    		exit(1);
		    		break;
		      	case 'A':
		      		Do_phase1_only_bool = true;
		      		break;
		      	case 'B':
		      		sharedThreadForIO_bool = true;
		      		break;
		      	case 'S':
		      		read_file_path_SE = optarg;
		      		option_S_set = true;
		      		break;
		      	case 'W':
		      		optargStr = optarg;
			        read_file_path_1 = optargStr + ".unmapped.read.1.fq";
			        option_1_set = true;
		        	read_file_path_2 = optargStr + ".unmapped.read.2.fq";
			        option_2_set = true;
			        alreadyMappedAlignmentFile = optargStr + ".mapped.alignment.sam";
			        extractUnmapAlignment2ReadFile_bool = true; 
			        inputSamFilePath = optargStr;  	      	
		      		break;
			    case '1':
			       //printf("option 1:'%s'\n",optarg);
			       read_file_path_1 = optarg;
			       option_1_set = true;
			       break;
			    case '2':
		         	//printf("option 2:'%s'\n",optarg);
		        	read_file_path_2 = optarg;
			        option_2_set = true;
			        break;
			    case 'G':
			        //printf("option G:'%s'\n",optarg);
			        global_index_file_path_prefix = optarg;
			        local_index_file_path_prefix = global_index_file_path_prefix + "/2ndLevelIndex/";
			        option_G_set = true;
			        option_L_set = true;
			        break;
			    // case 'L':
			    //     //printf("option L:'%s'\n",optarg);
			    //     local_index_file_path_prefix = optarg;
			    //     option_L_set = true;
			    //     break;
			    //case 'C':
			        //printf("option C:'%s'\n",optarg);
			        //chromsome_file_path_prefix = optarg;
			        //option_C_set = true;
			        //break;
			    case 'T':
			        //printf("option T:'%s'\n",optarg);
			        threadsStr = optarg;
			        threads_num = atoi(threadsStr.c_str());
			        option_T_set = true;
			        break;
			    case 'O':
			        //printf("option O:'%s'\n",optarg);
			        outputFolder_path = optarg;
			        option_O_set = true;
			        break;
			    case 'F':
			        //printf("option F:'%s'\n",optarg);
			        formatStr = optarg;
			        if((formatStr == "fasta")||(formatStr == "Fasta")||(formatStr == "FASTA"))
			        {
			        	fasta_or_fastq_bool = true;
			        }
			        else if((formatStr == "fastq")||(formatStr == "Fastq")||(formatStr == "FASTQ"))
			        {
			        	fasta_or_fastq_bool = false;
			        }
			        else
			        {
			        	cout << "input format error! Please claim '-F fasta' or '-F fastq' !" << endl;
			        	exit(1);
			        }
			        option_F_set = true;
			        break;
			    case '3':
			    	backSplice_search_bool = true;
			    	break;
			    case '4':
			    	fusion_search_bool = true;
			    	break;
			    case '5':
			    	report_sam_only_bool = true;
			    	break;
			    case '6':
			    	keep_tmp_bool = true;
			    	break;
			    case '7':
			    	fusion_post_sup_min_str = optarg;
					fusion_post_sup_min = atoi(fusion_post_sup_min_str.c_str());
			    	break;					
			    case '8':
					fusion_post_formatted_gtf_path = optarg;
			    	fusion_post_formatted_gtf_provided_bool = true;
			    	break;					
			    case '9':
					fusion_post_paralog_gene_path = optarg;
			    	fusion_post_paralog_gene_provided_bool = true;
			    	break;			
			    default:
			        printf("other option:%c\n",ch);
		    }
		    //printf("optopt+%c\n",optopt);
		}

		if(option_1_set && option_2_set && (!option_S_set))
		{	
			if(readFileIn_listFile_or_commandLine_bool)
			{
				ifstream r1_ifs(read_file_path_1);
				ifstream r2_ifs(read_file_path_2);
				string tmp_1stReadFilePath_read1;
				string tmp_1stReadFilePath_read2;
				getline(r1_ifs, tmp_1stReadFilePath_read1);
				getline(r2_ifs, tmp_1stReadFilePath_read2);
				r1_ifs.close(); r2_ifs.close(); 
				fasta_or_fastq_bool = this->checkInputFileFormat(tmp_1stReadFilePath_read1, tmp_1stReadFilePath_read2);
			}
			else
				fasta_or_fastq_bool = this->checkInputFileFormat(read_file_path_1, read_file_path_2);
			SE_or_PE_bool = false;
		}
		else if((!option_1_set) && (!option_2_set) && option_S_set)
		{
			if(readFileIn_listFile_or_commandLine_bool)
			{
				cout << "for now, batch mode only supports paired end reads" << endl;
				exit(1);
			}
			fasta_or_fastq_bool = this->checkInputFileFormat(read_file_path_SE);
			SE_or_PE_bool = true;
		}
		else
		{
			cout << "invalid read_path set" << endl;
			cout << "option 1 & option 2 or option S should be set" << endl;			
		}

		if(!option_G_set)
			cout << "option G unset" << endl;
		// if(!option_L_set)
		// 	cout << "option L unset" << endl;
		if(!option_T_set)
			cout << "option T unset" << endl;
		if(!option_O_set)
			cout << "option O unset" << endl;
		if( (!((option_1_set && option_2_set && (!option_S_set))||((!option_1_set) && (!option_2_set) && option_S_set)))
			||(!option_O_set)||(!option_T_set)||(!option_G_set)||(!option_L_set))
			exit(1);
		if(!((SNP_seq_index_provided_bool && SNP_provided_bool && segMap2SNPmer_phase_defined_bool)
			||((!SNP_seq_index_provided_bool)&&(!SNP_provided_bool)&&(!segMap2SNPmer_phase_defined_bool))))
		{
			if(updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool)
			{}
			else
			{		
				cout << "options P & Q & D should be set or unset at the same time" << endl;
				exit(1);
			}
		}
		if(updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool)
		{
			if(!SNP_provided_bool)
			{
				cout << "As option C has been set, updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool == true, option P should be set for SNP" << endl;
				exit(1);
			}
			if(SNP_seq_index_provided_bool || segMap2SNPmer_phase_defined_bool)
			{
				cout << "As option C has been set, updateChrSeqWithSNPonly_avoidSegMap2SNPmer_bool == true, option Q & D should not be set for segMap2SNPmerIndex" << endl;
				exit(1);
			}
		}
		if((threads_num == 1)&&(!sharedThreadForIO_bool))
		{
			cout << "IO and process stages have to share threads, when thread_num == 1" << endl;
			cout << "Please set option B" << endl;
			exit(1);
		}

		if(transcript_type_assigned_bool && transcript_provided_bool && realGenome_provided_bool)
		{}
		else if((!transcript_type_assigned_bool)&&(!transcript_provided_bool)&&(!realGenome_provided_bool))
		{}
		else
		{
			cout << "transcriptome prebuilt mode failed! " << endl;
			cout << "transcript_provided_bool: " << transcript_provided_bool << endl;
			cout << "realGenome_provided_bool: " << realGenome_provided_bool << endl;
			cout << "transcript_type_assigned_bool: " << transcript_type_assigned_bool << endl;
			cout << "all the above options should be set 0 or 1 !" << endl;
			exit(1);
		}
		
		#ifdef MPS_FUSION_POST_NEW
		if((!fusion_post_formatted_gtf_provided_bool)||(!fusion_post_paralog_gene_provided_bool))
		{
			cout << "fusion-post-gene-ann and/or fusion-post-paralog-gene not set!" << endl;
			cout << "Please set options " << endl;
			cout << "--fusion-post-gene-ann <fusion-post-gene-ann-path> and/or" << endl;
			cout << "--fusion-post-paralog-gene <fusion-post-paralog-gene-path>" << endl;
			exit(1);	
		}
		#endif
	}

	bool checkInputFileFormat(const string& read_file_path_1)
	{
		int read_file_path_1_len = read_file_path_1.length();
		bool file_1_fasta_or_fastq_bool;
		// check input file name;
		if(read_file_path_1_len > 6)
		{
			string tmpSuffix_len3 = read_file_path_1.substr(read_file_path_1_len-3);
			string tmpSuffix_len6 = read_file_path_1.substr(read_file_path_1_len-6);
			if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa")
				||(tmpSuffix_len6 == ".fasta")||(tmpSuffix_len6 == ".FASTA")||(tmpSuffix_len6 == ".Fasta"))
			{
				file_1_fasta_or_fastq_bool = true;
			}
			else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq")
				||(tmpSuffix_len6 == ".fastq")||(tmpSuffix_len6 == ".FASTQ")||(tmpSuffix_len6 == ".Fastq"))
			{
				file_1_fasta_or_fastq_bool = false;
			}
			else
			{
				cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
				exit(1);
			}
		}
		else if(read_file_path_1_len > 3)
		{
			string tmpSuffix_len3 = read_file_path_1.substr(read_file_path_1_len-3);
			if((tmpSuffix_len3 == ".fa")||(tmpSuffix_len3 == ".FA")||(tmpSuffix_len3 == ".Fa"))
			{
				file_1_fasta_or_fastq_bool = true;
			}
			else if((tmpSuffix_len3 == ".fq")||(tmpSuffix_len3 == ".FQ")||(tmpSuffix_len3 == ".Fq"))
			{
				file_1_fasta_or_fastq_bool = false;
			}
			else
			{
				cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
				exit(1);
			}
		}
		else
		{
			cout << read_file_path_1 << " is not a regular fasta or fastq format file name !\n Please use -H to see detailed information or use -F to specify it." << endl;
			exit(1);			
		}
		return file_1_fasta_or_fastq_bool;				
	}

	bool checkInputFileFormat(const string& read_file_path_1, const string& read_file_path_2) // true -- fasta, false -- fastq
	{
		bool file_1_format_bool = this->checkInputFileFormat(read_file_path_1);
		bool file_2_format_bool = this->checkInputFileFormat(read_file_path_2);
		if(file_1_format_bool != file_2_format_bool)
		{
			cout << read_file_path_1 << " and " << read_file_path_2 << " are in different formats, please check it." << endl;
			exit(1);
		}
		else
		{
			return file_1_format_bool;
		}
	}

	string optionInfoHelpStr()
	{
		string optionStr = "Processing Paired-end reads:\ncommand-line main arguments:\n-1 <string> path_to_read_file_end1\n-2 <string> path_to_read_file_end2\n-G <string> path_to_global_index\n-L <string> path_to_local_index\n-T <int> threads_num\n-O <string> path_to_output_folder\noptional arguments:\n-A Do_phase1_only\n-I annotation file\n-F <string> input_format_fasta_or_fastq\n-H command line help information\n";
		optionStr = optionStr + "Processing Single-end reads:\ncommand-line main arguments:\n-S <string> path_to_read_file_SE\n-G <string> path_to_global_index\n-L <string> path_to_local_index\n-T <int> threads_num\n-O <string> path_to_output_folder\noptional arguments:\n-A Do_phase1_only\n-I annotation file\n-F <string> input_format_fasta_or_fastq\n-H command line help information\n";
		
		return optionStr;
	}

	void outputOptStr(ofstream& output_ofs)
	{
		output_ofs << "read_file_path_1: " << read_file_path_1 << "\nread_file_path_2: " << read_file_path_2
			<< "\nread_file_path_SE: " << read_file_path_SE 
			<< "\nthreads_num: " << threads_num << "\ninput_format_fasta_or_fastq: " << fasta_or_fastq_bool 
			<< "\nglobal_index_file_path_prefix: " << global_index_file_path_prefix << "\nlocal_index_file_path_prefix: "
			<< local_index_file_path_prefix << "\noutputFolder_path: " << outputFolder_path << "\nmin_mapped_perc: " 
			<< mappedLength_perc_min << "\n\nannotation_provided_bool: " << annotation_provided_bool << "\nannotation_file_path: " 
			<< annotation_file_path << endl << "\nspliceJunctionAlignInferHash_provided_bool: "
			<< spliceJunctionAlignInferHash_provided_bool << "\nspliceJunctionAlignInferHash_file_path: "
			<< spliceJunctionAlignInferHash_file_path << endl;
	}

	void outputSwitchInfo(bool Do_Phase1_Only, bool outputAlignInfoAndSamForAllPairedAlignmentBool,
		bool removeAllIntermediateFilesBool, bool Do_cirRNA, bool outputDirectlyBool_Phase1Only, 
		int normalRecordNum_1stMapping, int normalRecordNum_fixOneEndUnmapped,
		int normalRecordNum_fixHeadTail, bool Do_extendHeadTail_phase1, 
		bool Do_extendHeadTail_fixOneEndUnmapped, bool Do_extendHeadTail_fixHeadTail, 
		bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
		bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
		ofstream& log_ofs)
	{
		log_ofs << endl << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
		log_ofs << endl << "Do_Phase1_Only: " << Do_Phase1_Only << endl << "outputAlignInfoAndSamForAllPairedAlignmentBool: " 
			<< outputAlignInfoAndSamForAllPairedAlignmentBool << endl << "removeAllIntermediateFilesBool: " 
			<< removeAllIntermediateFilesBool << endl << "Do_cirRNA: " << Do_cirRNA << endl
			<< "outputDirectlyBool_Phase1Only: " << outputDirectlyBool_Phase1Only << endl; ;
		log_ofs << endl << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
		log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
		log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;

		log_ofs << endl << "Do_extendHeadTail_phase1: " << Do_extendHeadTail_phase1 << endl;
		log_ofs << endl << "Do_extendHeadTail_fixOneEndUnmapped: " << Do_extendHeadTail_fixOneEndUnmapped << endl;
		log_ofs << endl << "Do_extendHeadTail_fixHeadTail: " << Do_extendHeadTail_fixHeadTail << endl; 

		log_ofs << endl << "Do_fixHeadTail_remapping: " << Do_fixHeadTail_remapping << endl;
		log_ofs << endl << "Do_fixHeadTail_greedyMapping: " << Do_fixHeadTail_greedyMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAndTargetMapping: " << Do_fixHeadTail_remappingAndTargetMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAgain: " << Do_fixHeadTail_remappingAgain << endl;

		log_ofs << endl << "minValSegLength: " << minValSegLength << endl;
		log_ofs << endl << "min_anchor_length: " << min_anchor_length << endl;
		log_ofs << endl << "CONFIDENT_SEG_LENGTH_FIX_LONG_END " << CONFIDENT_SEG_LENGTH_FIX_LONG_END << endl;
		log_ofs << endl << "DETECT_NONCANONICAL_SJ: " << DETECT_NONCANONICAL_SJ << endl;
	}

	void outputSwitchInfo(bool Do_Phase1_Only, bool outputAlignInfoAndSamForAllPairedAlignmentBool,
		bool removeAllIntermediateFilesBool, bool Do_cirRNA, bool outputDirectlyBool_Phase1Only, 
		//int normalRecordNum_1stMapping, 
		//int normalRecordNum_fixOneEndUnmapped,
		//int normalRecordNum_fixHeadTail, 
		bool Do_extendHeadTail_phase1, 
		bool Do_extendHeadTail_fixOneEndUnmapped, bool Do_extendHeadTail_fixHeadTail, 
		bool Do_fixHeadTail_remapping, bool Do_fixHeadTail_greedyMapping,
		bool Do_fixHeadTail_remappingAndTargetMapping, bool Do_fixHeadTail_remappingAgain,
		ofstream& log_ofs)
	{
		log_ofs << endl << "SE_or_PE_bool: " << SE_or_PE_bool << endl;
		log_ofs << endl << "Do_Phase1_Only: " << Do_Phase1_Only << endl << "outputAlignInfoAndSamForAllPairedAlignmentBool: " 
			<< outputAlignInfoAndSamForAllPairedAlignmentBool << endl << "removeAllIntermediateFilesBool: " 
			<< removeAllIntermediateFilesBool << endl << "Do_cirRNA: " << Do_cirRNA << endl
			<< "outputDirectlyBool_Phase1Only: " << outputDirectlyBool_Phase1Only << endl; ;
		//log_ofs << endl << "normalRecordNum_1stMapping: " << normalRecordNum_1stMapping << endl;
		//log_ofs << "normalRecordNum_fixOneEndUnmapped: " << normalRecordNum_fixOneEndUnmapped << endl;
		//log_ofs << "normalRecordNum_fixHeadTail: " << normalRecordNum_fixHeadTail << endl;

		log_ofs << endl << "Do_extendHeadTail_phase1: " << Do_extendHeadTail_phase1 << endl;
		log_ofs << endl << "Do_extendHeadTail_fixOneEndUnmapped: " << Do_extendHeadTail_fixOneEndUnmapped << endl;
		log_ofs << endl << "Do_extendHeadTail_fixHeadTail: " << Do_extendHeadTail_fixHeadTail << endl; 

		log_ofs << endl << "Do_fixHeadTail_remapping: " << Do_fixHeadTail_remapping << endl;
		log_ofs << endl << "Do_fixHeadTail_greedyMapping: " << Do_fixHeadTail_greedyMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAndTargetMapping: " << Do_fixHeadTail_remappingAndTargetMapping << endl;
		log_ofs << endl << "Do_fixHeadTail_remappingAgain: " << Do_fixHeadTail_remappingAgain << endl;

		log_ofs << endl << "minValSegLength: " << minValSegLength << endl;
		log_ofs << endl << "min_anchor_length: " << min_anchor_length << endl;
		log_ofs << endl << "CONFIDENT_SEG_LENGTH_FIX_LONG_END " << CONFIDENT_SEG_LENGTH_FIX_LONG_END << endl;
		log_ofs << endl << "DETECT_NONCANONICAL_SJ: " << DETECT_NONCANONICAL_SJ << endl;		
	}
};

#endif