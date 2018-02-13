#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 outputScript" << endl;
		exit(1);
	}
	string outputScript = argv[1];
	int dataSetNum = 74;
	ofstream script_ofs(outputScript.c_str());
	string cmd_cd_hisat2dir = "cd /home/lcph222/Xinan/tools/hisat2-2.0.5/";
	script_ofs << cmd_cd_hisat2dir << endl << endl;

	for(int tmpDataNO = 1; tmpDataNO <= dataSetNum; tmpDataNO ++)
	{	
		// hisat2 genome
		script_ofs << "./hisat2 -q -p 16 --no-mixed --no-discordant -x /home/xli262/chrom_Index/HISAT2_index/grch37/genome \\" << endl;
		script_ofs << "-1 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/1_paired.fastq \\" << endl;
		script_ofs << "-2 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/2_paired.fastq \\" << endl;
		script_ofs << "-S /scratch/lcph222/Xinan/personalGenome/1000genomeProject/HISAT2_results/sam_PEonly/" << tmpDataNO << ".sam" << endl << endl;
		// hisat2 genome + population SNPs
		script_ofs << "./hisat2 -q -p 16 --no-mixed --no-discordant -x /home/xli262/chrom_Index/HISAT2_index/grch37_snp/genome_snp \\" << endl;
		script_ofs << "-1 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/1_paired.fastq \\" << endl;
		script_ofs << "-2 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/2_paired.fastq \\" << endl;
		script_ofs << "-S /scratch/lcph222/Xinan/personalGenome/1000genomeProject/HISAT2_results/sam_PEonly/" << tmpDataNO << "_pSNP.sam" << endl << endl;
		// hisat2 genome + individual SNPs
		script_ofs << "./hisat2 -q -p 16 --no-mixed --no-discordant -x /scratch/lcph222/Xinan/personalGenome/1000genomeProject/HISAT2_index/" << tmpDataNO << "/KGP_" << tmpDataNO << " \\" << endl;
		script_ofs << "-1 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/1_paired.fastq \\" << endl;
		script_ofs << "-2 /scratch/lcph222/Xinan/personalGenome/1000genomeProject/fastq/" << tmpDataNO << "/2_paired.fastq \\" << endl;
		script_ofs << "-S /scratch/lcph222/Xinan/personalGenome/1000genomeProject/HISAT2_results/sam_PEonly/" << tmpDataNO << "_iSNP.sam" << endl << endl;
	}
	script_ofs.close();
	return 0;
}