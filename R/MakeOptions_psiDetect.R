#### ====================================== ####
####    Make/define the main code options   ####
#### ====================================== ####
#wdir = getwd( ) 
source(paste(wdir,"/R/HandlePackages.R",sep = ""))

option_list = list(
  make_option(c("-f", "--file"), type="character", default="", 
              help="path to direct seq. bam file, Ex: ~/Downloads/PsiNanopore-main/data/DirSeq.bam", metavar="character"),
  make_option(c("-g", "--file2"), type="character", default="", 
              help="path to control bam file, Ex: ~/Downloads/PsiNanopore-main/data/IVT.bam", metavar="character"),
  make_option(c("-k", "--file3"), type="character", default="", 
              help="path to kmer-dependent error file (i.e. ~/Downloads/PsiNanopore-main/data/kmer_summary.csv)", metavar="character"),
  make_option(c("-r", "--file4"), type="character", default="", 
              help="Reference genome fasta file Ex: ~/Downloads/GRCh38.p10.genome.fa (You can download this file here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.36/)", metavar="character"),
  make_option(c("-s", "--start_position"), type="character", default=NULL, 
              help="start position, Ex: 35599541", metavar="character"),
  make_option(c("-e", "--end_position"), type="character", default=NULL, 
              help="end position, Ex: 35641526", metavar="character"),
  make_option(c("-c", "--chromosome"), type="character", default=NULL, 
              help="chromosome, EX: chr1", metavar="character"),
  make_option(c("-m", "--max_p"), type="character", default=NULL, 
              help="Sites with a p-value this low or smaller will be outputted, Ex: 0.05", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="psi_candidates.csv", 
              help="output file name of psi candidate sites, Ex: ~/Desktop/psi_candidates.csv", metavar="character")
); 