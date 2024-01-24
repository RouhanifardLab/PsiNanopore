# Step 3

# 1 - add 5 mers to the merged dataframe of each gene
# 2 - for (-) strand genes, reverse target nucleotide, and RC 5mer
# 3 - for (-) strand genes, apply RC on pileup values too
# 3 - keep only uridines, calculate U-C mm 
# 4 - use IVT data of all genes, to calculate the expected error rates
# 5 - calculate p-values

rm(list = ls());cat("\014")  

####################################### path to merged raw pileups for each gene
data.dir = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/Merged/"

# path to save the processed data
####################### YOU NEED TO CREATE THIS FOLDER BEFORE RUNNING THE SCRIPT
out.dir = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/Processed_pileup/"

# get the list of the genes in the input directory
genes = list.files(data.dir)

print('step3')

# loop over all genes in the 'Merged' folder
for (k in c(1:length(genes))){
  gene = genes[k] # get the name of the gene
  
  # verbose: report progress in percentage
  if (k%%100 == 0) print(paste0(round(k/length(genes)*100, 2), "%"))
  
  # name of the processed file to save
  out.file.name = paste0(out.dir,gene)
  
  # proceed only if the current gene is not already processed
  if (file.exists(out.file.name) == F){
    # read input file (raw pileup) for the current gene
    raw.pileup.dir = paste0(data.dir, gene)
    raw.pileup.df = read.csv(raw.pileup.dir)
    
    # find the strand of the current gene
    gene.strand = raw.pileup.df$strand[1] 
    
    # on (-) strand, reverse nucleotides
    if (gene.strand == "-"){
      raw.nucleotides = raw.pileup.df$target.nucleotide
      corrected.nucleotides = raw.nucleotides
      corrected.nucleotides[which(raw.nucleotides == "A")] = "T"
      corrected.nucleotides[which(raw.nucleotides == "T")] = "A"
      corrected.nucleotides[which(raw.nucleotides == "C")] = "G"
      corrected.nucleotides[which(raw.nucleotides == "G")] = "C"
      raw.pileup.df$target.nucleotide = corrected.nucleotides
      
      
      # correct pileup values too
      raw.pileup.df.2 = raw.pileup.df
      
################################################################################      
     
      raw.pileup.df.2$A_rep1 = raw.pileup.df$T_rep1
      raw.pileup.df.2$T_rep1 = raw.pileup.df$A_rep1
      raw.pileup.df.2$C_rep1 = raw.pileup.df$G_rep1
      raw.pileup.df.2$G_rep1 = raw.pileup.df$C_rep1
       
      raw.pileup.df.2$A_rep2 = raw.pileup.df$T_rep2
      raw.pileup.df.2$T_rep2 = raw.pileup.df$A_rep2
      raw.pileup.df.2$C_rep2 = raw.pileup.df$G_rep2
      raw.pileup.df.2$G_rep2 = raw.pileup.df$C_rep2
      
      raw.pileup.df.2$A_rep3 = raw.pileup.df$T_rep3
      raw.pileup.df.2$T_rep3 = raw.pileup.df$A_rep3
      raw.pileup.df.2$C_rep3 = raw.pileup.df$G_rep3
      raw.pileup.df.2$G_rep3 = raw.pileup.df$C_rep3

      raw.pileup.df.2$A_IVT = raw.pileup.df$T_IVT
      raw.pileup.df.2$T_IVT = raw.pileup.df$A_IVT
      raw.pileup.df.2$C_IVT = raw.pileup.df$G_IVT
      raw.pileup.df.2$G_IVT = raw.pileup.df$C_IVT 

    }
    
    # define a new column, name it 'kmer'
    raw.pileup.df$kmer =""
    # add kmer for each position
    if (nrow(raw.pileup.df)>6){
      for (i in c(3:(nrow(raw.pileup.df) - 2))){
        #merged.df$target.nucleotide
        
        pos.n2 = raw.pileup.df$position[i-2]
        pos.p2 = raw.pileup.df$position[i+2]
        #  if we already have the required data, get the kmer for the current position
        if (pos.p2 == (pos.n2 + 4)){
          kmer = ""
          adjacent.bases = c(-2:2)
          if (gene.strand == "-") adjacent.bases = c(2:(-2))
          for (j in adjacent.bases){
            kmer = paste0(kmer, raw.pileup.df$target.nucleotide[i + j])
          }
          raw.pileup.df$kmer[i] = kmer
        }
      }
    }

    # once done, write the resulting dataframe
    write.csv(raw.pileup.df, out.file.name, row.names = F)
  }
}

