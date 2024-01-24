# Step 4

# use IVT data of all genes, to calculate the expected error rates


rm(list = ls());cat("\014")  

print('step4')


# find all possible kmers with a "T" in the middle
kmer.lst = c()
nuc.lst = c("A","C",'T',"G")

for (i1 in c(1:4)){
  for (i2 in c(1:4)){
    for (i3 in c(1:4)){
      for (i4 in c(1:4)){
        kmer = paste0(nuc.lst[i1], nuc.lst[i2], 'T', nuc.lst[i3], nuc.lst[i4])
        kmer.lst = c(kmer.lst, kmer)
      }
    }
  }
}
3

kmer.lst = kmer.lst[order(kmer.lst)]
U.count = rep(0, length(kmer.lst)); names(U.count) = kmer.lst
C.count = rep(0, length(kmer.lst)); names(C.count) = kmer.lst
found.sites = rep(0, length(kmer.lst)); names(found.sites) = kmer.lst

################################################################################
# path to merged raw pileups for each gene
data.dir                    = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/Processed_pileup/"
#CREATE kmer_analysis AND finalized_pileup FOLDER FIRST!
k.mer.summary.out.dir       = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/kmer_analysis/IVT_kmer_Analysis.csv"
merged.Uonly.Pileup.out.dir = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/finalized_pileup/merged_U_only.csv"

# get the list of the genes in the input directory
genes = list.files(data.dir)

################################################################################
 target.cols = c("Annotation", "chr",        "position",   "strand",
                 "target.nucleotide" ,
                 "T_rep1",    "C_rep1",
                 "T_rep2",    "C_rep2",
                 "T_rep3",    "C_rep3",
                 "T_IVT",      "C_IVT",    "kmer"   )


# loop over all genes in the 'Merged' folder
#
for (k in c(1:length(genes))){
  gene = genes[k] # get the name of the gene
  
  # verbose: report progress in percentage
  
  if (k%%100 == 0) {
    print(paste0(round(k/length(genes)*100, 2), "%"))
    kmer.analysis.df$mm = kmer.analysis.df$C / (kmer.analysis.df$U + kmer.analysis.df$C) * 100
    write.csv(kmer.analysis.df, k.mer.summary.out.dir, row.names = F)
    write.csv(merged.df, merged.Uonly.Pileup.out.dir, row.names = F)
    
  }
  
  
  raw.pileup.dir = paste0(data.dir, gene)
  raw.pileup.df = read.csv(raw.pileup.dir)
  
  gene.name = strsplit(gene,".csv")[[1]]
  raw.pileup.df$Annotation = gene.name
  raw.pileup.df = raw.pileup.df[, c(ncol(raw.pileup.df), c(1:(ncol(raw.pileup.df)-1)))]
  
  for (i in c(1:nrow(raw.pileup.df))){
    kmer = raw.pileup.df$kmer[i]
    if (is.na(kmer)) kmer =""
    if (kmer!="" & raw.pileup.df$target.nucleotide[i] == "T"){
      found.sites[kmer] = found.sites[kmer] + 1
      U.count[kmer] = U.count[kmer] + raw.pileup.df$T_JurkatIVT[i]
      C.count[kmer] = C.count[kmer] + raw.pileup.df$C_JurkatIVT[i]
    }
  }
  
  kmer.analysis.df = data.frame("kmer"=names(U.count), 
                                "U"=as.numeric(U.count), 
                                "C"=as.numeric(C.count), 
                                "n.sampled.sites"=as.numeric(found.sites))
  
  if (length(which(is.na(raw.pileup.df$kmer) == T))>0){
    raw.pileup.df$kmer[which(is.na(raw.pileup.df$kmer))] = ""
  }
  
  raw.pileup.df = raw.pileup.df[which(raw.pileup.df$kmer!= "" & raw.pileup.df$target.nucleotide == "T"),target.cols]
  if (nrow(raw.pileup.df)>0){ 
    if (exists('merged.df')){
      merged.df = rbind(merged.df, raw.pileup.df)
    } else {
      merged.df = raw.pileup.df
    }
  }
}

write.csv(kmer.analysis.df, k.mer.summary.out.dir, row.names = F)
write.csv(merged.df, merged.Uonly.Pileup.out.dir, row.names = F)




