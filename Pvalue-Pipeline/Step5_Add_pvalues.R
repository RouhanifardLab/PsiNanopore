# read the kmer summary data frame (calculated from IVT)

print("Reading kmer summary ... step5")
################################################################################
kmer.summary = read.csv("/CellLinesPVal/SH-SY5Y/init_gene_pileup/kmer_analysis/IVT_kmer_Analysis.csv")
kmer.summary$mm = kmer.summary$C / (kmer.summary$U + kmer.summary$C) * 100

merged.df = read.csv("/CellLinesPVal/SH-SY5Y/init_gene_pileup/finalized_pileup/merged_U_only.csv")

################################################################################
# calculate mismatches
 merged.df$mm.rep1   = merged.df$C_rep1 / (merged.df$T_rep1 + merged.df$C_rep1) * 100
 merged.df$mm.rep2   = merged.df$C_rep2 / (merged.df$T_rep2 + merged.df$C_rep2) * 100
 merged.df$mm.rep3   = merged.df$C_rep3 / (merged.df$T_rep3 + merged.df$C_rep3) * 100
 merged.df$mm.IVT  = merged.df$C_IVT / (merged.df$T_IVT + merged.df$C_IVT) * 100

################################################################################
# replace 'NA's with 0
 merged.df$mm.rep1[which(is.na(merged.df$mm.rep1))] = 0
 merged.df$mm.rep2[which(is.na(merged.df$mm.rep2))] = 0
 merged.df$mm.rep3[which(is.na(merged.df$mm.rep3))] = 0
 merged.df$mm.IVT[which(is.na(merged.df$mm.IVT))] = 0


# add the expected mm column (initialize with IVT mm %)
print("Calculating the expected error rate for each position ...")
#merged.df$expected.mm = merged.df$mm.ActivatedIVT 
merged.df$expected.mm = merged.df$mm.JurkatIVT 

#summary(merged.df)
# if kmer mm% is larger than IVT mm%, take that as the frame of reference, i.e. expected mm
for (i in c(1:nrow(kmer.summary))){
  expected.mm = kmer.summary$mm[i]
  kmer = kmer.summary$kmer[i]
  print(paste0(kmer, " - ", i, ' out of ', nrow(kmer.summary)))
  
  merged.df$expected.mm[
    which(merged.df$kmer == kmer & 
            merged.df$expected.mm<expected.mm)] = 
    expected.mm
}

#summary(merged.df)

# calculate p-values
print("Calculating p-values ...")

#  Function to calculate p-value
calc.p.val = function(n.read, mod.prob, subject.mm){
  n.err = round(subject.mm * n.read)
  n.err.lst = c(n.err:n.read)
  p.val = 0
  for (n.err in n.err.lst){
    comb = exp(lgamma(n.read+1) - lgamma(n.read-n.err+1)-lgamma(n.err+1))
    prob.analytical = comb * (mod.prob)^n.err *(1-mod.prob)^(n.read-n.err)
    p.val = p.val + prob.analytical
   
  }
  return (p.val)
}

################################################################################
# initialize a column to store calculated p-vlaues FOR ALL REPS EXCEPT IVT
 merged.df$p.value.rep1 = 1
 merged.df$p.value.rep2 = 1
 merged.df$p.value.rep3 = 1


min.acc.reads = 8

period = round(nrow(merged.df) / 1000)

print(nrow(merged.df))
for (i in c(1:nrow(merged.df))){#
  # verbose: report progress
  if (i%%period == 0) print(paste0(round(i/nrow(merged.df) * 100, 2),"%"))
  
  
################################################################################

   # rep1
   if ((merged.df$C_rep1[i] + merged.df$T_rep1[i]) >= min.acc.reads &  merged.df$mm.rep1[i]>merged.df$expected.mm[i]){
     merged.df$p.value.rep1[i] =
       calc.p.val(n.read = merged.df$C_rep1[i] + merged.df$T_rep1[i],
                  mod.prob = merged.df$expected.mm[i]/100, subject.mm = merged.df$mm.rep1[i]/100 )
     
  }
   
   # rep2
   if ((merged.df$C_rep2[i] + merged.df$T_rep2[i]) >= min.acc.reads &  merged.df$mm.rep2[i]>merged.df$expected.mm[i]){
     merged.df$p.value.rep2[i] =
       calc.p.val(n.read = merged.df$C_rep2[i] + merged.df$T_rep2[i],
                  mod.prob = merged.df$expected.mm[i]/100, subject.mm = merged.df$mm.rep2[i]/100 )
     
   }
   
  # # rep3
   if ((merged.df$C_rep3[i] + merged.df$T_rep3[i]) >= min.acc.reads &  merged.df$mm.rep3[i]>merged.df$expected.mm[i]){
     merged.df$p.value.rep3[i] =
       calc.p.val(n.read = merged.df$C_rep3[i] + merged.df$T_rep3[i],
                  mod.prob = merged.df$expected.mm[i]/100, subject.mm = merged.df$mm.rep3[i]/100 )
     
   }
  
 
  
  
  
}

#summary(merged.df[c(1:1000),])
write.csv(merged.df,"/CellLinesPVal/SH-SY5Y/init_gene_pileup/finalized_pileup/Merged_with_P_vals.csv", row.names = F)



 
