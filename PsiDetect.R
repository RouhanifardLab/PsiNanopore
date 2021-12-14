#### ========================================== ####
####           Main: p-value calculation        ####
#### ========================================== ####
#!/usr/bin/env Rscript
wdir = getwd( ) 
print(wdir)

print("Initializing the package...")
source(paste(wdir,"/R/MakeOptions_psiDetect.R",sep = ""))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

wdir = getwd( ) 

acc.bases = c("A","T","C","G")
print("------------------------------------")

print("Loading input files...")

kmer.summary.df = read.csv( as.character(opt$file3) )#"~/Rouhanifard lab Dropbox/Sepideh/DataSets_Nanopore/SH_SY5Y/kmer_summary_IVT_for_Pvalue/kmer_summary_5Y_IVT_for_Pvalu.csv")
fasta_file <- FaFile( as.character(opt$file4) )#file='~/Rouhanifard lab Dropbox/Sepideh/Genome_References/Genome_human/GRCh38.p10.genome.fa')
Direct1.bam.dir = as.character(opt$file) #"~/Desktop/Github_Material/PSMB2_hela_D1_hg38.bam"
ctrl.bam.dir    = as.character(opt$file2)#"~/Desktop/Github_Material/PSMB2_hela_IVT1_hg38.bam"
print("------------------------------------")
print("Data Preparation...")

D1b <- BamFile(Direct1.bam.dir)
ctrlb <- BamFile(ctrl.bam.dir)

target.chr = as.character(opt$chromosome)

pos1 = as.numeric(opt$start_position)
pos2 = as.numeric(opt$end_position)

param <- ScanBamParam(which=GRanges(target.chr,IRanges(start=pos1, 
                                                       end=pos2)))

pilup_params =  Rsamtools::PileupParam(max_depth = 10250,min_mapq = 1,distinguish_nucleotides = T,
                                       ignore_query_Ns = T)
PU = pileup(D1b, scanBamParam=param,pileupParam = pilup_params)

PU$nucleotide.new = PU$nucleotide
PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="A")] = "T"
PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="T")] = "A"
PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="C")] = "G"
PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="G")] = "C"
PU = PU[which(PU$nucleotide %in% acc.bases),]

PU.ctrl = pileup(ctrlb, scanBamParam=param,pileupParam = pilup_params)

PU.ctrl$nucleotide.new = PU.ctrl$nucleotide
PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="A")] = "T"
PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="T")] = "A"
PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="C")] = "G"
PU.ctrl$nucleotide.new[which(PU.ctrl$strand == "-" & PU.ctrl$nucleotide =="G")] = "C"

res.df = data.frame("pos" = unique(PU$pos), "A.count"=0,"C.count"=0,"G.count"=0,"T.count"=0,
                    "A.count.ctrl"=0,"C.count.ctrl"=0,"G.count.ctrl"=0,"T.count.ctrl"=0,
                    "reference"="","kmer"="","strand"="+")
print("------------------------------------")

print("Direct seq Pileup started...")
for (i in c(1:nrow(PU))){
  if (i%%1000 == 0){print(paste("Pileup: ",round(i/nrow(PU)*100,digits = 2),"%"))}
  ro.res.df = which(res.df$pos == PU$pos[i])
  
  col = 2
  if (PU$nucleotide.new[i] == "C"){col=3}
  if (PU$nucleotide.new[i] == "G"){col=4}
  if (PU$nucleotide.new[i] == "T"){col=5}
  res.df$strand[ro.res.df] = as.character(PU$strand[i])
  res.df[ro.res.df, col] = res.df[ro.res.df, col] + PU$count[i]
}
print("------------------------------------")
print("Control Pileup started...")
for (i in c(1:nrow(PU.ctrl))){
  if (i%%1000 == 0){print(paste("Control Pileup: ",round(i/nrow(PU.ctrl)*100,digits = 2),"%"))}
  
  ro.res.df = which(res.df$pos == PU.ctrl$pos[i])
  col = 6
  if (PU.ctrl$nucleotide.new[i] == "C"){col=7}
  if (PU.ctrl$nucleotide.new[i] == "G"){col=8}
  if (PU.ctrl$nucleotide.new[i] == "T"){col=9}
  res.df[ro.res.df, col] = res.df[ro.res.df, col] + PU.ctrl$count[i]
}

print("Done with pileup!")
print("------------------------------------")
print("Started assigning kmers...")

for (i in c(1:nrow(res.df))){
  if (i%%100 == 0){print(paste("Assigning kmers: ",round(i/nrow(res.df)*100,digits = 2),"%"))}
  
  kmer = ""
  for (j in c(-2:2)){
    gr1 <- GRanges(target.chr,IRanges(start=res.df$pos[i]+j, end=res.df$pos[i]+j))
    ### Extract the kmers
    refbase <- getSeq(fasta_file, gr1)
    refbase <- as.data.frame(refbase)$x
    refbase.new = refbase
    if (refbase == "A" & res.df$strand[i] == "-"){res.df.new = "T"}
    if (refbase == "T" & res.df$strand[i] == "-"){res.df.new = "A"}
    if (refbase == "C" & res.df$strand[i] == "-"){res.df.new = "G"}
    if (refbase == "G" & res.df$strand[i] == "-"){res.df.new = "C"}
    
    kmer = paste(kmer,res.df.new,sep = "")
    if (j == 0){
      res.df$reference[i] = res.df.new
    }
    
  }
  res.df$kmer[i] = kmer
}
print("Done with assigning kmers!")
print("------------------------------------")
print("Calculating p-values...")

res.df.U = res.df[which(res.df$reference == "T"),]
res.df.U$NReads = res.df.U$C.count + res.df.U$T.count
res.df.U = res.df.U[which(res.df.U$NReads >0),]
res.df.U$mm.perc = res.df.U$C.count / (res.df.U$T.count + res.df.U$C.count)*100
res.df.U$ctrl.err = res.df.U$C.count.ctrl / (res.df.U$T.count.ctrl + res.df.U$C.count.ctrl)*100

res.df.U = res.df.U[which(res.df.U$mm.perc >0 & res.df.U$ctrl.err>=0),]


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

res.df.U$kmer.err = 0
for (i in c(1:nrow(kmer.summary.df))) {
  res.df.U$kmer.err[which(res.df.U$kmer == kmer.summary.df$kmer[i])] = kmer.summary.df$avg.mm.ivt[i]
}

res.df.U$expected.err = res.df.U$kmer.err
res.df.U$expected.err[which(res.df.U$ctrl.err>res.df.U$kmer.err)] = res.df.U$ctrl.err[which(res.df.U$ctrl.err>res.df.U$kmer.err)]

res.df.U$p=1
for (i in c(1:nrow(res.df.U))){
  if (i %% 1000 == 0) {
    print(i/nrow(res.df.U)*100)
  }
  if (res.df.U$mm.perc[i]>res.df.U$expected.err[i]){
    res.df.U$p[i] = calc.p.val(n.read = res.df.U$NReads[i],
                                        mod.prob = res.df.U$expected.err[i]/100, 
                                        subject.mm = res.df.U$mm.perc[i]/100)
    counter = 1
    step = 100
    while (is.na(res.df.U$p[i])){
      res.df.U$p[i] = calc.p.val(n.read = res.df.U$NReads[i]-counter*step,
                                          mod.prob = res.df.U$expected.err[i]/100, 
                                          subject.mm = res.df.U$mm.perc[i]/100)
      counter =  counter+1
    }
  }
}

res.df.U$chr = target.chr
res.df.U = res.df.U[which(res.df.U$p <= as.numeric(opt$max_p)),]
res.df.U = res.df.U[,c(ncol(res.df.U),c(1:(ncol(res.df.U)-1)))]
write.csv(res.df.U, opt$out,row.names = F)

print("Done! Please find the output file here:")
print(opt$out)


