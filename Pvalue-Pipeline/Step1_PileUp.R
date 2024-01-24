# Objective: Analyze a bam file
##          - generate pile-up for all genes
##          - calculate p-value based on U->C mm for all uridines

rm(list = ls());cat("\014")  

############ LOAD LIBRARIES ############ 
library(Rsamtools)
library(rtracklayer)
library(BSgenome)
############ INPUTS ############ 
# So far, these are completed: Diff1, IVT
# Undiff 1 is running

#meta.data.df = data.frame("replicate.name"=c("PB1","PB2","Undiff1","Undiff2","Undiff3","IVT"),
#                          "bam.name" = c("SH-SY5Y_PB_Undiff_Direct_rep1.bam",
#                                         "SH-SY5Y_PB_Undiff_Direct_rep2.bam",
#                                         "SH-SY5Y_Undiff_Direct_rep1.bam",
#                                         "SH-SY5Y_Undiff_Direct_rep2.bam",
#                                         "SH-SY5Y_Undiff_Direct_rep3.bam",
#                                         "SH-SY5Y_Diff_IVT_rep1.bam"))

#meta.data.df = data.frame("replicate.name"=c("rep1","rep2","rep3"),
#                          "bam.name" = c("SH-SY5Y_Undiff_Direct_rep1.hg38v10.bam","SH-SY5Y_Undiff_Direct_rep2.hg38v10.bam","SH-SY5Y_Undiff_Direct_rep3.hg38v10.bam"))

meta.data.df = data.frame("replicate.name"=c("rep1"),
                          "bam.name" = c("SH-SY5Y_Undiff_Direct_rep1.hg38v10.bam"))

hg38.fa = FaFile("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/GRCh38.p10.genome.fa")
g = readGFF("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/gencode/gencode.v27.annotation.gff3")

# keep only protein_coding genes, and whole genes only
# dump chrY, because sample is female
g = g[which(g$gene_type == "protein_coding" & 
              g$type == "gene"  & 
              g$seqid!="chrY"),]

genes.df = data.frame("Annotation"=as.character(g$gene_name),
                      "chr"=as.character(g$seqid),
                      "start.point"=g$start,
                      "end.point"=g$end,
                      "strand"=as.character(g$strand),
                      stringsAsFactors = FALSE)

# remove duplicate rows from genes.df
genes.df$delete = F
for (i in c(1:nrow(genes.df))){
  if (genes.df$delete[i] == F){
    rows = which(genes.df$Annotation == genes.df$Annotation[i])
    if (length(rows) > 1){
      rows = rows[-1]
      genes.df$delete[rows] = T
    }
  }
}
genes.df = genes.df[which(genes.df$delete == F),]
genes.df = genes.df[, -ncol(genes.df)]



acceptable.bases = c("T", "A", "G", "+", "-", "C")

for (k in c(1:nrow(meta.data.df))){
  
  ############################################################################################################################################ modify
  ##################################create init_gene_pileup folder and all the replicate name folders########################################
  bam.dir = paste0("/CellLinesPVal/SH-SY5Y/bam_files/",meta.data.df$bam.name[k])
  #out.dir = paste0("~/repository/2023/Nanopore_workflow/data/init_gene_pileup/",meta.data.df$replicate.name[k],"/")
  out.dir = paste0("/CellLinesPVal/SH-SY5Y/init_gene_pileup/",meta.data.df$replicate.name[k],"/")
  
  Rsamtools::indexBam(bam.dir)
  min.acceptable.reads = 8 # minimum number of reads on a position to be eligible as significant
  min.mapq = 20
  
  bamfile <- BamFile(bam.dir)
  
  
  print("starting the loop over all genes")
  for (i in c(1:nrow(genes.df))){
    curr.strand = genes.df$strand[i]
    file.name = paste0(out.dir,"/", genes.df$Annotation[i],".csv")
    
    # skip this iteration if the output file already exists!
    if (file.exists(file.name) == F){
      # pile up using RSamtools functions
      param <- ScanBamParam(which=GRanges(strand = genes.df$strand[i],
                                          seqnames =  genes.df$chr[i],
                                          ranges = IRanges(start=genes.df$start.point[i], 
                                                           end=genes.df$end.point[i])))
      
      pilup_params =  Rsamtools::PileupParam(max_depth = 10250,
                                             min_mapq = min.mapq,
                                             min_base_quality = 13, #################
                                             distinguish_nucleotides = T,
                                             ignore_query_Ns = T, 
                                             include_insertions =  T,
                                             distinguish_strands = T)
      
      bam.pileup = pileup(bamfile, 
                          scanBamParam=param,
                          pileupParam = pilup_params)
      
      # remove NA rows
      bam.pileup = bam.pileup[which(is.na(bam.pileup$nucleotide) == F),]
      
      # if there are any valid reads for the specified gene, then proceed...
      if (nrow(bam.pileup) > 1){
        print(paste(meta.data.df$replicate.name[k],": ",genes.df$Annotation[i], "  (",i,"out of",nrow(genes.df),")"))
        
        # order the raw pileup based on positions
        bam.pileup = bam.pileup[order(bam.pileup$pos), ]
        
        # specify the row index based on position (because in the raw pileup dataframe, positons are not unique)
        bam.pileup$row.indx = 1
        for (j in c(2:nrow(bam.pileup))){
          if (bam.pileup$pos[j] == bam.pileup$pos[j-1]){
            bam.pileup$row.indx[j] = bam.pileup$row.indx[j-1]
          } else {
            bam.pileup$row.indx[j] = bam.pileup$row.indx[j-1] + 1
          }
        }
        
        # create a new data frame, with unique positions (the polished pileup dataframe)
        unique.positions = unique(bam.pileup$pos)
        unique.positions = unique.positions[order(unique.positions)]
        
        pos.df = data.frame("chr"=as.character(bam.pileup$seqnames[1]),
                            "position"=unique.positions, 
                            "strand"=curr.strand,
                            "target.nucleotide"="",
                            "N_reads" = 0,
                            "A"=0,
                            "T"=0,
                            "C"=0,
                            "G"=0,
                            "del"=0,
                            "ins"=0)
        
        
        # loop over the raw pileup dataframe, and complete the polished dataframe
        for (j in c(1:nrow(bam.pileup))){
          if (bam.pileup$nucleotide[j] != "N" & bam.pileup$nucleotide[j] != "="  & bam.pileup$strand[j] == curr.strand){# sanity check
            cur.nucleotide = bam.pileup$nucleotide[j]
            if (bam.pileup$nucleotide[j] == "-") cur.nucleotide = "del"
            if (bam.pileup$nucleotide[j] == "+") cur.nucleotide = "ins"
            
            col.indx = which(names(pos.df) == cur.nucleotide)
            row.indx = bam.pileup$row.indx[j]
            
            # update the called base on polished pileup
            pos.df[row.indx, col.indx] = pos.df[row.indx, col.indx] + bam.pileup$count[j]
            
            # update the total number of reads base on polished pileup
            if (cur.nucleotide != 'ins'){
              pos.df$N_reads[row.indx] = pos.df$N_reads[row.indx] + bam.pileup$count[j]
            }
            
          }
        }
        # keep only rows with enough coverage
        pos.df = pos.df[which(pos.df$N_reads >= min.acceptable.reads), ]
        
        # now, use the hg38 fasta file, to determine the exact nucleotide on each position
        
        
        if (nrow(pos.df)>0){
          GR = 
            GRanges(seqnames = as.character(pos.df$chr[1]), 
                    ranges = IRanges(pos.df$position[1], end = pos.df$position[nrow(pos.df)]))
          
          refBase =getSeq(hg38.fa, GR)
          refBase = as.character(refBase[[1]])
          refBase = strsplit(refBase, "")[[1]]
          
          pos.indx = pos.df$position - pos.df$position[1] + 1
          refBase = refBase[pos.indx]
          pos.df$target.nucleotide = refBase
          
          write.csv(pos.df, file.name, row.names = F)
        }
        
      }
    }
    
    
  }
  
}

  }
  
}
