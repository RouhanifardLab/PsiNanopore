rm(list = ls());cat("\014")  
#CREATE A FOLDER NAMED Merged IN THE init_gene_pileup FOLDER BEFORE RUNNING EVERYTHING, MODIFY ROW 4 AND THE DATA FRAME IN LINE 49 
# Find the list of genes for which there is at least 1 pileup file
data.dir = "/CellLinesPVal/SH-SY5Y/init_gene_pileup/" 

replicates = list.files(data.dir)
replicates = replicates[-which(replicates == "Merged")]
target.genes  = c()

print('step2')


for (rep in replicates){
  temp.dir = paste0(data.dir, rep)
  target.genes = c(target.genes, list.files(temp.dir))
}

target.genes = unique(target.genes)

# 

for (i in c(1:length(target.genes))){
  gene = target.genes[i]
  if (i%%100 == 0) print(paste0(round(i/length(target.genes)*100, 2), "%"))
  out.file.name = paste0(data.dir,"Merged/",gene)
  if (file.exists(out.file.name) == F){
    # loop over all replicates for the current gene, find chr, strand, and list of all covered positions
    gene.chr = ""
    gene.strand = ""
    positions = c()
    for (rep in replicates){
      temp.file.name = paste0(data.dir, rep,"/",gene )
      if (file.exists(temp.file.name)&(file.size(temp.file.name)!= 0)){
        #if (file.exists(temp.file.name)){
        gene.pileup.df = read.csv(temp.file.name)
        gene.chr = gene.pileup.df$chr[1]
        gene.strand = gene.pileup.df$strand[1]
        positions = c(positions, gene.pileup.df$position)
      }
    }
    
    
    
    
    
    # if there are any covered positions
    if (length(positions)>0){
      
      
      # get unique positions covered for the current gene
      positions = unique(positions)
      # order positions
      positions = positions[order(positions)]
      
      
      #THE NAMES HERE NEED TO BE THE SAME AS THE NAME OF THE FOLDERS CONTAINING THE DIFFERENT REPLICATES IN init_gene_pileup
       merged.df = data.frame("chr"=gene.chr, "position"=positions, "strand"=gene.strand,"target.nucleotide"="",
                              "N_reads_rep1"=0, "A_rep1"=0, "T_rep1"=0, "C_rep1"=0, "G_rep1"=0 ,"del_rep1"=0, "ins_rep1"=0,
                              "N_reads_rep2"=0, "A_rep2"=0, "T_rep2"=0, "C_rep2"=0, "G_rep2"=0 ,"del_rep2"=0, "ins_rep2"=0,
                              "N_reads_rep3"=0, "A_rep3"=0, "T_rep3"=0, "C_rep3"=0, "G_rep3"=0 ,"del_rep3"=0, "ins_rep3"=0,
                              "N_reads_IVT"=0, "A_IVT"=0, "T_IVT"=0, "C_IVT"=0, "G_IVT"=0 ,"del_IVT"=0, "ins_IVT"=0)

       
      
      
      for (rep in replicates){
        temp.file.name = paste0(data.dir, rep,"/",gene )
        if (file.exists(temp.file.name)&(file.size(temp.file.name)!= 0)){
          #if (file.exists(temp.file.name)){
          gene.pileup.df = read.csv(temp.file.name)
          
          row.on.merged.df = which(merged.df$position %in% gene.pileup.df$position)
          merged.df$target.nucleotide[row.on.merged.df] = gene.pileup.df$target.nucleotide
          
          gene.pileup.df = gene.pileup.df[,c("N_reads","A","T","C","G","del","ins")]
          names(gene.pileup.df) = paste0(names(gene.pileup.df), "_",rep)
          
          col.on.merge.df = which(names(merged.df) %in% names(gene.pileup.df))
          
          merged.df[row.on.merged.df, col.on.merge.df] = gene.pileup.df
        }
        
      }
      write.csv(merged.df, out.file.name, row.names = F)
    }
    
  }
  
}  



