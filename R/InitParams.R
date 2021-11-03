#### ====================================== ####
####        Initialize required params      ####
#### ====================================== ####
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
n.adj.nucleotides = 7
signal.transparency = .15
signal.transparency2 = signal.transparency
base.lst = c("A","T","C","G")
col.lst = c("darkgreen","red","blue","orange")
target.pos = as.numeric(opt$position)
target.chr = opt$chromosome

if (opt$file!=""){
  n_samples = 1; 
  sig.df0 = read.table(paste(wdir,'/',opt$file,sep = ""),header = T)
  sig.df_tart = sig.df0[which(sig.df0$position>(target.pos-n.adj.nucleotides) & 
                                sig.df0$position<(target.pos+n.adj.nucleotides) &
                                sig.df0$contig == target.chr),]
  LH = return_LH_limits(sig.df_tart)
  if (opt$file2!=""){
    n_samples = 2
  }
}