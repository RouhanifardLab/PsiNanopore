#### ====================================== ####
####           Main: Plot Raw Signals       ####
#### ====================================== ####
#!/usr/bin/env Rscript
wdir = getwd( ) 
print(wdir)
source(paste(wdir,"/R/MakeOptions.R",sep = ""))
source(paste(wdir,"/R/return_LH_limits.R",sep = ""))
source(paste(wdir,"/R/InitParams.R",sep = ""))
wdir = getwd( ) 

#/Users/nabi_137/repositories/2021/SignalView/Package/AGA_hela_IVT_pos_177439607_hg38.txt 
#chr4 75814112 +
pdf(file = opt$out,width = 7.5,height = 4)
par(mfrow=c(1,1),lend=1); 
plot(-100,-100,xlim=c(-n.adj.nucleotides,n.adj.nucleotides),
     ylim=c(LH[1],LH[2]),main="",
     xlab=NA,ylab=NA,ann = F,axes = F)
box(lwd=1) 
for (i in c(1:n_samples)){
  signal.transparency=signal.transparency2;col = c(alpha("red",signal.transparency))
  file.name = paste(wdir,"/",opt$file2,sep = "")
  if (i == 1){
    file.name = opt$file
    col = c(alpha("black",signal.transparency))
  }

  sig.df = read.table(file.name,header = T)
  sig.df = sig.df[which(sig.df$position>(target.pos-n.adj.nucleotides) & 
                     sig.df$position<(target.pos+n.adj.nucleotides) &
                     sig.df$contig == target.chr),]
  
  for (read.indx in unique(sig.df$read_index)){
    sig.df.tmp = sig.df[which(sig.df$read_index == read.indx),]
    sig.df.tmp$base = ""
    positions = unique(sig.df.tmp$position)
    for (pos in positions){
      tot.n.signals2 = 0
      x.tmp = pos - target.pos - .5
      sig.df.positions.temp = sig.df.tmp[which(sig.df.tmp$position == pos),]
      
      base.temp = strsplit(sig.df.positions.temp$reference_kmer[1],"")[[1]][3]
      sig.df.positions.temp$base = base.temp
      tot.n.signals = 0
      raw.signals.lst = c()
      text(x = x.tmp +.7*n.adj.nucleotides/10,y = 1.1*LH[1],labels = base.temp,
           col=col.lst[which(base.lst==base.temp)])
      for (row in c(1:nrow(sig.df.positions.temp))){
        raw.signals = as.numeric(strsplit(sig.df.positions.temp$samples[row],",")[[1]])
        raw.signals.lst = c(raw.signals.lst,raw.signals)
        tot.n.signals = tot.n.signals + length(raw.signals)
      }
      dx = 1/tot.n.signals
      x.lst = x.tmp + c(1:tot.n.signals)/tot.n.signals
      lines(x.lst,raw.signals.lst,col=col)
      
      ### plot border lines for bases
      lines(x = c(x.tmp,x.tmp),y = c(-1,1000),col=alpha("darkblue",.1),lwd=.3)
    }
    lines(x = c(x.tmp+1,x.tmp+1),y = c(-1,1000),col=alpha("darkblue",.1),lwd=.3)
  }
  axis(2,cex.axis=1,tck=.02,at=c(1:100)*20,labels = c(1:100)*20,las=2,lwd=1);
  x.ax.labs = unique(sig.df$position)
  axis(1,cex.axis=1,tck=.02,at=x.ax.labs+0.5,labels = x.ax.labs,las=2,lwd=1);
}
dev.off()
print("Done!")
openPDF(opt$out)

