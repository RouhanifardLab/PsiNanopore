#### ====================================== ####
####        Find the min/max of signals     ####
#### ====================================== ####
return_LH_limits = function(sig.df){
  max_sig = -Inf
  min_sig = Inf
  
  for (row in c(1:nrow(sig.df))){
    raw.signals = as.numeric(strsplit(sig.df$samples[row],",")[[1]])
    if (max(raw.signals, na.rm = T)>max_sig){max_sig = max(raw.signals, na.rm = T)}
    if (min(raw.signals, na.rm = T)<min_sig){min_sig = min(raw.signals, na.rm = T)}
  }
  limits.lst = c(min_sig, max_sig)
  return (limits.lst)
}


