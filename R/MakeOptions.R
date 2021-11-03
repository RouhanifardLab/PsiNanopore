#### ====================================== ####
####    Make/define the main code options   ####
#### ====================================== ####
#wdir = getwd( ) 
source(paste(wdir,"/R/HandlePackages.R",sep = ""))

option_list = list(
  make_option(c("-f", "--file"), type="character", default="", 
              help="signal file name", metavar="character"),
  make_option(c("-g", "--file2"), type="character", default="", 
              help="second signal file name", metavar="character"),
  make_option(c("-s", "--strand"), type="character", default=NULL, 
              help="strand", metavar="character"),
  make_option(c("-p", "--position"), type="character", default=NULL, 
              help="position", metavar="character"),
  make_option(c("-c", "--chromosome"), type="character", default=NULL, 
              help="chromosome", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.pdf", 
              help="output file name of signal plot", metavar="character")
); 