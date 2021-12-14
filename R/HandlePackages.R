#### ====================================== ####
####     Install/load required packages     ####
#### ====================================== ####

#  Check required packages  ...
#  .. and install the required/not installed ones. 

list.of.packages <- 
  c("ggplot2", "optparse","BiocManager")

new.packages <- 
  list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

if (!("Biobase" %in% installed.packages()[,"Package"])){
  BiocManager::install("Biobase")
}

if (!("Rsamtools" %in% installed.packages()[,"Package"])){
  BiocManager::install("Rsamtools")
}

if (!("BSgenome" %in% installed.packages()[,"Package"])){
  BiocManager::install("BSgenome")
}
#  Quietly Load the required packages
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Biobase)))
suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(BSgenome)))