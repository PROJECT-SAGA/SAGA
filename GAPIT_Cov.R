#!/usr/local/R-3.3.1/bin/Rscript

source("http://www.bioconductor.org/biocLite.R")
#biocLite("multtest")
#install.packages("gplots")
#install.packages("LDheatmap")
#install.packages("genetics")
library(multtest)
library(gplots)
library(scatterplot3d)
library(LDheatmap)
library(genetics)
library(compiler)

#source("http://www.zzlab.net/GAPIT/previous/gapit_functions_20150515.txt")
#source("http://zzlab.net/GAPIT/emma.txt")
#source("http://bioinfo-web.mpl.ird.fr/sources_gapit/gapit_functions.txt") 
#source("http://bioinfo-web.mpl.ird.fr/sources_gapit/emma.txt")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
#source("http://bioinfo-web.mpl.ird.fr/sources_gapit_111016/gapit_functions.txt") 
#source("http://bioinfo-web.mpl.ird.fr/sources_gapit_111016/emma.txt") 
####
# Rajouter controle de saisie des 2 arguments
#setwd("C:/Users/cherif/Dropbox/GWAS/Données GBS Qat/LG19_DP10_88")
args <- commandArgs(TRUE)

myY  <- read.table(args[1], head = TRUE)
myG <- read.table(args[2] , head = FALSE)
myCV <- read.table(args[3], head= TRUE)
#Using ECMLM by Li and et. al. (BMC Biology, 2014)
myGAPIT <- GAPIT( Y=myY,  G=myG, CV=myCV,PCA.total=0,Model.selection = TRUE)
