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

#source GAPIT v2
source ("http://zzlab.net/GAPIT/gapit_functions.txt")
source ("http://zzlab.net/GAPIT/emma.txt")
####
# Rajouter controle de saisie des 2 arguments
#setwd("C:/Users/cherif/Dropbox/GWAS/Données GBS Qat/LG19_DP10_88")
args <- commandArgs(TRUE)

# Trait arg
myY  <- read.table(args[1], head = TRUE)

# Genotype arg 
myG <- read.table(args[2] , head = FALSE)

# Covariate arg (PCA, Structure, etc.)
#myCV <- read.table(args[3], head= TRUE)

# Kinship arg 
#myKI <- read.table(args[4], head = FALSE)

#Using ECMLM by Li and et. al. (BMC Biology, 2014)
#myM <- GAPIT(Y=myY,G=myG,PCA.total=3,
#             Model.selection = TRUE, Geno.View.output=FALSE, Major.allele.zero = TRUE)

## Run GAPIT SUPER
myGAPIT_SUPER <- GAPIT(
  Y=myY[,c(1,2)],
  G=myG,
  #KI=myKI,
  #CV=myCV,
  PCA.total=3,
  sangwich.top=args[3], #options are GLM,MLM,CMLM, FaST and SUPER
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
  LD=0.1
)
