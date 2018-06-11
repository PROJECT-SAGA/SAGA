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
#setwd("C:/Users/cherif/Dropbox/GWAS/DonnéesGBS Qat/LG19_DP10_88")

args <- commandArgs(TRUE)
#Trait arg
myY  <- read.table(args[1], head = TRUE)

#Genotype arg
myG <- read.table(args[2] , head = FALSE)
#myCV <- read.table(args[3], head= TRUE)

#Using ECMLM by Li and et. al. (BMC Biology, 2014)
myGAPIT <- GAPIT(Y=myY,G=myG,PCA.total=3,Model.selection = TRUE)
