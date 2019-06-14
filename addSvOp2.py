#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
import os, sys, pathlib, argparse, re
import pandas as pd
# =============================================================================
#Verifying Steps
# =============================================================================
#verifying Argument
parser=argparse.ArgumentParser()
parser.add_argument("bed" ,help="Take a bed like file as argument/ file path needed if it's outside your working directory")
args=parser.parse_args()

# Captures path and verifies if it's an actual file 
from pathlib import Path
filename=Path(args.bed)

# Verifies if the file exists
if not filename.exists():
    print("Oops, file doesn't exist! \n")
    exit()
else:
    print(filename.name, "\n") 

if not Path.is_file(filename):
# debug   print(filename.name) 
    print("Oops, this is not a file! \n")
    exit()
else:
    print("file successfully loaded")
# =============================================================================    

print("Parsing File step\n")
# =============================================================================
#Opening file, extracting CIGAR pattern and remaining tag
# =============================================================================
fd = open(filename,"r")
fo = open("/scratch/test_EC/deletion/20M/hiseq_rearg_20M_BDGP6.sorted.bed.cov.tag.var.csv","w")
fo.write('chr	Start	End	Read_Id	Bed_Score	Sens	CIGAR	CovPerReadRegion	NM_dist2ref	MD_missmatchposi	MC_CIGARmate	AlignScore	Variant\n')


for line in fd :
    newline=line.rstrip()

    Novariant=re.search("(\+|-)\s+(\d{1,}M)\s+",newline)
    delet = re.search("(\d{1,}(H|S|M){1}\d{1,}(M|S|H))",newline)
    if Novariant:
        fo.write(newline+'\t'+'NoVariant'+'\n')
    elif delet:

        fo.write(newline+'\t'+'Deletion'+'\n')
    else:
        fo.write(newline+'\t'+'SNV'+'\n')

 
# # =============================================================================
fd.close()
fo.close()

# # =============================================================================
# # =============================================================================
