#!/bin/sh



############      SGE CONFIGURATION      ###################

# Ecrit les erreur dans le fichier de sortie standard 

#$ -j y 


# Shell que l'on veut utiliser 

#$ -S /bin/bash 


# Email pour suivre l'execution 

#$ -M emira.cherif@ird.fr


# Type de massage que l'on reçoit par mail

#    -  (b) un message au demarrage

#    -  (e) a la fin

#    -  (a)  en cas d'abandon

#$ -m bea 


# Queue que l'on veut utiliser

#$ -q bigmem.q@node3.alineos.net


# Nom du job

#$ -N iss

############################################################

 

#path_to_dir="/data/projects/illumina-51";

#path_to_tmp="/scratch/tranchant-$JOB_ID/" 

 

###### Creation du repertoire temporaire sur noeud

#mkdir $path_to_tmp

#scp -rp 91.203.34.157:/$path_to_dir/* $path_to_tmp

#echo "tranfert donnees master -> noeud";


###### Execution du programme
cmd1="module load system/python/3.6.5"

echo "Commande executee : $cmd1";

$cmd1;

path_ref="/scratch/test_EC/Drosophila_melanogaster.BDGP6.22_inv_rearranged.fasta"
path_ouput="/scratch/test_EC"
prefix_file="hiseq_ref_BDGP6_10M"
nreads="10M"
cmd="/home/cherif/.local/bin/iss generate --cpus 4 --draft $path_ref --model hiseq --gc_bias --n_reads $nreads --output $path_ouput/${prefix_file}reads" 

echo "Commande executee : $cmd";

$cmd;
