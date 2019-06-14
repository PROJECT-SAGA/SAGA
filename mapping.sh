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

#$ -N mapproc

############################################################

 

#path_to_dir="/data/projects/illumina-51";

#path_to_tmp="/scratch/tranchant-$JOB_ID/" 

 

###### Creation du repertoire temporaire sur noeud

#mkdir $path_to_tmp

#scp -rp 91.203.34.157:/$path_to_dir/* $path_to_tmp

#echo "tranfert donnees master -> noeud";


###### Execution du programme
path_ref="/scratch/test_EC/Drosophila_melanogaster.BDGP6.22.fa"
path_wdir="/scratch/test_EC/10M"
prefix_file="hiseq_inv_BDGP6_10Mreads"


#indexing
#cmd_Index="bwa index $path_ref" 

#echo "Commande executee : $cmd_Index";

#Mapping
rgstring='@RG\tID:foo\tSM:bar'
cmd_mem="bwa mem -R "$rgstring" $path_ref $path_wdir/${prefix_file}_R1.fastq $path_wdir/${prefix_file}_R2.fastq"
sub_var0="$path_wdir/${prefix_file}.sam"

echo "Commande executee : $cmd_mem > $sub_var0";

$cmd_mem > $sub_var0;

#SamToBam
cmd_view="samtools view -b $path_wdir/${prefix_file}.sam -o $path_wdir/${prefix_file}.bam"

echo "Commande executee : $cmd_view";

$cmd_view;

#BamOrder

cmd_order="samtools sort -l 0 -o $path_wdir/${prefix_file}.sorted.bam $path_wdir/${prefix_file}.bam"

echo "Commande executee : $cmd_order";

$cmd_order;

#bai

cmd_bai="samtools index -b $path_wdir/${prefix_file}.sorted.bam"

echo "Commande executee : $cmd_bai";

$cmd_bai;

#BamToBed

cmd_bamTobed="bamToBed -cigar -i $path_wdir/${prefix_file}.sorted.bam"
sub_var="$path_wdir/${prefix_file}.sorted.bed"

echo "Commande executee : $cmd_bamTobed > $sub_var";

$cmd_bamTobed > $sub_var;

#Samtools bedcov: bed Covrage
cmd_bedcov="samtools bedcov $path_wdir/${prefix_file}.sorted.bed $path_wdir/${prefix_file}.sorted.bam"
sub_var1="$path_wdir/${prefix_file}.sorted.bed.cov"

echo "Commande executee : $cmd_bedcov > $path_wdir/${prefix_file}.sorted.bed.cov";

$cmd_bedcov > $sub_var1;

# Extracting tags from bam file: NM MD MC AS
cmd_extract="samtools view $path_wdir/${prefix_file}.sorted.bam" 
cmd_extract1="grep -v '\*'"
cmd_extract2="cut -f 12,13,14,15"
sub_var2="$path_wdir/${prefix_file}.sorted.tag.bed"
echo "Commande executee : $cmd_extract | $cmd_extract1 | $cmd_extract2 > $sub_var2";

$cmd_extract | grep -v '\*' | $cmd_extract2 > $sub_var2;

# Fusion of bed files
cmd_fusion="paste  $path_wdir/${prefix_file}.sorted.bed.cov $path_wdir/${prefix_file}.sorted.tag.bed"
sub_var3="$path_wdir/${prefix_file}.sorted.bed.cov.tag"
echo "Commande executee : $cmd_fusion > $sub_var3";

$cmd_fusion > $sub_var3;

#Cleaning temp files
cmd_clean="rm $path_wdir/${prefix_file}.sorted.bed $path_wdir/${prefix_file}.sorted.bed.cov $path_wdir/${prefix_file}.sorted.tag.bed"

echo "Commande de nettoyage executee : $cmd_clean";

$cmd_clean;



##### Transfert des données du noeud vers master

#scp -rp $path_to_tmp/* 91.203.34.157:/$path_to_dir/

#echo "Transfert donnees node -> master";


#### Suppression du repertoire tmp noeud

#rm -rf $path_to_tmp

#echo "Suppression des donnees sur le noeud";
