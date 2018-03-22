#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2016 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD, UMR DIADE
# Written by Emira Cherif & Christine Tranchant
#
###################################################################################################################################


# perl perltoctoc/QuatarPillar3-3-ec.pl -vcf SAGA_DATA1/WITHDUPLICAT.GATKSELECTVARIANT.QUAL200.DP1000.SNP.vcf -conf saga_EC.config -ref Bank/genome_tassel-v2.fasta -trait perltoctoc/GAPIT_Traits_Table.txt -gff Bank/ref_DPV01_scaffolds.gff3

# copy from QuatarPillar-bgl-ec-debug-posi.pl 

# perl /data2/projects/gbsflowering/perltoctoc/QuatarPillar-bgl-ec-debug-posi.pl -vcf /scratch/cherif/sort-all-merge-LG1-19.GATKSV.SNP.vcf -conf /data2/projects/gbsflowering/perltoctoc/saga_EC.config -ref /data2/projects/gbsflowering/Bank/sort-LG1-19.fasta -trait /data2/projects/gbsflowering/perltoctoc/flo_206_dat_trait.txt -gff /data2/projects/gbsflowering/Bank/ref_DPV01_scaffolds.gff3


use strict;       
use warnings;
use Data::Dumper; 
use File::Copy;   
use Bio::SearchIO;
use Switch;

my %param = @ARGV;

#### UNTIL GAPIT
my $geno_DP=10; # Minimum Depth of coverage (DP) at the sample level required for a minimum of number of sample
my $geno_count=88; # Minimum sample's number required for the previous defined DP 
my $snp_maf=0.05; # Minor allele frequency cut-off
my $startrange=500; # sequence size range:  size before SNP
my $endrange=500; # range after SNP
my $pvalue=5E-2; # P-value
my $evalue=1e-10; # blast evalue
my $pid=97;       #blast: minimum hit percent_identity
my $cov=80;        #blast coverage 
my $maxwin=3000; # Annotation window:Upstream and Downstream max positions (pb) 
my $minwin=1000; # Annotation window:Upstream and Downstream min positions (pb) 

my $vcf_file=$param{'-vcf'};     # Input: VCF filename that will be read and analyzed
my $conf_file=$param{'-conf'};   # Config file containing settings    
my $ref_file=$param{'-ref'};     # Reference sequence file (fasta format)
my $trait_file=$param{'-trait'}; # Trait file used by GAPIT
my $gff_file=$param{'-gff'};     # Gff file
#my $ACP_file=$param{'-pca'};     # PCA file
#my $Kin_file=$param{'-kin'};     #Kinship file 
my $start_run=localtime();

open (STDOUT,"> SAGA.log.txt") or die "#### $0 Error: Cannot create the file : SAGA.log.txt. $!";
print "\n#### $0 Info: SAGA's job has started at $start_run \n\n";
open (STDERR,"> SAGA.error.txt") or die "#### $0 Error: Cannot create the file : SAGA.error.txt. $!";
print STDERR " Job has started at $start_run \n\n";


#Open of the input VCF file
open (CONF,$conf_file) or die "#### $0 Error: impossible to open $conf_file \n";

while (my $line=<CONF>)
{
   if ($line =~/^-([^\s]+)\s([^#\s]+).*$/)
   {
      switch ($1)
      {      
         case "dp"   { $geno_DP=$2 }
         case "c"    { $geno_count=$2 }
         case "maf"  { $snp_maf=$2 }
         case "up"   { $startrange=$2 }
         case "down" { $endrange=$2 }
         case "pvalue" { $pvalue=$2}
         case "evalue" {$evalue=$2}
         case "pid"    {$pid=$2}
         case "cov"    {$cov=$2}
         case "maxwin" {$maxwin=$2}
         case "minwin" {$minwin=$2}
      }
   }
}



print "#### $0 Info: Used command options: \n
      -dp      $geno_DP 
      -c       $geno_count
      -maf     $snp_maf 
      -up      $startrange pb
      -down    $endrange pb
      -pvalue  $pvalue \n\n";


# Check if all the files have been given by argument
#if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file) or not defined ($ACP_file) or not defined($ref_file) or not defined ($gff_file))
if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file)  or not defined($ref_file) or not defined ($gff_file))
#Checking if input file is provided 
{
   print STDERR " #### $0 Error: The vcf file, config file, trait file,PCA file, gff file and/or reference files must be provided \n\n";
   exit;
}

# Check if the files given by argument exist
#if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e  $ACP_file and -e $ref_file and -e $gff_file))
if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e $ref_file and -e $gff_file))
{
   print STDERR " #### $0 Error: The vcf file, config file, trait file, PCA file, gff file and/or reference files don't exist \n\n";
   exit;
}

# Get the absolute path of files given by argument
$vcf_file=`readlink -f $vcf_file` or die "#### $0 Error: readlink -f $vcf_file $!";
## DEBUG print "#### $vcf_file\n";
chomp $vcf_file;

# Get the absolute path of files given by argument
$conf_file=`readlink -f $conf_file` or die "#### $0 Error: readlink -f $conf_file $!";
## DEBUG print "#### $conf_file\n";
chomp $conf_file;

# Get the path oh the ref file analysed
$ref_file=`readlink -f $ref_file` or die "#### $0 Error: readlink -f $ref_file $!";
## DEBUG print "#### $ref_file\n";
chomp $ref_file;

# Get the path oh the gff file analysed
$gff_file=`readlink -f $gff_file` or die "#### $0 Error: readlink -f $gff_file $!";
chomp $gff_file;

# Get the path oh the trait file analysed
$trait_file=`readlink -f $trait_file` or die "#### $0 Error: readlink -f $trait_file $!";
chomp $trait_file;
# Get the path oh the PCA file analysed
#$ACP_file=`readlink -f $ACP_file` or die "#### $0 Error: readlink -f $ACP_file $!";
#chomp $ACP_file;

# Get the path oh the Kinship file analysed
#$Kin_file=`readlink -f $Kin_file` or die "#### $0 Error: readlink -f $Kin_file $!";
#chomp $Kin_file;


# Get just the name of the directory without the filename
my $dir=`dirname $vcf_file`;
## DEBUG print "#### $dir\n";
chomp $dir;

# Go to vcf directory
chdir $dir or die "#### $0 Error: Impossible to go into vcf directory : $dir $!";

#Create working directories
#Create SAGA directory
my $saga_dir=$dir."/SAGA";
mkdir $saga_dir or die "#### $0 Error: Impossible to create $saga_dir directory into $dir $!";
chdir $saga_dir or die "#### $0 Error: Impossible to go into $saga_dir directory $!";
#Create FilteredVcf directory
my $saga_filter=$dir."/SAGA/FilteredVcf";
mkdir $saga_filter or die "#### $0 Error: Impossible to create $saga_filter directory into $saga_dir $!";

#Create GAPIT directory
my $saga_gapit="GAPIT";
mkdir $saga_gapit or die "#### $0 Error: Impossible to create $saga_gapit directory into $saga_dir $!";


############################################################################################
#                          Step 1 VCF cleaning and filtration                              # 
############################################################################################
my $step1=localtime();
print "############     Step 1 VCF cleaning  filtratring and phasing started at $step1     ############  \n\n";

print "#### $0 Info: Analyzed VCF File : $vcf_file \n\n";

#Creation of a new VCF file with SNPs filtration
my @tmpVcf= split /\//,$vcf_file;
my $vcf_out= $tmpVcf[-1];
my $vcf_prefix= $vcf_out;
$vcf_prefix=~ s/\.vcf//;
$vcf_out=$vcf_prefix."-DP".$geno_DP."-COUNT".$geno_count.".vcf";

#print "$vcf_out \n";
#exit;

open (VCFOUT,">".$saga_filter."/".$vcf_out) or die "#### $0 Error: Impossible to open filtered vcf file $saga_filter/$vcf_out \n";

print "#### $0 Info: Filtered VCF File : $vcf_out \n\n"; 

#Open of the input VCF file
open (VCF,$vcf_file) or die "#### $0 Error: impossible to open $vcf_file \n";

print "#### $0 Info: First SNPs filtering based on DP in progress \n\n"; 
while (my$line=<VCF>)
{     
      if ($line =~ /^#/) # Select and print VCF header
      {
         print VCFOUT $line;
      }
      elsif ($line=~ /\.\/\./) # Removing missing data
      {
         #print "****MIERDA*** \n"; 
         next;   
      }
      else
      {
         if ($line=~ /:\.:\.:/ ) {$line=~ s/:\.:\.:/:0:0:/g;} # Substitution of GT:AD:DP:GQ:PL missing data (0/0:.:.:3:0,3,38 is replaced by 0/0:0:0:3:0,3,38)
         #print "****MIERDA*** \n"; 
         #Parsing VCF FORMAT field by sample and counting the number of genotype with the choosen $geno_DP
         my @tab= split /\t/, $line; 
         my $count = 0;
         for (my $indice=9; $indice < @tab; $indice++) 
         {
            chomp $tab[$indice];
            my @genotype= split /:/, $tab[$indice];
            if ($genotype[2] >=$geno_DP) {$count++;} 
         }#print "****MERDE1:$geno_count - $count ***** \n";
         if ($count>=$geno_count) # The VCF file lines are printed if the chosen $geno_count condition is reached  
         {
            print VCFOUT $line;
            #print "****MERDE######## \n" ;
         }
      } 
}
close VCF;
close VCFOUT;

print "#### $0 Info: First SNPs filtering based on DP was successfully done \n\n";
print "#### $0 Info: $vcf_out has been sucessfully created \n\n";
#exit;


#Running vcftools
my $vcftool=localtime();
print "### Vcftoools has started at $vcftool \n\n";

print "#### $0 Info: SNPs filtering based on minor allele frequency ( $snp_maf ) in progress \n\n";
my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --remove-filtered FILTER-QUAL --remove-filtered FILTER-DP --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
#my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --remove-filtered LOW-QUAL --remove-filtered SnpCluster --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
#my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --remove-indv PdMUE124 --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
#my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --remove /data2/projects/gbsflowering/Bank/datremov.txt --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
#my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --keep /data2/projects/gbsflowering/Bank/nigerKeep.txt --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
#my $cmd="vcftools --vcf $saga_filter"."/"."$vcf_out --maf $snp_maf  --max-missing 1 --remove-filtered-geno-all --keep /data2/projects/gbsflowering/Bank/DKeep.txt --out $saga_filter"."/maf".$snp_maf."-".$vcf_out." --recode ";
print "$cmd \n\n";
system ($cmd) and die ("#### $0 Error: vcftools step: $cmd");


#Rename vcftools output
my $vcf_pref= $vcf_out;
$vcf_pref=~ s/\.vcf//;
move ( $saga_filter."/maf".$snp_maf."-".$vcf_out.".recode.vcf",$saga_filter."/maf".$snp_maf."-".$vcf_pref.".recode.vcf")or die("Impossible de copier le fichier");

my $recode= "maf".$snp_maf."-".$vcf_pref.".recode.vcf"; 
my $tab= "maf".$snp_maf."-".$vcf_pref;
my $beaglout1= "maf".$snp_maf."-".$vcf_pref.".Bgl";
my $beaglout2= "maf".$snp_maf."-".$vcf_pref.".Bgt";
print "\n#### $0 Info: SNPs filtering based on minor allele frequency ( $snp_maf ) was successfully done \n\n";
print "#### $0 Info: $recode has been sucessfully created \n\n";

#Running Beagle
my $bgl=localtime();
print "### Beagle has started at $bgl \n\n";

#Imputing step
print "#### $0 Info: Imputing genotypes with gl argument in progress \n\n";

my $beaglecmd="java -Xmx20g -jar /usr/local/beagle-4.1/beagle.27Jul16.86a.jar gl=$saga_filter/$recode out=$saga_filter/$beaglout1";
print "$beaglecmd \n\n";
system ($beaglecmd) and die ("#### $0 Error: Beagle gl step: $beaglecmd");

print "\n#### $0 Info: Imputing genotypes with gl argument was successfully done \n\n";
print "#### $0 Info: $beaglout1.vcf.gz has been sucessfully created \n\n";

#Phasing and IBD steps
#print "#### $0 Info: Phasing genotypes and IBD segment detection with gt and ibd arguments in progress \n\n";
#my $beaglecmd2="java -Xmx20g -jar /usr/local/beagle-4.1/beagle.27Jul16.86a.jar gt=$saga_filter/$beaglout1.vcf.gz ibd=true out=$saga_filter/$beaglout2";

print "#### $0 Info: Phasing genotypes with gt argument in progress \n\n";

my $beaglecmd2="java -Xmx20g -jar /usr/local/beagle-4.1/beagle.27Jul16.86a.jar gt=$saga_filter/$beaglout1.vcf.gz out=$saga_filter/$beaglout2";

print "$beaglecmd2 \n\n";
system ($beaglecmd2) and die ("#### $0 Error: Beagle gt step: $beaglecmd2");

print "\n#### $0 Info: Phasing genotypes with gt was successfully done \n\n";
#print "\n#### $0 Info: Phasing genotypes and IBD segment detection with gt and ibd arguments were successfully done \n\n";
print "#### $0 Info: $beaglout1.vcf.gz has been sucessfully created \n\n";
#print "#### $0 Info: $beaglout1.ibd has been sucessfully created \n\n";
#print "#### $0 Info: $beaglout1.hbd has been sucessfully created \n\n";

#Decompressing step
print "#### $0 Info: Decompressing $beaglout2.vcf.gz \n\n";
my $gzipcmd="gzip -d $saga_filter/$beaglout2.vcf.gz";
print "$gzipcmd \n\n";
system ($gzipcmd) and die ("#### $0 Error: Beagle gt step: $gzipcmd");

print "#### $0 Info: Decompression of $beaglout2.vcf.gz done \n\n";




############################################################################################
#                              Step 2 Format conversion VCFtoHapMap                        # 
############################################################################################
my $step2=localtime();
print "############     Step 2 Format conversion VCFtoHapMap started at $step2     ############  \n\n";

#Creation of hapmap file from vcftools output (vcf)using GATK VariantsToTable
print "#### $0 Info: GATK VariantsToTable on ( $beaglout2.vcf ) is running \n\n";

my $GATKcmd="java  -Xmx20g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -R $ref_file -T VariantsToTable -V $saga_filter/$beaglout2.vcf -F CHROM -F POS -F ID -F REF -F ALT  -o $saga_filter/".$tab.".table -GF GT -AMD";
print "$GATKcmd \n\n";
system ($GATKcmd) and die ("#### $0 Error: GATK VariantsToTable step: $GATKcmd");


#print "#### $0 Info: GATK VariantsToTable on ( $recode ) is running \n\n";
#
#my $GATKcmd="java  -Xmx20g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -R $ref_file -T VariantsToTable -V $saga_filter/$recode -F CHROM -F POS -F ID -F REF -F ALT  -o $saga_filter/".$tab.".table -GF GT -AMD";
#print "$GATKcmd \n\n";
#system ($GATKcmd) and die ("#### $0 Error: GATK VariantsToTable step: $GATKcmd");

#Opening of the GATK VariantsToTable output and creation of the hapmap file
#GATK VariantsToTable output 
#CHROM   POS     ID      REF     ALT     PdFID59.GT      PdFID60.GT      PdFID61.GT      PdFIT10.GT      PdFIT11.GT      PdFIT12.GT      PdFIT13.GT      PdFIT14.GT      PdFIT15.GT      
#LG1     24485   .       C       A       C/C                 C/A            A/A            C/C             C/C            C/A              C/C              C/C            C/A     
#hapmap output
#rs             alleles chrom   pos     strand  assembly        center  protLSID        assayLSID       panelLSID       QCcode  PdFID59 PdFID60 PdFID61 PdFIT10 PdFIT11 PdFIT12 PdFIT13 PdFIT14
#LG1_24485       C/A     1       24485   .       NA               NA      NA             NA               NA             NA      CC      CA      AA      CC      CC      CA      CC      CC      
 

open(TAB, $saga_filter."/".$tab.".table") or die "#### $0 Error: Impossible to open $saga_filter"."/".$tab.".table \n";

my $hapmap_file=$saga_filter."/".$tab.".hmp.txt";

open(HAPMAP,">".$hapmap_file) or die "#### $0 Error: Impossible to create $saga_filter."/".$hapmap_file \n";

my @Header= split /\t/,<TAB>;
print HAPMAP "rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode";

for (my $indice=5; $indice <@Header; $indice++)
{
   $Header[$indice]=~ s/\.GT//;
   print HAPMAP "\t".$Header[$indice];
}

while (my $line=<TAB>)
{
  chomp $line;
  my @line= split/\t/,$line;
  $line[0]=~/(\d+)$/;

  print HAPMAP $line[0]."_".$line[1]."\t".$line[3]."/".$line[4]."\t".$1."\t".$line[1]."\t"."."."\t"."NA"."\t"."NA"."\t"."NA"."\t"."NA"."\t"."NA"."\t"."NA";
  for (my $i=5; $i<@line; $i++)
  {
      $line[$i]=~ s/\///;
      print HAPMAP "\t".$line[$i];
  }
  print HAPMAP "\n";
}

close TAB;
close HAPMAP;

print "\n#### $0 Info: $tab.hmp.txt has been sucessfully created \n\n";

############################################################################################
#                              Step 3 GWAS Analysis GAPIT                                  # 
############################################################################################
my $step3=localtime();
print "############     Step 3 GWAS Analysis GAPIT started at $step3     ############  \n\n";

chdir $saga_gapit or die "#### $0 Error: Impossible to go into $saga_gapit directory $!";

#Gapit analysis : the script GAPIT.R is executed
print "#### $0 Info: GAPIT analysis is running on $trait_file and $ hapmap_file \n\n";

#my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT_CMLM_SUP.R $trait_file $hapmap_file $ACP_file $Kin_file";
my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT_def.R $trait_file $hapmap_file";

system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");

print "#### $0 Info: GAPIT has sucessfully run \n\n";


############################################################################################
#                              Step 4 Annotation                                           # 
############################################################################################
my $step4=localtime();
print "############     Step 4 Annotation started at $step4     ############  \n\n";

print "#### $0 Info: SAGA is opening GFF file \n\n";
# open the gff file to extract annotation information
open(GFF,"$gff_file") or die "#### $0 Error: Cannot open the gff file : $gff_file. $!";

my %gfftab;  
   
while (my $line=<GFF>)
{
   chomp $line;
   next if $line =~ /^#/; # next line if we have a comment line  
   
   my @tab= split /\t/,$line;
       
   if ($tab[2]=~ m/RNA{1}/)
      {
         my @tab2=split /;/,$tab[8];
         my @gene= grep( /^gene/, @tab2);
         my @prod= grep(/^product/, @tab2);
           
         $gfftab{$tab[0]}{$tab[3]}{'end position'}=$tab[4];
           
           if (exists ($gene[0]))
           {
              $gfftab{$tab[0]}{$tab[3]}{'annot'}=$gene[0]."\t";
           }
           else
           {
              $gfftab{$tab[0]}{$tab[3]}{'annot'}="\t";
           }
            
           if (exists ($prod[0]))
           { 
               $gfftab{$tab[0]}{$tab[3]}{'annot'}.=$prod[0]."\t";
           }
           else
           {
               $gfftab{$tab[0]}{$tab[3]}{'annot'}.="\t";    
           }
      }
}

#Rename GAPIT file
my @gapitFiles=`ls GAPIT*` or die ("#### $0 Error: Can't list GAPIT files \n \n");

foreach my $file (@gapitFiles)
{
   if ($file =~ /\.\./)
   {
    chomp $file;
    my $fileRenamed = $file;
    $fileRenamed =~ s/\.\./\./;
    move($file, $fileRenamed) or die ("#### $0 Error: Can't rename GAPIT file \n \n");
   }
}

my @gapitCsvFiles=`ls GAPIT*.GWAS.Results.csv` or die ("#### $0 Error: Can't list csv GAPIT files \n \n");

foreach my $gwas (@gapitCsvFiles)
{
   
   my %gwastab;
   my @snp;
   my @posi;
   my @seq;
   
   my @tmp=split /\./, $gwas;
   
   #Create GAPIT.TRAIT sub-directory on GAPIT directory
   print "#### $0 Info: Starting $tmp[1] annotation step \n\n";
   my $saga_Trait_Annot=$tmp[1]."_Annot";
      mkdir "$saga_Trait_Annot" or die "#### $0 Error: Impossible to create $saga_Trait_Annot directory into saga_gapit $!";

  
   # Create the BLAST directory where the files generated by blastn will be saved
   my $blastpath=$saga_Trait_Annot."/BLAST";
   system ("mkdir $blastpath") and die "#### $0 Error: mkdir $blastpath. $!";
   chomp $blastpath;
   
   # open the file generated by GAPIT and parsed
   open (GWAS,"$gwas") or die "#### $0 Error: Cannot open the file : $gwas. $!";
   
   # open the file that will contain the SNP selected after filtering and blast annotation,
   open (SNPOUT,">".$saga_Trait_Annot."/".$tmp[1]."_snp.select.csv") or die "#### $0 Error: Cannot create the file : ".$tmp[1]."_snp.select.csv. $!";
   print SNPOUT "SNP_ID\tLG\tSNP_Position\tP.value" ."\t". "maf" ."\t". "Rsquare.of.Model.without.SNP\tRsquare.of.Model.with.SNP\tFDR Adjusted P-values\tExtract_start\tExtract_end\tFilter\tAln_scaf\tNum Hits\tAln_start" ."\t". "Aln_end\tAln_cov\tPid\tEvalue\tGene_start" ."\t". "Gene_stop" ."\t". "Gene_ID" ."\t". "Gene_product" ."\t". "Hit_position". "\n";
   
   # open the file that will contain the SNP non selecting after filtering step
   open (discSNP,">".$saga_Trait_Annot."/".$tmp[1]."_discard.snp.csv") or die "#### $0 Error: Cannot create the file : ".$tmp[1]."_discard.snp.csv. $!";
   
 
    ###############################################
    #      Step 4a Annotation: SNPs selection     # 
    ###############################################
   #my $step4a=localtime();
   ## DEBUG print "##########     Step 4a Annotation: SNPs selection started at $step4a     ##########  \n\n";
   
   ## DEBUGG my $count=1;
   <GWAS>;
   while (my $line=<GWAS>)
   {
      chomp $line;
      my @tab= split /,/,$line;
      ## DEBUG print "\n----------- $count : $tab[0] $tab[3]\n";
      if ($tab[3]<=$pvalue) #Select SNPs according the defined association pvalue and extract corresponding features 
      {
         $gwastab{$tab[0]}{'LG'}="LG".$tab[1];
         $gwastab{$tab[0]}{'posi'}=$tab[2];
         chomp $tab[8];
         $gwastab{$tab[0]}{'file'}=$tab[3]."\t".$tab[4]."\t".$tab[6]."\t".$tab[7]."\t".$tab[8];
         ## DEBUG print "---------- SELECT $tab[0]\n";
         ## DEBUG $count++;
      }
      else
      {
         print discSNP $line."\n";
         ## DEBUG print "-------- DISCARD $tab[0]";
         ## DEBUG $count++;
      }
      
   }
   
   #print Dumper(\%gwastab);



   ###############################################
   #     Step 4b Annotation: SNPs annotation     # 
   ###############################################
   #my $step4b=localtime();
   ## DEBUG print print "##########     Step 4b Annotation: SNPs annotation started at $step4b     ##########  \n\n";
   
   print "#### $0 Info: Sequence extraction and blast for $tmp[1]_trait annotation are running  \n\n";
   
   #Extract sequences surrounding signficant associated SNPs
   foreach my $snpID(keys %gwastab)
   {
      my $start=$gwastab{$snpID}{'posi'}-$startrange;                                                                           # Defines sequence start position 
      if ($start <1)
      {
         $start=1;
      }
      my $end=$gwastab{$snpID}{'posi'}+$endrange;                                                                                # Defines sequence end position
       
      my $cmd="blastdbcmd -db $ref_file -entry '$gwastab{$snpID}{'LG'}'  -range $start-$end -out $blastpath/$snpID.fasta";       # Extracts sequences from reference file
      system ($cmd) and die ("#### $0 Error: blastdbcmd: $cmd");
      my $cmd2="blastn -db /data2/projects/gbsflowering/Bank/42345_ref_DPV01_chrUn.fa -query $blastpath/$snpID.fasta -out $blastpath/$snpID.blastn -evalue $evalue "; #-outfmt 6 "; # Blasts extracted sequences 
      #my $cmd2="megablast -d /data2/projects/gbsflowering/Bank/42345_ref_DPV01_chrUn.fa -i $blastpath/$snpID.fasta -o $blastpath/$snpID.blastn -m 9 -D 2"; #-evalue $evalue"; # Blasts extracted sequences 
      
      #megablast -d /data2/projects/gbsflowering/Bank/42345_ref_DPV01_chrUn.fa -i /data2/projects/gbsflowering/SAGAV/SAGA/GAPIT/Flo_Annot/BLAST/LG3_4854345.fasta -o /data2/projects/gbsflowering/SAGAV/SAGA/GAPIT/Flo_Annot/BLAST/LG3_4854345.megablast -m 9
      system ($cmd2) and die ("blastn: $cmd2");
      
      # Starting of  Bio::SearchIO module : http://www.bioperl.org/wiki/HOWTO:SearchIO 
      # Extracts information from blast files 
      my $in = new Bio::SearchIO(-format => 'blast', 
                                 -file   => "$blastpath/$snpID.blastn");
  
      while( my $result = $in->next_result )
      {
         my $num_hits=$result->num_hits();
         my $count_hits=0;         
         #print "GRRRRRRRRRRRRRRRRRR result = in->next_result\n";
         if ($num_hits==0)
         {
            #print  "No blast result";
            print SNPOUT "$snpID\t$gwastab{$snpID}{'LG'}\t$gwastab{$snpID}{'posi'}\t$gwastab{$snpID}{'file'}\tNA\tNA\tNO HIT\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
         }  
         else   
         {
               #print "GRRRRRRRRRRRRRRRRRR num_hits!=0\n";

               ## $result is a Bio::Search::Result::ResultI compliant object
               ## $hit is a Bio::Search::Hit::HitI compliant object
               while (my $hit = $result->next_hit) 
               {
                  my $subj_accesid= $hit->name;
                  my ($ref1,$scaf)=split /\|/,$subj_accesid; # Gets the ref id from the gff : NM_xxxx.1
                  #print "###  $snpID $subj_accesid \n";
                  ## $hsp is a Bio::Search::HSP::HSPI compliant object
                  my $hsp = $hit->next_hsp;
                  my $aln_length=0;
                  my $aln_pid=0; 
                  my $aln_evalue=0;
                  my $aln_start=0;
                  my $aln_end=0;
                  my $query_length=0;
                  my $aln_cov=0;
                  if (defined $hsp)
                  {
                     $aln_length=$hsp->length('total');
                     $aln_pid=$hsp->percent_identity; 
                     $aln_evalue=$hsp->evalue;
                     $aln_start=$hsp->start('hit');
                     $aln_end=$hsp->end('hit');
                     $query_length=$result->query_length;
                     $aln_cov=sprintf("%.2f",($aln_length*100)/$query_length); #Calculates the percentage of query alignment
                     $count_hits++;
                     #print "--- ".$scaf."\t".$count_hits."\t".$aln_start."\t".$aln_end."\t".$aln_cov."\t$aln_pid\t$aln_evalue\n";
                  }
                  else {#print "--- aie".$scaf."\n";
                        goto PASS; }
                  
                  if ( $aln_pid >= $pid and $aln_cov >= 80) #Defines the minimum $aln_pid value to extract required information from the blastn file
                  {
                     $count_hits++;
                     $gwastab{$snpID}{'blast'}="PASS\t".$scaf."\t".$count_hits."\t".$aln_start."\t".$aln_end."\t".$aln_cov."\t$aln_pid\t$aln_evalue";
                     print SNPOUT "$snpID\t$gwastab{$snpID}{'LG'}\t$gwastab{$snpID}{'posi'}\t$gwastab{$snpID}{'file'}\t$start\t$end\t$gwastab{$snpID}{'blast'}\t";
                     delete $gwastab{$snpID}{'blast'};
                  }   
                  elsif ( $aln_pid >= ($pid-12) and $aln_cov >= 60) ########### SEN RAPPELER DE CA, TIENS!!!!!!!!!       ................. Tu as déjà oublié hein !!!!
                  {
                     $count_hits++;
                     $gwastab{$snpID}{'blast'}="LOWCOV\t".$scaf."\t".$count_hits."\t".$aln_start."\t".$aln_end."\t".$aln_cov."\t$aln_pid\t$aln_evalue";
                     print SNPOUT "$snpID\t$gwastab{$snpID}{'LG'}\t$gwastab{$snpID}{'posi'}\t$gwastab{$snpID}{'file'}\t$start\t$end\t$gwastab{$snpID}{'blast'}\t";
                     delete $gwastab{$snpID}{'blast'};
                  } 
                  else
                  {
                     $gwastab{$snpID}{'blast'}="BADCOV\t".$scaf."\tNA\t".$aln_start."\t".$aln_end."\t".$aln_cov."\t$aln_pid";
                     print SNPOUT "$snpID\t$gwastab{$snpID}{'LG'}\t$gwastab{$snpID}{'posi'}\t$gwastab{$snpID}{'file'}\t$start\t$end\t$gwastab{$snpID}{'blast'}\tNA\tNA\tNA\tNA\tNA\tNA\n";
                     delete $gwastab{$snpID}{'blast'};
                  }
                  
                  #print "GRRRRRRRRRRRRRRRRRR $snpID\t$gwastab{$snpID}{'LG'}\t$gwastab{$snpID}{'posi'}\t$gwastab{$snpID}{'file'}\t$count_hits\t$aln_start\t$aln_end\t$gwastab{$snpID}{'blast'}\n";                             
                  if (( $aln_pid >= $pid and $aln_cov >= 80) or ( $aln_pid >= ($pid-12) and $aln_cov >= 60))    ########### SEN RAPPELER DE CA, TIENS!!!!!!!!!       ................. Tu as déjà oublié hein !!!!
                  {
                     if (not exists $gfftab{$scaf})
                     {
                        $gwastab{$snpID}{'gff'}=" NA \t NA \t NA \t NA \t NO GFF ANNOTATION";
                     }
                     else
                     {
                        #Filling the annotation outputted file 
                        foreach my $gene_start (sort { $a <=> $b }keys(%{$gfftab{$scaf}}) )# Sorting { $a <=> $b} according to numeric order
                        {
                            my $gene_stop= $gfftab{$scaf}{$gene_start}{'end position'};
                            my $annotation= $gfftab{$scaf}{$gene_start}{'annot'};
                            #Setting up the annotation positions 
                            my $upmax=1;
                            $upmax=$gene_start-$maxwin if ($gene_start-$maxwin > 0); #Upstream position max
                            my $upmin=1;
                            $upmin=$gene_start-$minwin if ($gene_start-$minwin > 0); #Upstream position min
                            my $downmax=$gene_stop+$maxwin;#Downstream position max
                            my $downmin=$gene_stop+$minwin;#Downstream position min
                            # SNP annotation
                            if ($gene_start<=$aln_start and $gene_stop>=$aln_end) 
                            {
                               $gwastab{$snpID}{'gff'}=$gene_start."\t".$gene_stop."\t".$annotation."Gene"; # Gene   
                            }
                            elsif($gene_start>$aln_start and $gene_start<=$aln_end) # 
                            {
                               $gwastab{$snpID}{'gff'}=$gene_start."\t".$gene_stop."\t".$annotation."Up _Gene"; # Up _Gene
                            }
                            elsif($gene_stop>=$aln_start and $gene_stop<$aln_end)
                            {
                               $gwastab{$snpID}{'gff'}=$gene_start."\t".$gene_stop."\t".$annotation."Down _Gene"; # Down _Gene 
                            }    
                            elsif ($upmax<=$aln_start and $upmin>$aln_end) 
                            {  
                               $gwastab{$snpID}{'gff'}=$gene_start."\t".$gene_stop."\t".$annotation."Env_promot_Gene"; # Env_promot_Gene  
                            }
                            elsif ($downmax>=$aln_end and $downmin<$aln_start) 
                            { 
                              $gwastab{$snpID}{'gff'}=$gene_start."\t".$gene_stop."\t".$annotation."Env_utr_Gene"; # Env_utr_Gene   
                            }
                                 
                            $gwastab{$snpID}{'gff'}=" NA \t NA \t NA \t NA \t NON CODING" if (not exists ($gwastab{$snpID}{'gff'}));
                      
                        } #foreach my $gene_start (sort { $a <=> $b }keys(%{$gfftab{$scaf}}) )# Sorting { $a <=> $b} according to numeric order
                        
                     }  #fin  else if (not exists $gfftab{$scaf})
                     print SNPOUT "$gwastab{$snpID}{'gff'}\n";
                  } # END if (( $aln_pid >= $pid and $aln_cov >= 80) or ( $aln_pid >= ($pid-2) and $aln_cov >= 60))
                  
                  
                  
                  
                  ######### TO COMMENT AFTER RESOLVINBG BUG                 
                  #last;
                  
                  
                  
                  
                  
               } # END while (my $hit = $result->next_hit) 
         } # END ELSE if ($num_hits==0) 
      } # END while( my $result = $in->next_result )
      
      PASS :
 
   }# END foreach my $snpID(keys %gwastab) 
   
   print "#### $0 Info: $tmp[1]_snp.select.csv has been sucessfully created \n\n";
   print "#### $0 Info: $tmp[1]_discard.snp.csv has been sucessfully created \n\n"; 
   
   close GWAS;
   close SNPOUT;
   close discSNP;
   close GFF;
   
} # END foreach my $gwas (@gapitCsvFiles)
my $end_run=localtime ();
print "#### $0 Info: Job ended: $end_run \n\n";
print STDERR " #### $0 Info: Job ended: $end_run \n\n";

close STDOUT;
close STDERR;
close CONF;

############ A DECOMMENTER QUAND DEBUGGUER
#####
system ("chgrp gbsflowering $saga_dir -R") and die "#### $0 Error: chgrp $saga_dir. $!";
#####
system ("chmod 770 $saga_dir -R") and die "#### $0 Error: chmod 770 $saga_dir. $!";
