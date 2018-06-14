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



use strict;       
use warnings;
use Data::Dumper; 
use File::Copy;   
use Bio::SearchIO;
use Switch;

my %param = @ARGV;

##########
########## CD : Verifier les parametres  a conserver ou pas in fine
##########
my $phasing=0;
my $geno_DP=10; # Minimum Depth of coverage (DP) at the sample level required for a minimum of number of sample
my $geno_count=88; # Minimum sample's number required for the previous defined DP 
my $snp_maf=0.05; # Minor allele frequency cut-off
my $pvalue=5E-2; # P-value
my $maxwin=10000; # Annotation window:Upstream and Downstream max positions (pb) 
my $minwin=2000; # Annotation window:Upstream and Downstream min positions (pb) 

my $conf_file=$param{'-conf'};   # Config file containing settings    
my $ref_file=$param{'-ref'};     # Reference sequence file (fasta format)
my $vcf_file=$param{'-vcf'};     # Input: VCF filename that will be read and analyzed
my $trait_file=$param{'-trait'}; # Trait file used by 

my $gff_file= defined($param{'-gff'}) ? $param{'-gff'} : 0;     # Gff file
my $snpeff_db= defined ($param{'-snpeff'})? $param{'-snpeff'} : 0;

$phasing=1 if defined ($param{'-phasing'});     # Gff file???????

#my $ACP_file=$param{'-pca'};     # PCA file
#my $Kin_file=$param{'-kin'};     #Kinship file

my $start_run=localtime();


### DEBUG
print "---- $gff_file ----- $snpeff_db -------";

############################################
# Open log files where script is executed
############################################
# standard output
open (STDOUT,"> SAGA.log.txt") or die "#### $0 Error: Cannot create the file : SAGA.log.txt. $!";
print STDOUT "\n#### $0 Info: SAGA's job has started at $start_run \n\n";
# error output
open (STDERR,"> SAGA.error.txt") or die "#### $0 Error: Cannot create the file : SAGA.error.txt. $!";
print STDERR " Job has started at $start_run \n\n";





############################################
# Check if all the files have been given by argument
############################################

#if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file) or not defined ($ACP_file) or not defined($ref_file) or not defined ($gff_file))
if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file)  or not defined($ref_file))
#Checking if input file is provided 
{
   print STDERR " #### $0 Error: The vcf file, config file, trait file,PCA file, gff file and/or reference files must be provided \n\n";
   exit;
}

# Check if the files given by argument exist
#if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e  $ACP_file and -e $ref_file and -e $gff_file))
if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e $ref_file))
{
   print STDERR " #### $0 Error: The vcf file, config file, trait file, PCA file, gff file and/or reference files don't exist \n\n";
   exit;
}
else
{
   # Get the absolute path of files given by argument
   $vcf_file=`readlink -f $vcf_file` or die "#### $0 Error: readlink -f $vcf_file $!";
   chomp $vcf_file;

   $conf_file=`readlink -f $conf_file` or die "#### $0 Error: readlink -f $conf_file $!";
   chomp $conf_file;

   $ref_file=`readlink -f $ref_file` or die "#### $0 Error: readlink -f $ref_file $!";
   chomp $ref_file;

   $trait_file=`readlink -f $trait_file` or die "#### $0 Error: readlink -f $trait_file $!";
   chomp $trait_file;
}


# Check if the script could perform annotation 
if ($gff_file eq "0")
{
   # If no gff and no snpeff database providen
   if ($snpeff_db eq "0")
   {
      print STDERR " #### $0 Error: The gff file or a snpeff database must be provided \n\n";
      exit;
   }
   #if snpeff database given, let snpeff check that the database exists   

}
elsif (not(-e $gff_file))
{
   print STDERR " #### $0 Error: The gff file doenn't exist \n\n";
   exit;
}
else
{
   # Get the path oh the gff file analysed
   $gff_file=`readlink -f $gff_file` or die "#### $0 Error: readlink -f $gff_file $!";
   chomp $gff_file;
}



############################################
#Get the file path
############################################

##########
########## CD 2/06 : Si on decommente apres, tester l existence du fichier
##########
# Get the path oh the PCA file analysed
#$ACP_file=`readlink -f $ACP_file` or die "#### $0 Error: readlink -f $ACP_file $!";
#chomp $ACP_file;

# Get the path oh the Kinship file analysed
#$Kin_file=`readlink -f $Kin_file` or die "#### $0 Error: readlink -f $Kin_file $!";
#chomp $Kin_file;


############################################
#Open of the input conf file to get additional parameters
############################################
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
         case "pvalue" { $pvalue=$2}
      }
   }
}



##########
########## CD : Maj tous les parametres  a afficher une fois fixŽs
##########
print "#### $0 Info: Used command options: \n
      -dp      $geno_DP 
      -c       $geno_count
      -maf     $snp_maf 
      -pvalue  $pvalue \n\n";






############################################
#Move into the working directory and preparing it
############################################

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

#Creation of a new VCF file with SNPs cleaning
my @tmpVcf= split /\//,$vcf_file;
my $vcf_out= $tmpVcf[-1];
my $vcf_prefix= $vcf_out;
$vcf_prefix=~ s/\.vcf//;
$vcf_out=$vcf_prefix."-DP".$geno_DP."-COUNT".$geno_count.".vcf";

#Open the vcf files (input & output)
print "#### $0 Info: Filtered VCF File : $vcf_out \n\n"; 
open (VCFOUT,">".$saga_filter."/".$vcf_out) or die "#### $0 Error: Impossible to open filtered vcf file $saga_filter/$vcf_out \n";

print "#### $0 Info: First SNPs filtering based on DP in progress \n\n"; 
open (VCF,$vcf_file) or die "#### $0 Error: impossible to open $vcf_file \n";
while (my$line=<VCF>)
{     
      if ($line =~ /^#/) # Select and print VCF header
      {
         print VCFOUT $line;
      }
      elsif ($line=~ /\.\/\./) # Removing missing data
      {
         next;   
      }
      else
      {
         $line=~ s/:\.:\.:/:0:0:/g if ($line=~ /:\.:\.:/ ); # Substitution of GT:AD:DP:GQ:PL missing data (0/0:.:.:3:0,3,38 is replaced by 0/0:0:0:3:0,3,38)
  
         #Parsing VCF FORMAT field by sample and counting the number of genotype with the choosen $geno_DP
         my @tab= split /\t/, $line; 
         my $count = 0;
         for (my $indice=9; $indice < @tab; $indice++) 
         {
            chomp $tab[$indice];
            my @genotype= split /:/, $tab[$indice];
            $count++ if ($genotype[2] >=$geno_DP); 
         }
         
         print VCFOUT $line if ($count>=$geno_count); # The VCF file lines are printed if the chosen $geno_count condition is reached  
      } 
}
close VCF;
close VCFOUT;

print "#### $0 Info: First SNPs filtering based on DP was successfully done \n\n";
print "#### $0 Info: $vcf_out has been sucessfully created \n\n";



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
print "################################ ".$saga_filter."/maf".$snp_maf."-".$vcf_out.".recode.vcf   ---- ".$saga_filter."/maf".$snp_maf."-".$vcf_pref.".recode.vcf";

move ($saga_filter."/maf".$snp_maf."-".$vcf_out.".recode.vcf",$saga_filter."/maf".$snp_maf."-".$vcf_pref.".recode.vcf")or die("Impossible de copier le fichier ");

my $recode= "maf".$snp_maf."-".$vcf_pref.".recode.vcf"; 
my $tab= "maf".$snp_maf."-".$vcf_pref;
my $beaglout1= "maf".$snp_maf."-".$vcf_pref.".Bgl";
my $beaglout2= "maf".$snp_maf."-".$vcf_pref.".Bgt";
print "\n#### $0 Info: SNPs filtering based on minor allele frequency ( $snp_maf ) was successfully done \n\n";
print "#### $0 Info: $recode has been sucessfully created \n\n";


if ($phasing)
{
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
   $beaglout2.=".vcf";
}
else
{
   $beaglout2=$recode;
}

############################################################################################
#                              Step 2 Format conversion VCFtoHapMap                        # 
############################################################################################

if (not ($snpeff_db eq "0"))
{

   print "############     Step 2a SNP annotation with SNPeff started at ".localtime()."    ############  \n\n";

   #Launch SNPeff
   print "#### $0 Info: SNPeff on ( $beaglout2 ) is running \n\n";

   my $snpeffCmd="java -jar /usr/local/snpEff-4.2/snpEff.jar $snpeff_db $saga_filter"."/"."$beaglout2 > $saga_filter"."/".$beaglout2.".SNPeff.vcf";
   print "$snpeffCmd \n\n";
   system ($snpeffCmd) and die ("#### $0 Error: SNPeff step: $snpeffCmd");

   $beaglout2=$beaglout2.".SNPeff.vcf";
}



############################################################################################
#                              Step 2 Format conversion VCFtoHapMap                        # 
############################################################################################
my $step2=localtime();
print "############     Step 3 Format conversion VCFtoHapMap started at $step2     ############  \n\n";

#Creation of hapmap file from vcftools output (vcf)using GATK VariantsToTable
print "#### $0 Info: GATK VariantsToTable on ( $beaglout2 ) is running \n\n";

my $GATKcmd="java  -Xmx20g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -R $ref_file -T VariantsToTable -V $saga_filter/$beaglout2 -F CHROM -F POS -F ID -F REF -F ALT  -o $saga_filter/".$tab.".table -GF GT -AMD";
print "$GATKcmd \n\n";
system ($GATKcmd) and die ("#### $0 Error: GATK VariantsToTable step: $GATKcmd");

###
##Opening of the GATK VariantsToTable output and creation of the hapmap file
###

# GATK VariantsToTable output 
# CHROM   POS     ID      REF     ALT     PdFID59.GT      PdFID60.GT      PdFID61.GT      PdFIT10.GT      PdFIT11.GT      PdFIT12.GT      PdFIT13.GT      PdFIT14.GT      PdFIT15.GT      
# LG1     24485   .       C       A       C/C                 C/A            A/A            C/C             C/C            C/A              C/C              C/C            C/A     

# hapmap output
# rs             alleles chrom   pos     strand  assembly        center  protLSID        assayLSID       panelLSID       QCcode  PdFID59 PdFID60 PdFID61 PdFIT10 PdFIT11 PdFIT12 PdFIT13 PdFIT14
# LG1_24485       C/A     1       24485   .       NA               NA      NA             NA               NA             NA      CC      CA      AA      CC      CC      CA      CC      CC      
 

my $hapmap_file=$saga_filter."/".$tab.".hmp.txt";
open(HAPMAP,">".$hapmap_file) or die "#### $0 Error: Impossible to create $saga_filter."/".$hapmap_file \n";
open(TAB, $saga_filter."/".$tab.".table") or die "#### $0 Error: Impossible to open $saga_filter"."/".$tab.".table \n";

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
#my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT_def.R $trait_file $hapmap_file";
#my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_SUPER.R $trait_file $hapmap_file"; 
#my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_Cov.R $trait_file $hapmap_file /data2/projects/gbsflowering/RICE/GAPIT.PCA.txt";  # CD
my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT.R $trait_file $hapmap_file /data2/projects/gbsflowering/RICE/GAPIT.PCA.txt";

system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");

print "#### $0 Info: GAPIT has sucessfully run \n\n";


############################################################################################
#                              Step 4 Annotation                                           # 
############################################################################################
my $step4=localtime();
print "############     Step 4 Annotation started at $step4     ############  \n\n";


######################### VCF PARSING & SNPEFF EXTRACTING INTO %ANNOT 
my %annot;
my %gff;

if (not ($snpeff_db eq "0"))
{
   open(VCF,"<",$saga_filter."/".$beaglout2) or die "#### $0 Error: Cannot open the file : $saga_filter."/".$beaglout2. $!";
   while (my $line = <VCF>)
   {
         chomp $line;
         next if ( $line =~/^#/ );
         my @line = split /\t/, $line;
         $annot{$line[0]."_".$line[1]}=$line[7];
   }
   #print "!!!!!!!!!!!!!! %annot\n";
   #print Dumper(\%annot);
   #print "!!!!!!!!!!!!!! FIN %annot\n";
}
#####################################################################
else
{
   open(GFF,"<",$gff_file) or die "#### $0 Error: Cannot open the file : $gff_file. $!";
   while (my $line = <GFF>)
   {
         chomp $line;
         next if ( $line =~/^#/ );
         my @line = split /\t/, $line;
         if ($line[2] eq "gene")
         {
            $gff{$line[0]}{$line[3]}{'stop'}=$line[4];
            $gff{$line[0]}{$line[3]}{'annot'}=$line[8];
         }
   }
   #print "!!!!!!!!!!!!!! %gff\n";
   #print Dumper(\%gff);
   #print "!!!!!!!!!!!!!! FIN %gff\n";
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

   # open the file generated by GAPIT and parsed
   open (GWAS,"$gwas") or die "#### $0 Error: Cannot open the file : $gwas. $!";
   
   # open the file that will contain the SNP selected after filtering and blast annotation,
   open (SNPOUT,">".$saga_Trait_Annot."/".$tmp[1]."_snp.select.csv") or die "#### $0 Error: Cannot create the file : ".$tmp[1]."_snp.select.csv. $!";
   #print SNPOUT "LG\tSNP_Position\tSNP_Position\tSNP_ID\tP.value" ."\t". "maf" ."\t". "Rsquare.of.Model.without.SNP\tRsquare.of.Model.with.SNP\tFDR Adjusted P-values\n";
   
   # open the file that will contain the SNP non selecting after filtering step
   open (discSNP,">".$saga_Trait_Annot."/".$tmp[1]."_discard.snp.csv") or die "#### $0 Error: Cannot create the file : ".$tmp[1]."_discard.snp.csv. $!";
   
 
    ###############################################
    #      Step 4a Annotation: SNPs selection     # 
    ###############################################
   #my $step4a=localtime();
   ## DEBUG print "##########     Step 4a Annotation: SNPs selection started at $step4a     ##########  \n\n";
   
   <GWAS>;
   while (my $line=<GWAS>)
   {
      chomp $line;
      my @tab= split /,/,$line;
      
      $gwastab{$tab[0]}{'LG'}="LG".$tab[1];
      $gwastab{$tab[0]}{'posi'}=$tab[2];
      chomp $tab[8];
      $gwastab{$tab[0]}{'file'}=$tab[3]."\t".$tab[4]."\t".$tab[6]."\t".$tab[7]."\t".$tab[8];

   
      #### DEBUG PB CHR GAPITWT
      my $annot;
      my $lastGene=0; 
      
      if (not ($snpeff_db eq "0") and not defined $annot)
      {
          #print "!!!!!!!!!!!!!! SNPEFF GWAS $tab[0] ---- \n";  
          $annot=$annot{$tab[0]};
      }
      elsif (not ($gff_file eq "0"))
      {
         #print "!!!!!!!!!!!!!! GFF GWAS1 $tab[1]\n";  
         foreach my $start (sort {$a <=> $b} keys % { $gff{"Chr".$tab[1]} })
         {
            #print "!!!!!!!!!!!!!! GFF GWAS1 $tab[1] - $start ---- ".$gff{"Chr".$tab[1]}{$start}{'stop'}." ----\n";  
            if ($tab[2] > $start and $tab[2] < $gff{"Chr$tab[1]"}{$start}{'stop'})
            {
               $annot="Gene $start-".$gff{"Chr$tab[1]"}{$start}{'stop'}." (".$gff{"Chr$tab[1]"}{$start}{'annot'};
               last;
            }
            elsif ($lastGene !=0 and $tab[2] < $start and $tab[2] > $lastGene) # add window
            {
               $annot="Up Gene $start-".$gff{"Chr$tab[1]"}{$start}{'stop'}." (".$gff{"Chr$tab[1]"}{$start}{'annot'};
            }
            $lastGene=$gff{"Chr$tab[1]"}{$start}{'stop'};
         }
      }
      else { $annot = "na"; } #print}
         
      if ($tab[3]<=$pvalue) #Select SNPs according the defined association pvalue and extract corresponding features 
      {
         
         print SNPOUT "$tab[0]\tChr$tab[1]\t$tab[2]\t$tab[3]\t$tab[4]\t$tab[6]\t$tab[7]\t$tab[8]\t ---- $annot ----\n";
      }
      else
      {
         print discSNP $line."\t ---- $annot -----\n"; ########## A Enlever $annot qd ok
      }
      
   }
   
   #print Dumper(\%gwastab);



   ###############################################
   #     Step 4b Annotation: SNPs annotation     # 
   ###############################################
   #my $step4b=localtime();
   ## DEBUG print print "##########     Step 4b Annotation: SNPs annotation started at $step4b     ##########  \n\n";

   #my $bedCmd="bedtools window -a $gff_file -b $saga_Trait_Annot"."/".$tmp[1]."_snp.select.csv -w 8000 > $saga_Trait_Annot"."/".$tmp[1].".gene";
   #system ($bedCmd) and die ("#### $0 Error: bedtools window step: $bedCmd");

   #print "#### $0 Info: bedtools window has sucessfully run \n\n";
   print "#### $0 Info: $tmp[1]_snp.select.csv has been sucessfully created \n\n";
   print "#### $0 Info: $tmp[1]_discard.snp.csv has been sucessfully created \n\n"; 
   
   close GWAS;
   close SNPOUT;
   close discSNP;
   
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
