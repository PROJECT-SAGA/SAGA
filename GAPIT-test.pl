#!/usr/bin/perl 

use strict;
use warnings;
use Data::Dumper; 
use File::Copy;   
use Bio::SearchIO;
use Switch;

my %param = @ARGV;

#### UNTIL GAPIT
#my $geno_DP=10; # Minimum Depth of coverage (DP) at the sample level required for a minimum of number of sample
#my $geno_count=88; # Minimum sample's number required for the previous defined DP 
#my $snp_maf=0.05; # Minor allele frequency cut-off
#my $startrange=500; # sequence size range:  size before SNP
#my $endrange=500; # range after SNP
#my $pvalue=5E-2; # P-value
#my $evalue=1e-10; # blast evalue
#my $pid=97;       #blast: minimum hit percent_identity
#my $cov=80;        #blast coverage 
#my $maxwin=3000; # Annotation window:Upstream and Downstream max positions (pb) 
#my $minwin=1000; # Annotation window:Upstream and Downstream min positions (pb)
#my $gapit="/data2/projects/gbsflowering/perltoctoc/GAPIT_def.R \n"; # gérer le path par défaut?
my $gapit="DEF"; # Run default GAPIT method using ECMLM by Li and et. al. (BMC Biology, 2014) 
my $covar_file;  #Mendatory file tu run GAPIT with provided COVARIATE file
my $super="MLM"; #Only for GAPIT SUPER GAWS method 
#print $gapit;
#exit;
#
my $hapmap_file=$param{'-hap'}; 
#my $vcf_file=$param{'-vcf'};     # Input: VCF filename that will be read and analyzed
my $conf_file=$param{'-conf'};   # Config file containing settings    
#my $ref_file=$param{'-ref'};     # Reference sequence file (fasta format)
my $trait_file=$param{'-trait'}; # Trait file used by GAPIT
#my $gff_file=$param{'-gff'};     # Gff file
##my $ACP_file=$param{'-pca'};     # PCA file
##my $Kin_file=$param{'-kin'};     #Kinship file

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
         #case "dp"   { $geno_DP=$2 }
         #case "c"    { $geno_count=$2 }
         #case "maf"  { $snp_maf=$2 }
         #case "up"   { $startrange=$2 }
         #case "down" { $endrange=$2 }
         #case "pvalue" { $pvalue=$2}
         #case "evalue" {$evalue=$2}
         #case "pid"    {$pid=$2}
         #case "cov"    {$cov=$2}
         #case "maxwin" {$maxwin=$2}
         #case "minwin" {$minwin=$2}
         case "GAPIT" {$gapit=$2}
         case "cov_file" {$covar_file=$2}
         case "SUPER_Option" {$super=$2}
      }
   }
}
# Check if the files given by argument exist
#if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e  $ACP_file and -e $ref_file and -e $gff_file))
#if (not( -e $vcf_file and -e $conf_file and -e  $trait_file and -e $ref_file and -e $gff_file))
#{
#   print STDERR " #### $0 Error: The vcf file, config file, trait file, PCA file, gff file and/or reference files don't exist \n\n";
#   exit;
#}



#print "#### $0 Info: Used command options: \n
#      -dp      $geno_DP 
#      -c       $geno_count
#      -maf     $snp_maf 
#      -up      $startrange pb
#      -down    $endrange pb
#      -pvalue  $pvalue \n\n"
#      -GAPIT   $gapit";

print "#### $0 Info: Used command options: \n
       -GAPIT   $gapit";
if ($gapit eq "SUPER")
{
print "
       -SUPER_Option $super \n\n";
} else {next;}
exit;
## Check if all the files have been given by argument
##if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file) or not defined ($ACP_file) or not defined($ref_file) or not defined ($gff_file))
#if (not defined($vcf_file) or not defined ($conf_file) or not defined ($trait_file)  or not defined($ref_file) or not defined ($gff_file))
##Checking if input file is provided 
#{
#   print STDERR " #### $0 Error: The vcf file, config file, trait file,PCA file, gff file and/or reference files must be provided \n\n";
#   exit;
#}

############################################################################################
#                              Step 3 GWAS Analysis GAPIT                                  # 
############################################################################################
my $step3=localtime();
print "############     Step 3 GWAS Analysis GAPIT started at $step3     ############  \n\n";

#chdir $saga_gapit or die "#### $0 Error: Impossible to go into $saga_gapit directory $!";

#Gapit analysis : the script GAPIT.R is executed
print "#### $0 Info: GWAS Analysis is running on $trait_file and $hapmap_file \n\n";

print " $gapit \n";
#exit;

#my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT_CMLM_SUP.R $trait_file $hapmap_file $ACP_file $Kin_file";

#for ($gapit=$2)
#{
#    
#   if ($gapit=="COV")
#   {
#        print "#### $0 Info: GAPIT with $gapit is running \n\n";
#    
#        my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_Cov.R $trait_file $hapmap_file $covar_file";
#        print "$Rcmd \n";
#        exit;
#        system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");
#
#        print "#### $0 Info: GAPIT with $gapit has sucessfully run \n\n";
#   } elsif ($gapit=="SUPER"){
#        print "#### $0 Info: GAPIT with $gapit is running \n\n";
#
#        my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_SUPER.R $trait_file $hapmap_file $super";
#        print "$Rcmd \n";
#        exit;
#        system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");
#
#        print "#### $0 Info: GAPIT $gapit has sucessfully run \n\n";
#    } else {
#        print "#### $0 Info: Default script $gapit is running \n\n";
#
#        my $Rcmd="Rscript $gapit $trait_file $hapmap_file";
#        system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");
#    
#    print "#### $0 Info: Default script $gapit has sucessfully run \n\n";
#    } 
#}

#
if ($gapit eq "COV") {
print "#### $0 Info: GAPIT Using ECMLM by Li and et. al. (BMC Biology, 2014) with provided"." $gapit"."ARIATE file is running \n\n";
    
    my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_Cov.R $trait_file $hapmap_file $covar_file";
    
    if (not defined($covar_file))
    #Checking if input file is provided 
    {
        print STDERR " #### $0 Error: The Covariate file is missing from the config file and must be provided \n\n";
        exit;
    }
    # Get the absolute path of files given by argument
    $covar_file=`readlink -f $covar_file` or die "#### $0 Error: readlink -f $covar_file $!";
    ## DEBUG print "#### $covar_file\n";
    chomp $covar_file;
    
    print "#### $0 R command line: $Rcmd \n\n";
    system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");

print "#### $0 Info: GAPIT with $gapit has been sucessfully run \n\n";

} elsif ($gapit eq "SUPER"){
print "#### $0 Info: GAPIT using $gapit GWAS method is running \n\n";

    my $Rcmd="Rscript /data2/projects/gbsflowering/SAGAGIT/GAPIT_SUPER.R $trait_file $hapmap_file $super";
    print "#### $0 R command line: $Rcmd \n\n";
    system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");

print "#### $0 Info: GAPIT $gapit GWAS method has been sucessfully run \n\n";

} else {
print "#### $0 Info: Default GAPIT $gapit Using ECMLM by Li and et. al. (BMC Biology, 2014) is running \n\n";

    my $Rcmd="Rscript /data2/projects/gbsflowering/perltoctoc/GAPIT_def.R $trait_file $hapmap_file";
    print "#### $0 R command line: $Rcmd \n\n";
    
    system ($Rcmd) and die ("#### $0 Error: GAPIT step: $Rcmd");
    
print "#### $0 Info: Default script $gapit has been sucessfully run \n\n";
}

print "#### $0 Info: GWAS Analysis has been sucessfully done \n\n";
#else {print "toto \n";}
