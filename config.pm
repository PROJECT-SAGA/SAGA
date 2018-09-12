package config;

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
use Exporter;

our @ISA=qw(Exporter);

our @EXPORT=qw($java $vcftools $gatk $snpeff $RScripts ); #$gapit);


# Directory where SAGA code has been cloned from github
our $saga="putSagaPath";

# External bioinformatic tools path 
our $java = "java -Xmx12g -jar";

our $vcftools = "/usr/local/vcftools-0.1.13/bin/vcftools";                                                                # = "putVcfToolsPath";
our $beagle= "$java -Xmx20g -jar /usr/local/beagle-4.1/beagle.27Jul16.86a.jar ";    # = "$java putBeaglePath";
our $snpeff = "$java -jar /usr/local/snpEff-4.2/snpEff.jar";                        # = "putVctToolsPath";
our $gatk = "$java -Xmx20g -jar /usr/local/GenomeAnalysisTK.jar";                   # = "$java putGatkPath"; #/usr/local/gatk-3.6/GenomeAnalysisTK.jar 
our $RScripts = "/data2/projects/gbsflowering/SAGAGIT/Rscripts/";                   # = "putRScriptsPath"; 
#our $gapit = " ";

