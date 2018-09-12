# requirements

### perl version recommanded 5.24
> Data::Dumper 
> File::Copy   
> Bio::SearchIO ?
> Switch

### Softwares
> vcftools (v0.1.13)
> beagle (?)
> snpEff (1.4)
> gatk (>3.3)
> gapit (?)

### Installation

> create the directory `path/to/saga` where SAGA will be installed 

> Go into the SAGA directory and clone the git repository
`git clone https://github.com/PROJECT-SAGA/SAGA.git`

> Modify the file /path/to/SAGA/config.pm
This file will provide to SAGA the paths for all the softwares you will use in your workflows.

> Add the module config.pm  path to the PERL5LIB environment variable
export PERL5LIB=$PERL5LIB:/path/to/SAGA/

> Add the SAGA directory to the PATH environment variable
export PATH=$PATH:/path/to/SAGA

Note: you can add this to your ~/.bashrc to make it always available when you log in, as follows:

`
echo"""
export PERL5LIB=$PERL5LIB:/path/to/saga
export PATH=$PATH:/path/to/saga
""" >> ~/.bashrc
`

