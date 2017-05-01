#!/bin/bash

## This script will install the tools required for the STRetch pipeline.
## It will fetched each tool from the web and placed into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required. Note that R must be installed manually
##

mkdir -p tools/bin
cd tools

#a list of which programs need to be installed
commands="bpipe conda bwa samtools"

#installation method
function bpipe_install {
   wget -O bpipe-0.9.9.2.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.2/bpipe-0.9.9.2.tar.gz
   tar -zxvf bpipe-0.9.9.2.tar.gz ; rm bpipe-0.9.9.2.tar.gz
   ln -s $PWD/bpipe-0.9.9.2/bin/* $PWD/bin/
}

# Installs miniconda, Python 3 + required packages, BedTools and goleft
# (and any other dependancies listed in environment.yml)
function conda_install {
    wget -O miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash miniconda.sh -b -p $PWD/miniconda
    rm miniconda.sh
    $PWD/miniconda/bin/conda env create -f environment.yml
    ln -s $PWD/miniconda/envs/STR/bin/* $PWD/bin/
#    source activate STR
}

function bwa_install {
    wget --no-check-certificate https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2
    tar -jxvf bwa-0.7.15.tar.bz2
    rm bwa-0.7.15.tar.bz2
    make prefix=$PWD install -C bwa-0.7.15/
}

function samtools_install {
   wget --no-check-certificate https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2
   tar -jxvf samtools-1.3.1.tar.bz2
   rm samtools-1.3.1.tar.bz2
   make prefix=$PWD install -C samtools-1.3.1/
}

echo "// Path to tools used by the pipeline" > ../tools.groovy

for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which $PWD/bin/$c 2>/dev/null`
    fi
    echo "$c=\"$c_path\"" >> ../tools.groovy
done

#finally check that R is installed
R_path=`which R 2>/dev/null`
if [ -z $R_path ] ; then
    echo "R not found!"
    echo "Please go to http://www.r-project.org/ and follow the installation instructions."
    echo "Please also install the required R packages:"
    echo "install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'))"
fi
echo "R=\"$R_path\"" >> ../tools.groovy

#loop through commands to check they are all installed
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo -n "WARNING: $c could not be found!!!! "
	echo "You will need to download and install $c manually, then add its path to tools.groovy"
	Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running STRetch."
    else
        echo "$c looks like it has been installed"
    fi
done
echo "**********************************************************"
echo $Final_message
