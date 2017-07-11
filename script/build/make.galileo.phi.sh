
#! /bin/bash

module purge
module load intel
module load intelmpi
module load mkl
source $INTEL_HOME/bin/compilervars.sh intel64

make galileo-phi

