
#! /bin/bash


module purge
module load compilers/gcc-4.8.2
module load compilers/intel-parallel-studio-2016
source /shared/software/compilers/intel/bin/compilervars.sh intel64

make cnaf-phi

