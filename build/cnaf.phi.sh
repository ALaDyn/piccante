
#! /bin/bash


module purge
module load compilers/gcc-4.8.2
module load compilers/intel-parallel-studio-2016
source /shared/software/compilers/intel/compilers_and_libraries_2016/linux/bin/compilervars.sh intel64
source /shared/software/compilers/intel/compilers_and_libraries_2016/linux/mpi/intel64/bin/mpivars.sh

make cnaf-phi

