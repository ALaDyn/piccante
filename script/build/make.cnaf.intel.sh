
#! /bin/bash


module purge
module load compilers/gcc-4.8.2
module load compilers/intel-parallel-studio-2016

make cnaf-intel

