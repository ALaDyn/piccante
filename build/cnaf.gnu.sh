#! /bin/bash

module purge
module load compilers/gcc-4.9.0
module load boost_1_56_0_gcc4_9_0
module load compilers/openmpi-1.8.1_gcc-4.8.2

make cnaf
mv piccante.exe piccante.gcc
