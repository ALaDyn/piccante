#! /bin/bash

if [ $# -lt 1 ]
then
 echo "usage: ./compile_at_cineca.sh path_to_main.cpp"
 exit
fi

rm -f main-1.cpp
ln -s $1 main-1.cpp

EXE_NAME=($(basename $1 .cpp))


module purge
module load compilers/gcc-4.8.2
module load compilers/openmpi-1.8.1_gcc-4.8.2

make all

mv piccante piccante.${EXE_NAME}
