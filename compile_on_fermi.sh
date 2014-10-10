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
module load profile/advanced
module load bgq-xl
module load gsl/1.15--bgq-xl--1.0
module load boost/1.51.0--bgq-xl--1.0
module load scalasca/1.4.2

make -f makefile.cineca.xl all

mv piccante piccante.${EXE_NAME}

