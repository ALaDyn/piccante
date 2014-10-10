#! /bin/bash

if [ $# -lt 1 ]
then
 echo "usage: ./compile_on_juqueen.sh path_to_main.cpp"
 exit
fi

rm -f main-1.cpp
ln -s $1 main-1.cpp

EXE_NAME=($(basename $1 .cpp))


module purge
module load gsl/1.15_O3g
module load boost/1.47.0

make -f makefile.juqueen.xl all

mv piccante piccante.${EXE_NAME}

