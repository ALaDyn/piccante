#! /bin/bash

if [ $# -lt 1 ]
then
 echo "usage: ./compile_on_galileo.sh path_to_main.cpp"
 exit
fi

rm -f main-1.cpp
ln -s $1 main-1.cpp

EXE_NAME=($(basename $1 .cpp))


module purge
module load profile/advanced
module load intel intelmpi/5.0.2--binary
module load gsl/1.16--intel--cs-xe-2015--binary
module load boost/1.57.0--intel--cs-xe-2015--binarymodule

make -f makefile.galileo

make clean

