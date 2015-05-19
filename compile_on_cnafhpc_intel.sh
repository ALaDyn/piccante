
#! /bin/bash

#if [ $# -lt 1 ]
#then
# echo "usage: ./compile_at_cineca.sh path_to_main.cpp"
# exit
#fi

#rm -f main-1.cpp
#ln -s $1 main-1.cpp

#EXE_NAME=($(basename $1 .cpp))


module purge
module load compilers/ips-xe-2013-sp1
module load compilers/intel-mpi

make boost
#mv piccante piccante.${EXE_NAME}
mv piccante piccante.intel

