#!/bin/bash
#PBS -A my_cineca_computing_account
#PBS -l walltime=4:0:00
#PBS -l select=4:ncpus=16:mpiprocs=16:mem=120GB
#PBS -j eo

##PBS -o opic.txt   #commented to have output in real-time via redirection and not all at the end of the job
##PBS -e epic.txt


cd ${PBS_O_WORKDIR}

module purge

module load gnu
module load intel
module load intelmpi
module load boost

mpirun ./piccante >> ./opic.txt



