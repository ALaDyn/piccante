#!/bin/bash
#PBS -A IscrB_PULPICS                           
#PBS -l select=1:ncpus=68:mpiprocs=68:mem=108GB:mcdram=flat:numa=snc2
#PBS -l walltime=0:10:00
#PBS -o opic.out
#PBS -e epic.err

cd /marconi_scratch/userexternal/asgatton/TEST_1

module purge
module load intel
module load intelmpi
module load boost
module load env-knl
#module load gnu

mpirun ./piccante.exe > myoutput 
