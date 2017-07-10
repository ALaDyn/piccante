#! /bin/bash

module purge

module load gnu
module load intel
module load intelmpi
module load boost

make marconi

