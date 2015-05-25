module purge
module load profile/advanced
module load bgq-xl
module load gsl/1.15--bgq-xl--1.0
module load boost/1.51.0--bgq-xl--1.0
module load scalasca/1.4.2

make -f makefile.fermi.xl debug


