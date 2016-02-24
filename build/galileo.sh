#! /bin/bash

module purge
module load profile/advanced
module load intel/cs-xe-2015--binary
module load intelmpi/5.0.2--binary
module load boost/1.57.0--intel--cs-xe-2015--binary
module load gsl/1.16--intel--cs-xe-2015--binary

make galileo
