#!/bin/bash

# @ account_no = my_cineca_computing_account
# @ shell = /bin/bash
# @ job_type = serial
# @ job_name = building_piccante.$(jobid)
# @ output = building_piccante.out.$(jobid)
# @ error = building_piccante.err.$(jobid)
# @ wall_clock_limit = 0:30:00
# @ class = serial
# @ notification = always
# @ notify_user = myemail@address
# @ queue


module purge
module load profile/advanced
module load bgq-xl
module load boost/1.51.0--bgq-xl--1.0

cd /path/to/piccante/

make fermi-perf

