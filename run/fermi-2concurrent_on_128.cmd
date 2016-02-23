#!/bin/bash

# @ account_no = my_cineca_computing_account
# @ job_name = sub_block.$(jobid)
# @ output = opic.$(jobid)
# @ error = epic.$(jobid)
# @ shell = /bin/bash
# @ environment = COPY_ALL
# @ wall_clock_limit = 00:30:00
# @ job_type = bluegene
# @ notification = always
# @ bg_size = 8
# @ notify_user = myemail@address
# @ queue

#################################################################
### USER SECTION
### Please modify only the following variables
#################################################################
export WDR=path/to/working/master/directory

export JOB_1=subfolder_of_WDR_containing_job1_data
export WDR_1=$WDR/${JOB_1}
export EXE_1=/path/to/executable/forthe/first/job

export JOB_2=subfolder_of_WDR_containing_job2_data
export WDR_2=$WDR/${JOB_2}
export EXE_2=/path/to/executable/forthe/second/job


####################################################################
#### end of user defined variables
#### do not touch from here
export N_BGSIZE=8           ### Dimension of bg_size, the same set in the LoadLeveler keyword

export N_SUBBLOCK=2           ### No. of sub-block you want.

export NPROC=64              ### No. of MPI tasks in each sub-block.

export RANK_PER_NODE=16       ### No. of MPI tasks in each node.
export EXECUTABLES="$EXE_1 $EXE_2"
export WORKINGDIRS="$WDR_1 $WDR_2"
export JOBNAMES="$JOB_1 $JOB_2"
n_exe () { echo $EXECUTABLES | awk "{print \$$1}"; }
n_dir () { echo $WORKINGDIRS | awk "{print \$$1}"; }
job_name () { echo $JOBNAMES | awk "{print \$$1}"; }

module load subblock
source ${SUBBLOCK_HOME}/bgsize_${N_BGSIZE}/npart_${N_SUBBLOCK}.txt

for i in `seq 1 $N_PART`;
do
  cd $(n_dir $i)
  echo $(n_exe $i)
  runjob --block $BLOCK_ALLOCATED --corner $(n_cor $i) --shape $SHAPE_SB  --np $NPROC --ranks-per-node $RANK_PER_NODE : $(n_exe $i) >> opic.$(job_name $i) 2>> epic.$(job_name $i) &
  cd ..
done
wait

