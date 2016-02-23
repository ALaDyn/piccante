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
# @ bg_size = 64
# @ notify_user = myemail@address
# @ queue

#################################################################
### USER SECTION
### Please modify only the following variables
#################################################################
export EXECUTABLE_NAME=piccante
export WDR=path/to/working/master/directory
export JOB_1=subfolder_of_WDR_containing_job1_data
export JOB_2=subfolder_of_WDR_containing_job2_data
export JOB_3=subfolder_of_WDR_containing_job3_data
export JOB_4=subfolder_of_WDR_containing_job4_data


####################################################################
#### end of user defined variables
#### do not touch from here
export N_BGSIZE=64           ### Dimension of bg_size, the same set in the LoadLeveler keyword
                              ### Choose between 64 or 128

export N_SUBBLOCK=4           ### No. of sub-block you want.
                              ### Choose between 2, 4, 8, 16, 32, 64

export NPROC=256              ### No. of MPI tasks in each sub-block.

export RANK_PER_NODE=16       ### No. of MPI tasks in each node.
export WDR_1=$WDR/${JOB_1}
export WDR_2=$WDR/${JOB_2}
export WDR_3=$WDR/${JOB_3}
export WDR_4=$WDR/${JOB_4}
export EXE_1=${WDR_1}/${EXECUTABLE_NAME}
export EXE_2=${WDR_2}/${EXECUTABLE_NAME}
export EXE_3=${WDR_3}/${EXECUTABLE_NAME}
export EXE_4=${WDR_4}/${EXECUTABLE_NAME}
export EXECUTABLES="$EXE_1 $EXE_2 $EXE_3 $EXE_4"
export WORKINGDIRS="$WDR_1 $WDR_2 $WDR_3 $WDR_4"
export JOBNAMES="$JOB_1 $JOB_2 $JOB_3 $JOB_4"
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

