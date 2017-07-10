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
export JOB_5=subfolder_of_WDR_containing_job5_data
export JOB_6=subfolder_of_WDR_containing_job6_data
export JOB_7=subfolder_of_WDR_containing_job7_data
export JOB_8=subfolder_of_WDR_containing_job8_data
export JOB_9=subfolder_of_WDR_containing_job9_data
export JOB_10=subfolder_of_WDR_containing_job10_data
export JOB_11=subfolder_of_WDR_containing_job11_data
export JOB_12=subfolder_of_WDR_containing_job12_data
export JOB_13=subfolder_of_WDR_containing_job13_data
export JOB_14=subfolder_of_WDR_containing_job14_data
export JOB_15=subfolder_of_WDR_containing_job15_data
export JOB_16=subfolder_of_WDR_containing_job16_data



####################################################################
#### end of user defined variables
#### do not touch from here
export N_BGSIZE=64            ### Dimension of bg_size, the same set in the LoadLeveler keyword
                              ### Choose between 64 or 128

export N_SUBBLOCK=16          ### No. of sub-block you want.
                              ### Choose between 2, 4, 8, 16, 32, 64

export NPROC=64               ### No. of MPI tasks in each sub-block.

export RANK_PER_NODE=16       ### No. of MPI tasks in each node.
export WDR_1=$WDR/${JOB_1}
export WDR_2=$WDR/${JOB_2}
export WDR_3=$WDR/${JOB_3}
export WDR_4=$WDR/${JOB_4}
export WDR_5=$WDR/${JOB_5}
export WDR_6=$WDR/${JOB_6}
export WDR_7=$WDR/${JOB_7}
export WDR_8=$WDR/${JOB_8}
export WDR_9=$WDR/${JOB_9}
export WDR_10=$WDR/${JOB_10}
export WDR_11=$WDR/${JOB_11}
export WDR_12=$WDR/${JOB_12}
export WDR_13=$WDR/${JOB_13}
export WDR_14=$WDR/${JOB_14}
export WDR_15=$WDR/${JOB_15}
export WDR_16=$WDR/${JOB_16}
export EXE_1=${WDR_1}/${EXECUTABLE_NAME}
export EXE_2=${WDR_2}/${EXECUTABLE_NAME}
export EXE_3=${WDR_3}/${EXECUTABLE_NAME}
export EXE_4=${WDR_4}/${EXECUTABLE_NAME}
export EXE_5=${WDR_5}/${EXECUTABLE_NAME}
export EXE_6=${WDR_6}/${EXECUTABLE_NAME}
export EXE_7=${WDR_7}/${EXECUTABLE_NAME}
export EXE_8=${WDR_8}/${EXECUTABLE_NAME}
export EXE_9=${WDR_1}/${EXECUTABLE_NAME}
export EXE_10=${WDR_2}/${EXECUTABLE_NAME}
export EXE_11=${WDR_3}/${EXECUTABLE_NAME}
export EXE_12=${WDR_4}/${EXECUTABLE_NAME}
export EXE_13=${WDR_5}/${EXECUTABLE_NAME}
export EXE_14=${WDR_6}/${EXECUTABLE_NAME}
export EXE_15=${WDR_7}/${EXECUTABLE_NAME}
export EXE_16=${WDR_8}/${EXECUTABLE_NAME}
export EXECUTABLES="$EXE_1 $EXE_2 $EXE_3 $EXE_4 $EXE_5 $EXE_6 $EXE_7 $EXE_8 $EXE_9 $EXE_10 $EXE_11 $EXE_12 $EXE_13 $EXE_14 $EXE_15 $EXE_16"
export WORKINGDIRS="$WDR_1 $WDR_2 $WDR_3 $WDR_4 $WDR_5 $WDR_6 $WDR_7 $WDR_8 $WDR_9 $WDR_10 $WDR_11 $WDR_12 $WDR_13 $WDR_14 $WDR_15 $WDR_16"
export JOBNAMES="$JOB_1 $JOB_2 $JOB_3 $JOB_4 $JOB_5 $JOB_6 $JOB_7 $JOB_8 $JOB_9 $JOB_10 $JOB_11 $JOB_12 $JOB_13 $JOB_14 $JOB_15 $JOB_16"
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

