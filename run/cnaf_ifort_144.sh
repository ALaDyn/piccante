#!/bin/bash


NUMERO_TOTALE_CORE_DA_USARE=144

NOME_ESEGUIBILE="./piccante"
stderr_file=epic.txt
stdout_file=opic.txt

job=job_test_ifort.cmd
job_name=impi
queue=hpc_inf

###########################
rm -f $job
touch $job
chmod 755 $job

touch ${stderr_file}
touch ${stdout_file}

echo "#BSUB -J ${job_name}" > $job # Job name
echo "#BSUB -o %J.out" >> $job # Job standard output
echo "#BSUB -e %J.err" >> $job # Job standard error
echo "#BSUB -q ${queue}" >> $job
echo "#BSUB -n ${NUMERO_TOTALE_CORE_DA_USARE}" >> $job
echo "module load compilers/gcc-4.8.2" >> $job
echo "module load compilers/intel-parallel-studio-2016" >> $job
echo "module load boost_1_56_0_gcc4_9_0" >> $job
echo "/shared/software/compilers/impi/intel64/bin/mpirun -np ${NUMERO_TOTALE_CORE_DA_USARE} -genv PSM_SHAREDCONTEXTS_MAX 8 -genv I_MPI_FABRICS shm:tmi ${NOME_ESEGUIBILE} >> ${stdout_file} 2>> ${stderr_file}" >> $job


echo "Lanciare il job con il seguente comando: "
echo "bsub < $job"
