export PAHT=${PATH}:/shared/software/project/aladyn/cmake-3.8.2-Linux-x86_64/bin/

module purge

module load compilers/gcc-6.2.0
module load compilers/intel-parallel-studio-2017
module load boost_1_56_0_gcc4_9_0 

mkdir -p build ; cd build ; cmake .. -DCMAKE_LINKER=/shared/software/compilers/intel_2017/compilers_and_libraries_2017.0.098/linux/bin/intel64/icpc  -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_C_COMPILER=/shared/software/compilers/intel_2017/compilers_and_libraries_2017.0.098/linux/bin/intel64/icc -DCMAKE_CXX_COMPILER=/shared/software/compilers/intel_2017/compilers_and_libraries_2017.0.098/linux/bin/intel64/icpc ; cmake --build . --target install ; cd ..

