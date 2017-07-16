module purge

module load compilers/gcc-6.2.0
module load compilers/openmpi-2.1.1_gcc-6.2
module load boost_1_56_0_gcc4_9_0

mkdir -p build_gcc ; cd build_gcc ; cmake .. -DMPI_C_COMPILER=/shared/software/compilers/openmpi-2.1.1/bin/mpicc -DMPI_CXX_COMPILER=/shared/software/compilers/openmpi-2.1.1/bin/mpicxx  -DCMAKE_LINKER=/shared/software/compilers/openmpi-2.1.1/bin/mpicxx  -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_C_COMPILER=/shared/software/compilers/gcc-6.2.0/bin/gcc -DCMAKE_CXX_COMPILER=/shared/software/compilers/gcc-6.2.0/bin/g++ ; cmake --build . --target install ; cd ..

