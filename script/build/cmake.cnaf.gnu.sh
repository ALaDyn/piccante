module purge

module load compilers/gcc-4.9.0
module load compilers/openmpi-1.8.1_gcc-4.9.0_with_cuda6.5
module load boost_1_56_0_gcc4_9_0

mkdir -p build_gcc ; cd build_gcc ; cmake .. -DMPI_C_COMPILER=/shared/software/compilers/openmpi-1.8.4_gcc4.9.0_cuda6.5/bin/mpicc -DMPI_CXX_COMPILER=/shared/software/compilers/openmpi-1.8.4_gcc4.9.0_cuda6.5/bin/mpicxx  -DCMAKE_LINKER=/shared/software/compilers/gcc-4.9.0/bin/gcc  -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_C_COMPILER=/shared/software/compilers/gcc-4.9.0/bin/gcc -DCMAKE_CXX_COMPILER=/shared/software/compilers/gcc-4.9.0/bin/g++ ; cmake --build . --target install ; cd ..

