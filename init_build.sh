pushd /home/export/online1/cpc152/finals/ljx/code/src/gromacs/mdlib/nbnxn_kernels/sw_subcore
make
popd
rm -rf /home/export/online1/cpc152/finals/ljx/code/build
mkdir /home/export/online1/cpc152/finals/ljx/code/build
pushd /home/export/online1/cpc152/finals/ljx/code/build
cmake .. -DCMAKE_CXX_COMPILER=mpiCC -DCMAKE_C_COMPILER=mpicc -DCMAKE_LINKER=mpiCC -DGMX_CYCLE_SUBCOUNTERS=ON -DCMAKE_EXE_LINKER_FLAGS=-lm\ -lm_slave\ -msimd -DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=on -DGMX_DOUBLE=OFF -DGMX_OPENMP=OFF -DGMX_BUILD_MDRUN_ONLY=ON -DBUILD_SHARED_LIBS=off -LH
make -j 20
popd
