#!/bin/bash

# Number of processors for making and running
npmake=2
nprun=1

gcmdir="/Users/glwagner/Software/MITgcm"
optfile="$gcmdir/tools/build_options/neve"

# Clean up build and run directory
rm -r ./build
rm -r ./run
mkdir build
mkdir run

# Compile
cd ./build

if [ $nprun -eq 1 ]
then # do not compile with mpi:
    $gcmdir/tools/genmake2 \
        -optfile=$optfile \
        -mods=../code/ \
        -enable=mnc \
        -rootdir=$gcmdir/
else # compile with mpi:
    $gcmdir/tools/genmake2 \
        -optfile=$optfile \
        -enable=mnc \
        -mpi \
        -mods=../code/ \
        -rootdir=$gcmdir/
fi


make depend
make -j$npmake

# Run
cd ../run
ln -s ../input/* .

if [ $nprun -eq 1 ]
then
    ../build/mitgcmuv > output.txt
else
    mpiexec -n $nprun ../build/mitgcmuv > output.txt
fi

