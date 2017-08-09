#!/bin/bash

nprun=1

# Run
cd ./run
rm -r ./*
ln -s ../input/* .

if [ $nprun -eq 1 ]
then
    ../build/mitgcmuv > output.txt
else
    mpiexec -n $nprun ../build/mitgcmuv > output.txt
fi
