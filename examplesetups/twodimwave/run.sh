#!/bin/bash

cd ./run
#rm ./run/*

ln -s ../input/* .

../build/mitgcmuv > output.txt
#mpiexec -n 4 ../build/mitgcmuv > output.txt
