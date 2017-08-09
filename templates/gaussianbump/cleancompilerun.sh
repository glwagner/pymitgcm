#!/bin/bash

# Clean up build and run directory
rm -r ./build/*
rm -r ./run/*

# Compile
gcmdir='/data5/glwagner/gcm/MITgcm'
optfile="$gcmdir/tools/build_options/sverdrup_glwagner"

cd ./build

$gcmdir/tools/genmake2 \
    -optfile=$optfile \
    -mods=../code/ \
    -enable=mnc \
    -rootdir=$gcmdir/

make depend
make -j8

# Run
cd ../
./run.sh
