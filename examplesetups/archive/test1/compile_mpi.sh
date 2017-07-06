#!/bin/bash

gcmdir='/data5/glwagner/gcm/MITgcm'
optfile="$gcmdir/tools/build_options/sverdrup_glwagner"

cd ./build

$gcmdir/tools/genmake2 \
    -optfile=$optfile \
    -mpi \
    -mods=../code/ \
    -rootdir=$gcmdir/

make depend
make
