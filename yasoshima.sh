#!/bin/bash

num=32
for i in `seq -f %03g 90 -10 0`
do

cd MELT_$i
mkdir aveERmod_01
mkdir aveERmod_01/refs
rm -f aveERmod_01/refs/weight*

cp ../script/yasoshima.f90 script
gfortran script/yasoshima.f90 -o ave.exe
./ave.exe

cp ../MELT_000/ERmod/parameters_fe aveERmod_01/
cd aveERmod_01
slvfe > slvfe.tt

cd ../../

done
