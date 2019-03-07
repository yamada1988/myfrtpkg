#!/bin/bash

num=32
for i in 100
do

cd MELT_$i
mkdir aveERmod_01
mkdir aveERmod_01/refs
rm -f aveERmod_01/refs/weight*

cp ../script/ave_cor.f90 script
gfortran script/ave_cor.f90 -o ave.exe
./ave.exe

cd ../
for j in `seq 1 1`
do
  j1=`expr $num \* $j - $num + 1`
  j2=`expr $j1 + $num - 1`
  python script/make_avengr.py $i $i  $j1 $j2 $j
done



done
