#!/bin/bash


cd output/simulation/random_ini

for i in copy_*
do
cd $i
cp ../../../../src/convert.sh .
cp ../../../../src/fortran2bd.py .
./convert.sh
rm convert.sh fortran2bd.py
cp -r ../../../org_sys/bd/ ../../BD/$i
mv restart_* ../../BD/$i
cp ../../../ini_struct/$i/ex_force.txt ../../BD/$i
cp ../../../ini_struct/data_mod ../../BD/$i
cp ../../../ini_struct/LH.in ../../BD/$i
cp ../../../../bin/code ../../BD/$i
cd ../
done

cd ../BD
for i in copy_*
do
cd $i
sbatch run-code.sbatch
cd ..
done


