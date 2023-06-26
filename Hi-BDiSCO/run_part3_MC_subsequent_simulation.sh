#!/bin/bash

cd output/simulation/BD/

for i in copy_*
do
cd $i
cp ../../../../src/bd2fortran.py .
cp ../../../../src/run_bd2f.sh .
./run_bd2f.sh
rm run_bd2f.sh bd2fortran.py
cp -r ../../../org_sys/mc ../../MC/$i
cp 1.dat ../../MC/$i
cp ../../../../bin/chrom_ap1.x ../../MC/$i
cp ../../../ini_struct/1.dat ../../MC/$i/
cp ../../../ini_struct/dim.in ../../MC/$i/
cp ../../../ini_struct/fold_elig.dat ../../MC/$i/
cp ../../../ini_struct/LHbound.0.in ../../MC/$i/
cp ../../../ini_struct/ll_elig.dat ../../MC/$i/
cd ..
done

cd ../MC
for i in copy_*
do
cd $i
sbatch sub_batch.s
cd ..
done





