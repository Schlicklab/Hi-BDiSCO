#!/bin/bash

cd output/ini_struct

for i in copy_*
do
cp -r ../org_sys/mc ../simulation/random_ini/$i
cp 1.dat ../simulation/random_ini/$i/
cp dim.in ../simulation/random_ini/$i/
cp fold_elig.dat ../simulation/random_ini/$i/
cp LHbound.0.in ../simulation/random_ini/$i/
cp ll_elig.dat ../simulation/random_ini/$i/
done

cd ../simulation/random_ini/
for i in copy_*
do
cd $i
cp ../../../../bin/chrom_ap1.x .
sbatch sub_batch.s
cd ../
done
