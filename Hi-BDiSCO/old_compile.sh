#!/bin/bash
#module purge
#module load openmpi/intel/4.1.1
#module load cuda/11.3.1
#module load gcc/6.3.0

mpicxx -std=c++11 -c -Wall main.cpp readfile.h readfile.cpp constants.h func.h func.cpp mt.h mt.cpp utilities.h utilities.cpp
nvcc -c func_cuda.cu -lcusolver -lcublas
mpicxx *.o -lcudart -lcusolver -lcublas -o code

rm *.o
