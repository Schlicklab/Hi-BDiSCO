#!/bin/bash

a=`head -1 out.xyz`
tail -$a out.xyz > out_f.xyz

python3 bd2fortran.py
