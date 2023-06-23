#!/bin/bash

a=`head -2 dim.in |tail -1`
b=`expr $a \* 82`
c=`tail -$a dim.in | awk '{Total=Total+$1} END{print 4*Total}'`
line_per_frame=`expr $b + $c`

tail -$line_per_frame fort.10 > out.fort

python3 fortran2bd.py
