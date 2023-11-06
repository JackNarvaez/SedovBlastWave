#!/usr/bin/env bash

files=('sedov.f90' 'main.f90')
fflags='-O3'
out='sedov'

clear

gfortran ${fflags} -o ${out} ${files[@]}
./${out}

gnuplot -p plot.gp

rm -f sedov *.mod  *.gp