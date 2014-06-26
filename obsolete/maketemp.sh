#!/bin/bash
gfortran -o temp.x src/kinds.f90 src/constants.f90 src/strings.f90 src/maths.f90 src/temp.f90 
sleep 1 
rm *.mod