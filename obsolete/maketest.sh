#!/bin/bash
gfortran -o test.x src/kinds.f90 src/constants.f90 src/strings.f90 src/maths.f90 src/test.f90 
sleep 1 
rm *.mod