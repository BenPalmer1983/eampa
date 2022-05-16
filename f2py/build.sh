#!/bin/bash
python3 -m numpy.f2py --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
kinds/kinds.f90 \
math/math.f90 \
sls/sls.f90 \
transforms/transforms.f90 \
polyfit/polyfit.f90 \
interp/interp.f90 \
splines/splines.f90 \
optimizer/optimizer.f90 \
fnc/fnc.f90 \
atom/atom.f90 \
-m f_toolbox
cp *.so ../src

# python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math"