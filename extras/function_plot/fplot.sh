#!/bin/bash
export PYTHONPATH=$PYTHONPATH:"/cloud/Code/python/eampapy/f2py"

# Set the number of threads to 1
export OMP_NUM_THREADS=5

# Run
python3 fplot.py "$@"

