#!/usr/bin/env bash

QE_HOME=/home/jonl/CODES/qe-git
QE_VERSION=$(basename $PWD | sed 's/QE-//')

# Run all necessary QE calculations
make calculation
if [ $? -ne 0 ]; then
  exit 1
fi

# Compare WS and H(R), INTW vs. w90
python compare_hr.py
if [ $? -ne 0 ]; then
  exit 1
fi
