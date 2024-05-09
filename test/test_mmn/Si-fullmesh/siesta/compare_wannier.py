#!/usr/bin/env python
import numpy as np
import re
import sys

# Set accuracy threshold for comparison
acc = 1.0E-3

# Read INTW+W90 generated file from IBZ calculation
out_file = "intw_si.wout"
#
w_center=[]
with open(out_file, "rt") as f:
    data = f.readlines()
j=0
for line in data:
    j+=1
    if "Final State" in line:
        for line in data[j:len(data)-1]:
            # Save sum of wannier centers:
            if "Sum of centres and spreads" in line:
                x = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])
                y = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1])
                z = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[2])
            # Save sum of spreads:
            if "Omega Total" in line:
                spread = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])



# Read QE+W90 generated file
out_file = "siesta_si.wout"
#
with open(out_file, "rt") as f:
    data = f.readlines()
j=0
for line in data:
    j+=1
    if "Final State" in line:
        for line in data[j:len(data)-1]:
            # Save sum of wannier centers:
            if "Sum of centres and spreads" in line:
                ref_x = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])
                ref_y = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1])
                ref_z = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[2])
            # Save sum of spreads:
            if "Omega Total" in line:
                ref_spread = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])


# Compare
ierr=0
# Compare Wannier centers and spreads:
if (np.sqrt( (x - ref_x)**2 + (y - ref_y)**2 + (z - ref_z )**2 ) > acc
    or np.abs( spread - ref_spread) > acc):
    print('Error in Wannier functions')
    print('Sum of Wannier function centers:', x, y, z)
    print('Sum of Wannier function spreads', spread)
    print('(reference) Sum of Wannier function centers:', ref_x, ref_y, ref_z)
    print('(reference) Sum of Wannier function spreads', ref_spread)
    ierr=1
    sys.exit(1)
#
print(ierr)
