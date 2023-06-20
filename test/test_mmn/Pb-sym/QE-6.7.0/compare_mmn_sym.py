#!/usr/bin/env python
import numpy as np
import re
import sys

# Set accuracy threshold for comparison
acc = 1.0E-4

# Read INTW generated file from IBZ calculation
out_file = "pb.mmn"
#
intw_remmn=[]
intw_immmn=[]
intw_modmmn=[]
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==2):
        r = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])
        i = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1])
        intw_remmn.append(r)
        intw_immmn.append(i)
        intw_modmmn.append(np.sqrt(r**2+i**2))
#
intw_remmn=np.array(intw_remmn)
intw_immmn=np.array(intw_immmn)
intw_modmmn=np.array(intw_modmmn)

# Read previously generated benchmark file
# Tested that it gives same Wannier functions as full mesh calculation
out_file = "reference/pb.mmn"
#
ref_remmn=[]
ref_immmn=[]
ref_modmmn=[]
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==2):
        r = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0])
        i = float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1])
        ref_remmn.append(r)
        ref_immmn.append(i)
        ref_modmmn.append(np.sqrt(r**2+i**2))
#
ref_remmn=np.array(ref_remmn)
ref_immmn=np.array(ref_immmn)
ref_modmmn=np.array(ref_modmmn)


# Compare
ierr=0
for i in range(len(intw_remmn)):
    ## Compare real and imaginary parts separately
    #if (np.abs(intw_remmn[i] - ref_remmn[i]) > 1.0E-6) or (np.abs(intw_immmn[i] - ref_immmn[i]) > 1.0E-6) :
    # Compare magnitude
    if (np.abs(intw_modmmn[i] - ref_modmmn[i]) > acc):
        print('Error in mmn element:', i)
        print('INTW value:', intw_remmn[i],'+ i', intw_immmn[i], np.abs(intw_remmn[i]+1j*intw_immmn[i]))
        print('pw2wannier90 value:', ref_remmn[i],'+ i', ref_immmn[i], np.abs(ref_remmn[i]+1j*ref_immmn[i]))
        print('Difference in magnitude:', np.abs(intw_modmmn[i] - ref_modmmn[i]) )
        ierr=1
        sys.exit(1)
        break
#
print(ierr)
