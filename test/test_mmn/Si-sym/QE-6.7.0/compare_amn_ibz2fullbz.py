#!/usr/bin/env python
import numpy as np
import re
import sys

#
prefix = 'si'

# Set accuracy for comparison
acc = 1.0E-3

# Detect degenerate bands
out_file = prefix+".eig"
band, kpoint, eigen = np.loadtxt(out_file, unpack=True)
nb = np.int(np.max(band))
nk = np.int(np.max(kpoint))
#
degen = np.zeros((nk, nb), dtype=bool)
j=0
for ik in range(nk):
    ep=1.0E8
    for ib in range(nb):
        e=eigen[j]
        if (np.abs(e-ep)<1.0E-6):
            degen[ik, ib] = True
            degen[ik, ib-1] = True
        #
        ep=e
        j=j+1
 

# Read INTW generated file starting from IBZ wave functions
out_file = prefix+".amn"
#
band_1, proj_1, kpoint_1, re_amn_1, im_amn_1 = np.loadtxt(out_file, unpack=True, skiprows=2)
#
band_1 = band_1.astype(int)
proj_1 = proj_1.astype(int)
kpoint_1 = kpoint_1.astype(int)



# Read INTW generated file using full BZ wave functions
out_file = "../../Si-fullmesh/QE-6.7.0/"+prefix+".amn"
band_2, proj_2, kpoint_2, re_amn_2, im_amn_2 = np.loadtxt(out_file, unpack=True, skiprows=2)
#
band_2 = band_2.astype(int)
proj_2 = proj_2.astype(int)
kpoint_2 = kpoint_2.astype(int)


# Compare
ierr=0
for i in range(len(re_amn_1)):
    # Only compare non-degenerate nk-s
    if not degen[kpoint_1[i]-1, band_1[i]-1]: # Lists start at 0
        mod_1 = np.abs(re_amn_1[i]+1j*im_amn_1[i])
        mod_2 = np.abs(re_amn_2[i]+1j*im_amn_2[i])
        # Check magnitude
        if (np.abs(mod_1 - mod_2) > acc):
            print('Error in mmn element:', i)
            print('k-point:', kpoint_1[i])
            print('band:', band_1[i])
            print('projection:', proj_1[i])
            print('INTW value:', re_amn_1[i],'+ i', im_amn_1[i], "; Modulus:", mod_1)
            print('pw2wannier90 value:', re_amn_2[i],'+ i', im_amn_2[i], "; Modulus:", mod_2)
            ierr=1
            sys.exit(1)
            break
#
print(ierr)
