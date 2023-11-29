#!/usr/bin/env python
import numpy as np
import re
import sys

# Set accuracy for comparison
acc = 1.0E-4

# Detect degenerate bands
out_file = "si.eig"
band, kpoint, eigen = np.loadtxt(out_file, unpack=True)
nb = int(np.max(band))
nk = int(np.max(kpoint))
#
degen = np.zeros((nk, nb), dtype=bool)
j = 0
for ik in range(nk):
    ep=1.0E8
    for ib in range(nb):
        e=eigen[j]
        if (np.abs(e-ep)<1.0E-3):
            degen[ik, ib] = True
            degen[ik, ib-1] = True
        if ib == nb-1:
            degen[ik, ib] = True
        #
        ep=e
        j=j+1

# Read INTW generated file
out_file = "intw_si.mmn"
#
with open(out_file, "rt") as f:
    data = f.readlines()
#
intw_nbands, intw_nkpts, intw_nprojs = map(int, data[1].split())
#
intw = []
for i in range(len(data[2:])):
    line = data[i+2]
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)) != 2):
        kpt = list(map(int, line.split()[0:2]))
        remmn = []
        immmn = []
        for j in range(intw_nbands**2):
            remmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', data[i+j+3])[0]))
            immmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', data[i+j+3])[1]))
        intw.append([kpt, np.array([remmn, immmn]).T])


# Read siesta generated file
out_file = "siesta_si.mmn"
#
with open(out_file, "rt") as f:
    data = f.readlines()
#
siesta_nbands, siesta_nkpts, siesta_nprojs = map(int, data[1].split())
#
siesta = []
for i in range(len(data[2:])):
    line = data[i+2]
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)) != 2):
        kpt = list(map(int, line.split()[0:2]))
        remmn = []
        immmn = []
        for j in range(siesta_nbands**2):
            remmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', data[i+j+3])[0]))
            immmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', data[i+j+3])[1]))
        siesta.append([kpt, np.array([remmn, immmn]).T])


# Compare
if intw_nbands != siesta_nbands:
    print("ERROR: wrong nbands")
    sys.exit(1)
if intw_nkpts != siesta_nkpts:
    print("ERROR: wrong nkpts")
    sys.exit(1)
if intw_nprojs != siesta_nprojs:
    print("ERROR: wrong nprojs")
    sys.exit(1)

ierr = 0
for i in range(len(siesta)):
    siesta_kpt = siesta[i][0]
    intw_kpt = intw[i][0]
    if intw_kpt != siesta_kpt:
        print("ERROR: wrong kpt")
        sys.exit(1)
    for j, (siesta_mmn, intw_mmn) in enumerate(zip(siesta[i][1], intw[i][1])):
        iband = j%siesta_nbands
        jband = int(j%(siesta_nbands**2)/siesta_nbands)
        # Only compare non-degenerate nk-s
        if (not degen[intw_kpt[0]-1, iband]) and (not degen[intw_kpt[1]-1, jband]):
            siesta_modmmn = np.abs(siesta_mmn[0] + 1j*siesta_mmn[1])
            intw_modmmn = np.abs(intw_mmn[0] + 1j*intw_mmn[1])
            if (np.abs(intw_modmmn - siesta_modmmn) > acc):
                print(intw_kpt)
                print(degen[intw_kpt[0]-1])
                print(degen[intw_kpt[1]-1])
                print('SIESTA value:', siesta_mmn[0],'+ i', siesta_mmn[1], siesta_modmmn)
                print('INTW value:', intw_mmn[0],'+ i', intw_mmn[1], intw_modmmn)
                print('Error in mmn element:', i, j)
                ierr = 1
                sys.exit(1)
                break

print(ierr)
