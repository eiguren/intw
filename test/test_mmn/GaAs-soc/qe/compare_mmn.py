#!/usr/bin/env python
import numpy as np
import re
import sys

# Read INTW generated file
out_file = "intw_fullmesh_gaas.mmn"
#
intw_remmn=[]
intw_immmn=[]
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==2):
        intw_remmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0]))
        intw_immmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1]))
#
intw_remmn=np.array(intw_remmn)
intw_immmn=np.array(intw_immmn)

# Read pw2wannier generated file
out_file = "qe_gaas.mmn"
#
qe_remmn=[]
qe_immmn=[]
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==2):
        qe_remmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0]))
        qe_immmn.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1]))
#
qe_remmn=np.array(qe_remmn)
qe_immmn=np.array(qe_immmn)


# Compare
ierr=0
for i in range(len(intw_remmn)):
    if (np.abs(intw_remmn[i] - qe_remmn[i]) > 1.0E-6) or (np.abs(intw_immmn[i] - qe_immmn[i]) > 1.0E-6) :
        print('Error in mmn element:', i)
        print('INTW value:', intw_remmn[i],'+ i', intw_immmn[i], ', modulus =', np.abs(intw_remmn[i]+1j*intw_immmn[i]))
        print('pw2wannier90 value:', qe_remmn[i],'+ i', qe_immmn[i], ', modulus =', np.abs(qe_remmn[i]+1j*qe_immmn[i]))
        ierr=1
        sys.exit(1)
        break
#
print(ierr)
