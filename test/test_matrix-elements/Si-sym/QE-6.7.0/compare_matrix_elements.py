#!/usr/bin/env python3

from scipy.io import FortranFile
import numpy as np
import sys

# Find some variables
NAT = 2
NMODE = 3*NAT
NBANDS = 8
NKPTS = 4*4*4

MIN_MAT_EL = 0.01
MAX_DIFF_PERCENTAGE = 5.00



non_deg = []
with open("non_deg_file", "r") as f:
    for line in f.readlines():
        if line=="\n":
            non_deg.append([None])
        else:
            non_deg.append(list(map(int, line.replace('\n', '').split(' '))))




mat_r = np.fromfile("reference/ref_ep_mat.dat_1",  dtype=np.dtype('c16'))
mat_t = np.fromfile("ep_mat.dat_1",  dtype=np.dtype('c16'))


print("%5s %5s %5s %5s %20s %20s %20s %20s" %
      ("ikpt", "iband", "jband", "imode", "mat1", "mat2", "diff", "percentage"))
for im in range(NMODE):
    for ik in range(NKPTS):
        for ib in range(NBANDS-1):
            for jb in range(NBANDS-1):
                iel = NBANDS*NBANDS*NKPTS*im + NBANDS*NKPTS*ib + NKPTS*jb + ik
                me_r = np.abs(mat_r[iel])
                me_t =np.abs(mat_t[iel])
                if (ib+1 in non_deg[ik] and jb+1 in non_deg[ik] and ib+1):
                    if (me_r > MIN_MAT_EL):
                        print("%5i %5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % (
                        ik+1, ib+1, jb+1, im+1, me_r,
                        me_t, me_r-me_t, (me_r-me_t)/me_r*100)
                        )
                        if (abs(me_r-me_t)/me_r*100 > MAX_DIFF_PERCENTAGE):
                            print("MAX_DIFF_PERCENTAGE exceeded!")
                            sys.exit(1)
