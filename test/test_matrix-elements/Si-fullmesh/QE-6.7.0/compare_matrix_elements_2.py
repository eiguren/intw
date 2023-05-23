#!/usr/bin/env python3

import numpy as np
import sys

# Find some variables
NAT = 2
NMODE = 3*NAT
NBANDS = 8
NKPTS = 4*4*4

MIN_MAT_EL = 0.00001
MAX_DIFF_PERCENTAGE = 5.00


mat_r = np.fromfile("reference/ref_ep_mat.dat_1",  dtype=np.dtype('c16'))
mat_t = np.fromfile("ep_mat.dat_1",  dtype=np.dtype('c16'))


print("%5s %5s %20s %20s %20s %20s" %
      ("ikpt", "imode", "Tr[mat1]", "Tr[mat2]", "diff", "percentage"))
for im in range(NMODE):
    for ik in range(NKPTS):
        ep_mat_r = np.zeros((NBANDS,NBANDS))
        ep_mat_t = np.zeros((NBANDS,NBANDS))
        for ib in range(NBANDS):
            for jb in range(NBANDS):
                iel = NBANDS*NBANDS*NKPTS*im + NBANDS*NKPTS*ib + NKPTS*jb + ik
                ep_mat_r[ib,jb] = np.abs(mat_r[iel])
                ep_mat_t[ib,jb] = np.abs(mat_t[iel])
                #
        Tr_r = np.trace(ep_mat_r)
        Tr_t = np.trace(ep_mat_t)
        Tr_diff = Tr_r-Tr_t
        Tr_perc = Tr_diff/Tr_r*100
        if (Tr_r > MIN_MAT_EL):
            print("%5i %5i %20.15f %20.15f %20.15f %20.5f" % (ik+1, im+1, Tr_r, Tr_t, Tr_diff, Tr_perc) )
            if Tr_perc > MAX_DIFF_PERCENTAGE:
                print("MAX_Tr_DIFF exceded!")
                sys.exit(1)
