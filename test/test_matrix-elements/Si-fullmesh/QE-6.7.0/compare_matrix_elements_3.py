#!/usr/bin/env python3

import numpy as np
import sys

# Find some variables
NAT = 2
NMODE = 3*NAT
NBANDS = 8
NKPTS = 4*4*4

MIN_MAT_EL = 0.00001
MAX_DIFF_PERCENTAGE = 10.00


mat_r = np.fromfile("reference/ref_ep_mat.dat_1",  dtype=np.dtype('c16'))
mat_t = np.fromfile("ep_mat.dat_1",  dtype=np.dtype('c16'))


print("%5s %5s %5s %20s %20s %20s %20s" %
      ("ikpt", "iband", "imode", "Eig[mat1]", "Eig[mat2]", "diff", "percentage"))
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
        Eig_r = np.sort(np.linalg.eigvals(ep_mat_r))
        Eig_t = np.sort(np.linalg.eigvals(ep_mat_t))
        for ib in range(NBANDS):
            if Eig_r[ib] > MIN_MAT_EL:
                print("%5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % ( ik+1, ib+1, im+1, Eig_r[ib], Eig_t[ib], Eig_r[ib]-Eig_t[ib], (Eig_r[ib]-Eig_t[ib])/Eig_r[ib]*100 ))
                if (Eig_r[ib]-Eig_t[ib])/Eig_r[ib]*100 > MAX_DIFF_PERCENTAGE:
                    print("MAX_Eig_DIFF exceded!")
                    sys.exit(1)

