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
for ik in range(NKPTS):
    non_deg_k = []
    eig_file = "/home/haritz/Codes/intw/intw3.3/test/test_matrix-elements/Si-fullmesh_backup/QE-6.7.0/reference2/si.save.intw/wfc%05i.dat" % (ik+1)
    with FortranFile(eig_file, "r") as eig_file:
        eig_file.read_ints(np.int32)
        eig_file.read_ints(np.int32)
        eig = np.around(eig_file.read_reals(float), decimals=5)
        for ib in range(NBANDS):
            if any(abs(eig[np.arange(len(eig))!=ib]-eig[ib]) < 0.01):
                pass
            else:
                non_deg_k.append(ib+1)
    non_deg.append(non_deg_k)

print(non_deg_k[0])


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
