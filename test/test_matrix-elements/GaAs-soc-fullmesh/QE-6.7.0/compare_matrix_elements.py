#!/usr/bin/env python3

from scipy.io import FortranFile
import numpy as np
import sys

# Find some variables
NAT = 2
NMODE = 3*NAT
NBANDS = 36
NKPTS = 4*4*4

MIN_MAT_EL = 0.01
MAX_DIFF_PERCENTAGE = 5.00


def QE2intw(k):
    return list(map(lambda x: x % 1.0, k))


def k2ik(k, klist):
    ik = []
    for i, kpt in enumerate(klist):
        if sum(np.abs(np.array(k) - np.array(kpt))) < 0.001:
            ik.append(i)
    if len(ik) > 1:
        print("ERROR(k2ik): more than one k-point found!")
        print(k)
        sys.exit(1)
    if len(ik) == 0:
        print("ERROR(k2ik): no k-point found!")
        print(k)
        sys.exit(1)
    return ik[0]


# Read k-point list and non-degenreated band indices
NON_DEG = []
KPT_LIST = []
with open("non_deg_file", "r") as f:
    for line in f.readlines():
        print(line)
        kpt_text, non_deg_text = line.split(":")
        kpt_QE = list(map(float, kpt_text.replace("k = (", "").replace(")", "").split(",")))
        kpt_intw = QE2intw(kpt_QE)
        KPT_LIST.append(kpt_intw)

        if non_deg_text == "\n":
            NON_DEG.append([None])
        else:
            NON_DEG.append(list(map(int, non_deg_text.strip().replace('\n', '').split(' '))))
print(NON_DEG)
for ik, kpt in enumerate(KPT_LIST):
    print(ik, kpt)



mat_r = np.fromfile("reference/ref_ep_mat.dat_1",  dtype=np.dtype('c16'))
mat_t = np.fromfile("ep_mat.dat_1",  dtype=np.dtype('c16'))


print("%5s %5s %5s %5s %20s %20s %20s %20s" %
      ("ikpt", "iband", "jband", "imode", "mat1", "mat2", "diff", "percentage"))
for im in range(NMODE):
    for ik in range(NKPTS):
        for ib in range(NBANDS):
            for jb in range(NBANDS):
                iel = NBANDS*NBANDS*NKPTS*im + NBANDS*NKPTS*ib + NKPTS*jb + ik
                me_r = np.abs(mat_r[iel])
                me_t = np.abs(mat_t[iel])
                if (ib+1 in NON_DEG[ik] and jb+1 in NON_DEG[ik]):
                    if (me_r > MIN_MAT_EL):
                        print("%5i %5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % (
                        ik+1, ib+1, jb+1, im+1, me_r,
                        me_t, me_r-me_t, (me_r-me_t)/me_r*100)
                        )
                        if (abs(me_r-me_t)/me_r*100 > MAX_DIFF_PERCENTAGE):
                            print("MAX_DIFF_PERCENTAGE exceeded!")
                            sys.exit(1)
