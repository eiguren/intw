#!/usr/bin/env python3

from scipy.io import FortranFile
import numpy as np
import sys
import re

# Find some variables
NAT = 1
NMODE = 3*NAT
NBANDS = 10
NK1 = 4
NK2 = 4
NK3 = 4
NKPTS = NK1*NK2*NK3
NQ1 = 2
NQ2 = 2
NQ3 = 2
NQPTS = NQ1*NQ2*NQ3

MIN_MAT_EL = 0.01
MAX_DIFF_PERCENTAGE = 5.00


def string2number(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


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


# Read unit cell vectors
with open("cu.save.intw/crystal.dat", "r") as f:

    text = re.sub("\s+", ",", f.readline().strip())
    if text != "ALAT":
        print("Error reading crystal.dat")
        sys.exit(1)
    alat = string2number(f.readline())

    text = re.sub("\s+", ",", f.readline().strip())
    if text != "AT":
        print("Error reading crystal.dat")
        sys.exit(1)
    a1 = list(map(float, re.sub("\s+", ",", f.readline().strip()).split(",")))
    a2 = list(map(float, re.sub("\s+", ",", f.readline().strip()).split(",")))
    a3 = list(map(float, re.sub("\s+", ",", f.readline().strip()).split(",")))
    at = np.array([a1, a2, a3]).T

# Read k-point list and non-degenreated band indices
NON_DEG = []
KPT_LIST = []
with open("non_deg_file", "r") as f:
    for line in f.readlines():
        kpt_text, non_deg_text = line.split(": ")
        kpt_QE = list(map(float, kpt_text.replace("k = (", "").replace(")", "").split(",")))
        kpt_intw = QE2intw(kpt_QE)
        KPT_LIST.append(kpt_intw)

        if non_deg_text == "\n":
            NON_DEG.append([None])
        else:
            NON_DEG.append(list(map(int, non_deg_text.replace('\n', '').split(' '))))
print(NON_DEG)
for ik, kpt in enumerate(KPT_LIST):
    print(ik, kpt)

# Read q-point list
QPT_LIST = []
with open("qlist.txt", "r") as f:
    for iq in range(NQPTS):
        line = re.sub("\s+", ",", f.readline().strip())
        iqpt_r, qx, qy, qz = map(string2number, line.split(','))
        qpt_QE = np.matmul(at, np.array([qx, qy, qz]))
        qpt_intw = QE2intw(qpt_QE)
        QPT_LIST.append(qpt_intw)

# Compare the matrix elements

for iq in range(NQPTS):

    print('\nQ-point %i' % (iq+1))

    q = np.array(QPT_LIST[iq])
    print(q)

    mat_el_file_r = 'reference/ref_ep_mat.dat_%i' % (iq+1)
    mat_r = np.fromfile(mat_el_file_r,  dtype=np.dtype('c16'))

    mat_el_file_t = 'ep_mat.dat_%i' % (iq+1)
    mat_t = np.fromfile(mat_el_file_t,  dtype=np.dtype('c16'))

    print("%5s %5s %5s %5s %20s %20s %20s %20s" %
          ("ikpt", "iband", "jband", "imode", "mat1", "mat2", "diff", "percentage"))

    for im in range(NMODE):
        for ik in range(NKPTS):
            k = np.array(KPT_LIST[ik])
            # print(k)
            kq = QE2intw(k + q)
            # print(kq)
            ikq = k2ik(kq, KPT_LIST)
            # print(ikq)
            for ib in range(NBANDS):
                for jb in range(NBANDS):
                    iel = NBANDS*NBANDS*NKPTS*im + NBANDS*NKPTS*jb + NKPTS*ib + ik
                    me_r = np.abs(mat_r[iel])
                    me_t = np.abs(mat_t[iel])
                    if (ib+1 in NON_DEG[ikq] and jb+1 in NON_DEG[ik]):
                        if (me_r > MIN_MAT_EL):
                            print("%5i %5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % (
                                  ik+1, ib+1, jb+1, im+1, me_r,
                                  me_t, me_r-me_t, (me_r-me_t)/me_r*100)
                                  )
                            if (abs(me_r-me_t)/me_r*100 > MAX_DIFF_PERCENTAGE):
                                pass
                                print("MAX_DIFF_PERCENTAGE exceeded!")
                                print("ik =", ik+1)
                                print("KPT_LIST[ik] =", KPT_LIST[ik])
                                print("k =", k)
                                print("NON_DEG[ik] =", NON_DEG[ik])
                                print("ikq =", ikq+1)
                                print("KPT_LIST[ikq] =", KPT_LIST[ikq])
                                print("kq =", kq)
                                print("NON_DEG[ikq] =", NON_DEG[ikq])
                                sys.exit(1)
