#!/usr/bin/env python3

from scipy.io import FortranFile
import numpy as np
import sys
import re

PREFIX = "gaas-nscf"

NAT = 2
NMODE = 3*NAT
NBANDS = 18
NSPIN = 1
NK1 = 4
NK2 = 4
NK3 = 4
NKPTS = NK1*NK2*NK3
NQ1 = 1
NQ2 = 1
NQ3 = 1
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


def find_non_deg_list(eig):
    if len(eig) != NBANDS:
        print("ERROR(find_non_deg_list): wrong dimensions!")
        sys.exit(1)
    non_deg = []
    for ib in range(NBANDS-1):
        ei_deg = False
        for jb in range(NBANDS):
            if (ib != jb and abs(eig[ib] - eig[jb]) < 0.01):
                ei_deg = True
                break
        if not ei_deg:
            non_deg.append(ib+1)
    if len(non_deg) == 0:
        non_deg.append(None)
    return non_deg


def kmesh(n1, n2, n3):
    k_list = []
    for i1 in range(n1):
        for i2 in range(n2):
            for i3 in range(n3):
                kpt = QE2intw([i1/n1, i2/n2, i3/n3])
                k_list.append(kpt)
    return k_list


# Find some variables
with open(PREFIX+".save.intw/crystal.dat", "r") as f:
    lines = f.readlines()
    for i, line in enumerate(lines):

        text = line.strip()

        if text == "ALAT":
            alat = float(lines[i+1].strip())

        if text == "NAT":
            nat = int(lines[i+1].strip())
            nmode = 3*nat

        if text == "NBAND":
            nband = int(lines[i+1].strip())

        if text == "NKS":
            nkpts = int(lines[i+1].strip())

        if text == "LSPIN":
            lspin = str(lines[i+1].strip())
            if lspin == "T":
                nspin = 2
            else:
                nspin = 1

        if text == "AT":
            a1 = list(map(float, re.sub("\s+", ",", lines[i+1].strip()).split(",")))
            a2 = list(map(float, re.sub("\s+", ",", lines[i+2].strip()).split(",")))
            a3 = list(map(float, re.sub("\s+", ",", lines[i+3].strip()).split(",")))
            at = np.array([a1, a2, a3]).T

# print(alat)
# print(at)
# print(nat)
# print(nmode)
# print(nband)
# print(nkpts)
# print(nspin)
# sys.exit()

# Check the parameters
if nat != NAT:
    print("ERROR: wrong NAT")
    sys.exit(1)
if nband != NBANDS:
    print("ERROR: wrong NBANDS")
    sys.exit(1)
if nkpts != NKPTS:
    print("ERROR: wrong NKPTS")
    sys.exit(1)
if nspin != NSPIN:
    print("ERROR: wrong NSPIN")
    sys.exit(1)

# Read k-point list and find non-degenreated band indices
NON_DEG = []
KPT_LIST = []
with open(PREFIX+".save.intw/kpoints.dat", "r") as f:
    lines = f.readlines()
    for ik, line in enumerate(lines):
        kpt_QE_cart = list(map(float, re.sub("\s+", ",", line.strip()).split(",")))
        kpt_QE = np.matmul(at, kpt_QE_cart)
        kpt_intw = QE2intw(kpt_QE)
        KPT_LIST.append(kpt_intw)

        eig_file = PREFIX+".save.intw/wfc%05i.dat" % (ik+1)
        with FortranFile(eig_file, "r") as eig_file:
            eig_file.read_ints(np.int32)
            eig_file.read_ints(np.int32)
            eig = eig_file.read_reals(float)
        non_deg_k = find_non_deg_list(eig)
        NON_DEG.append(non_deg_k)

if len(KPT_LIST) != NKPTS:
    print("ERROR: wrong KPTS")
    sys.exit(1)

# for ik, kpt in enumerate(KPT_LIST):
#     print(ik, kpt, ":", NON_DEG[ik])

# Find q-mesh
QPT_MESH = kmesh(NQ1, NQ2, NQ3)

# Check q-points list
QPT_LIST = []
with open("qlist.txt", "r") as f:
    lines = f.readlines()
    for ik, line in enumerate(lines):
        iqpt_r, qx, qy, qz = map(string2number, re.sub("\s+", ",", line.strip()).split(','))
        qpt_QE_cart = np.array([qx, qy, qz])
        qpt_QE = np.matmul(at, qpt_QE_cart)
        qpt_intw = QE2intw(qpt_QE)
        if qpt_intw not in QPT_MESH:
            print("ERROR: q-point is not in the q-mesh")
            sys.exit(1)
        QPT_LIST.append(qpt_intw)

# if len(QPT_LIST) != NQPTS:
#     print("ERROR: wrong QPTS")
#     sys.exit(1)

#
# Compare the matrix elements
#

i1 = NKPTS
i2 = NKPTS*NBANDS
i3 = NKPTS*NBANDS*NBANDS
i4 = NKPTS*NBANDS*NBANDS*NSPIN
i5 = NKPTS*NBANDS*NBANDS*NSPIN*NSPIN
for iq in range(NQPTS):

    print('\nQ-point %i' % (iq+1))

    q = np.array(QPT_MESH[iq])
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
                    me_r = 0
                    me_t = 0
                    for is1 in range(NSPIN):
                        for is2 in range(NSPIN):
                            iel = i5*im + i4*is2 + i3*is1 + i2*jb + i1*ib + ik
                            me_r += mat_r[iel]
                            me_t += mat_t[iel]
                    me_r = np.abs(me_r)
                    me_t = np.abs(me_t)
                    if (ib+1 in NON_DEG[ikq] and jb+1 in NON_DEG[ik]):
                        if (me_r > MIN_MAT_EL):
                            print("%5i %5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % (
                                ik+1, ib+1, jb+1, im+1, me_r,
                                me_t, me_r-me_t, (me_r-me_t)/me_r*100)
                                )
                            if np.isnan(me_t) or (abs(me_r-me_t)/me_r*100 > MAX_DIFF_PERCENTAGE):
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
