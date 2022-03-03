#!/usr/bin/env python3

import sys

# Find some variables
NAT = 1
NMODE = 3*NAT
NBANDS = 10
NKPTS = 4*4*4

MIN_MAT_EL = 0.01
MAX_DIFF_PERCENTAGE = 5.00


def string2number(s):
    try:
        return int(s)
    except:
        return float(s)


NON_DEG = []
with open("non_deg_file", "r") as f:
    for line in f.readlines():
        NON_DEG.append(list(map(int, line.replace('\n', '').split(' '))))

print("%5s %5s %5s %5s %20s %20s %20s %20s" %
      ("ikpt", "iband", "jband", "imode", "mat1", "mat2", "diff", "percentage"))
with open("reference_file", "r") as reference:
    with open("test_file", "r") as test:
        for ik in range(NKPTS):
            for ib in range(NBANDS):
                for jb in range(NBANDS):
                    for im in range(NMODE):
                        ikpt_r, iband_r, jband_r, imode_r, me_r = map(
                            string2number, reference.readline().split(' '))
                        ikpt_t, iband_t, jband_t, imode_t, me_t = map(
                            string2number, test.readline().split(' '))
                        # Some checks
                        if (ikpt_r == 223):
                            pass
                        else:
                            if (ikpt_r != ikpt_t):
                                print("Error on ikpt index")
                                sys.exit(1)
                            if (iband_r != iband_t):
                                print("Error on iband index")
                                sys.exit(1)
                            if (jband_r != jband_t):
                                print("Error on jband index")
                                sys.exit(1)
                            if (imode_r != imode_t):
                                print("Error on imode index")
                                sys.exit(1)
                            if (iband_r in NON_DEG[ikpt_r-1] and jband_r in NON_DEG[ikpt_r-1]):
                                if (me_r > MIN_MAT_EL):
                                    print("%5i %5i %5i %5i %20.15f %20.15f %20.15f %20.5f" % (
                                    ikpt_r, iband_r, jband_r, imode_r, me_r,
                                    me_t, me_r-me_t, (me_r-me_t)/me_r*100)
                                    )
                                    if (abs(me_r-me_t)/me_r*100 > MAX_DIFF_PERCENTAGE):
                                        print("MAX_DIFF_PERCENTAGE exceeded!")
                                        #sys.exit(1)
