#!/usr/bin/env python
import numpy as np
import re
import sys

# Set accuracy for comparison
acc = 1.0E-5

# Read INTW generated file
out_file = "si_hr_intw.dat"
#
ndegen=[]
re_ham_r=[]
im_ham_r=[]
irvec=[]
#
i=0
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    i=i+1
    # Skip header
    if i==1:
       continue
    # Read general info
    elif i==2:
       num_wann = int(re.findall(r'[+-]?[0-9]+', line)[0])
       print('Number of Wannier functions:', num_wann)
    elif i==3:
       nrpts = int(re.findall(r'[+-]?[0-9]+', line)[0])
       print('Number of R-vectors within WS cell:', nrpts)

    else:
        # Read ndegen
        if (len(re.findall(r'[+-]?[0-9]+', line))>0 and len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==0):
            for ir in range(len(re.findall(r'[+-]?[0-9]+', line))):
                ndegen.append(int(re.findall(r'[+-]?[0-9]+', line)[ir]))

        # Read irvec and hr
        if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))>0):
            # Read band-pair
            ib1 = int(re.findall(r'[+-]?[0-9]+', line)[3])
            ib2 = int(re.findall(r'[+-]?[0-9]+', line)[4])
            # Read irvec for first band-pair only
            if ib1==1 and ib2==1:
                jb=0
                re_ham_r_jb=[]
                im_ham_r_jb=[]
                #
                irvec_ir=[]
                for ia in range(3):
                    irvec_ir.append(int(re.findall(r'[+-]?[0-9]+', line)[ia]))
                irvec.append(irvec_ir)
            # Read hr
            jb=jb+1
            re_ham_r_jb.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0]))
            im_ham_r_jb.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1]))
            if jb==num_wann**2:
                re_ham_r.append(re_ham_r_jb)
                im_ham_r.append(im_ham_r_jb)

#
ndegen=np.array(ndegen)
irvec=np.array(irvec)
re_ham_r=np.array(re_ham_r)
im_ham_r=np.array(im_ham_r)



# ----------------------------------------

# Read w90 generated file
out_file = "si_hr.dat"
#
ndegen_w90=[]
re_ham_r_w90=[]
im_ham_r_w90=[]
irvec_w90=[]
#
i=0
with open(out_file, "rt") as f:
    data = f.readlines()
for line in data:
    i=i+1
    # Skip header
    if i==1:
       continue
    # Read general info
    elif i==2:
       num_wann = int(re.findall(r'[+-]?[0-9]+', line)[0])
       print('Number of Wannier functions:', num_wann)
    elif i==3:
       nrpts = int(re.findall(r'[+-]?[0-9]+', line)[0])
       print('Number of R-vectors within WS cell:', nrpts)

    else:
        # Read ndegen
        if (len(re.findall(r'[+-]?[0-9]+', line))>0 and len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))==0):
            for ir in range(len(re.findall(r'[+-]?[0-9]+', line))):
                ndegen_w90.append(int(re.findall(r'[+-]?[0-9]+', line)[ir]))

        #
        if (len(re.findall(r'[+-]?[0-9]+\.[0-9]+', line))>0):
            # Read band-pair
            ib1 = int(re.findall(r'[+-]?[0-9]+', line)[3])
            ib2 = int(re.findall(r'[+-]?[0-9]+', line)[4])
            # Read irvec for first band-pair only
            if ib1==1 and ib2==1:
                jb=0
                re_ham_r_jb=[]
                im_ham_r_jb=[]
                #
                irvec_ir=[]
                for ia in range(3):
                    irvec_ir.append(int(re.findall(r'[+-]?[0-9]+', line)[ia]))
                irvec_w90.append(irvec_ir)
            # Read hr
            jb=jb+1
            re_ham_r_jb.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[0]))
            im_ham_r_jb.append(float(re.findall(r'[+-]?[0-9]+\.[0-9]+', line)[1]))
            if jb==num_wann**2:
                re_ham_r_w90.append(re_ham_r_jb)
                im_ham_r_w90.append(im_ham_r_jb)

#
ndegen_w90=np.array(ndegen_w90)
irvec_w90=np.array(irvec_w90)
#for ir in range(len(irvec_w90)):
#    print(irvec_w90[ir])
re_ham_r_w90=np.array(re_ham_r_w90)
im_ham_r_w90=np.array(im_ham_r_w90)

#for ir in range(len(irvec)):
#    print(irvec[ir])

# ----------------------------------------

# Compare
ierr=0
print(len(irvec))
for ir in range(nrpts):
    irvec_found=False    
    for jr in range(nrpts):
        # Detect same irvec
        irvec_diff = irvec[ir] - irvec_w90[jr]
        if (irvec_diff[0]==0 and irvec_diff[1]==0 and irvec_diff[2]==0):
            ## Print irvecs (for testing)
            #print(irvec[ir], irvec_w90[jr], irvec_diff)
            # Compare ndegen
            if (ndegen[ir] != ndegen_w90[jr]):
                print('ndegen not the same!')
                print(irvec[ir], irvec_w90[jr], irvec_diff)
                print(ndegen[ir], ndegen[jr])
                sys.exit(1)
                break
            # Compare H(R) for all band pairs
            nb=0
            for ib in range(num_wann):
                for jb in range(num_wann):
                    re_diff = np.abs(re_ham_r[ir][nb]-re_ham_r_w90[jr][nb])
                    im_diff = np.abs(im_ham_r[ir][nb]-im_ham_r_w90[jr][nb])
                    if (re_diff > acc or im_diff > acc):
                        print('hr not the same!')
                        print(irvec[ir], irvec_w90[jr], irvec_diff)
                        print(ib+1, jb+1)
                        print(re_ham_r[ir][nb], im_ham_r[ir][nb])
                        print(re_ham_r_w90[jr][nb], im_ham_r_w90[jr][nb])
                        ierr=1
                        sys.exit(1)
                        break
                    #else:
                        #print('hr OK!')
                        #print(irvec[ir], irvec_w90[jr], irvec_diff)
                        #print(ib+1, jb+1)
                        #print(re_ham_r[ir][nb], re_ham_r_w90[jr][nb])
                    #
                    nb=nb+1
            #
            irvec_found = True
            break

    # Check if irvec has not been found in w90 list
    if not (irvec_found):
        print('irvec not found in w90!')
        print(irvec[ir])                    
        sys.exit(1)
        break

#
print(ierr)
