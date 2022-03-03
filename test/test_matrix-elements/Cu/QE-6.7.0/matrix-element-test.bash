#!/usr/bin/env bash


# Find non degenerate energies on each k point and save to a file

OUTPUT=qe/cu.scf.out

NAT=$(cat $OUTPUT | grep "number of atoms/cell" | tail -n 1 | awk '{print $5}')
NMODE=$(echo 3*$NAT | bc)
NBANDS=$(cat $OUTPUT | grep "number of Kohn-Sham states" | awk '{print $5}')
NKPTS=$(cat $OUTPUT | grep "number of k points=" | head -n 1 | awk '{print $5}')


rm -rf non_deg_file

for i in $(seq $NBANDS);
do
  l=$(echo 8*$i | bc)
  if [ $l -ge $NBANDS ];
  then
    LIN=$(echo $i+1 | bc)
    break
  fi
done


for i in $(seq $NKPTS);
do
  l=$(echo $i+1 | bc)
  kpt=$(cat $OUTPUT | grep -A $l "number of k points=" | tail -n 1 | tr ',' ' ' | tr ')' ' ' | awk '{printf("%7.4f%7.4f%7.4f\n",$5,$6,$7)}')
  ENERGIES=$(cat $OUTPUT | grep -A $LIN "\\$kpt" | tail -n $LIN | tr '\n' ' ')
  NON_DEG_ENERGIES=$(echo $ENERGIES | tr ' ' '\n' | uniq -u | tr '\n' ' ')
  #
  ib=0
  NON_DEG=""
	for e in $ENERGIES;
	do
		ib=$((ib+1))
		echo $NON_DEG_ENERGIES | grep -w -q "$e" && NON_DEG=$(echo "$NON_DEG $ib")
	done
	echo $NON_DEG >> non_deg_file
done


# create refenrece and test files

REFERENCE_MAT_FILE=qe/cu.ph.out_haritz
TEST_MAT_FILE=ep_mat.dat_1

rm -rf reference_file
touch reference_file
cat $REFERENCE_MAT_FILE | grep haritz | awk '{print $2,$3,$4,$5,$8}' > reference_file

rm -rf test_file
touch test_file
for NU in $(seq $NMODES);
do
  for K in $(seq $NKPTS);
  do
    for M in $(seq $NBANDS);
    do
      for N in $(seq $NBANDS);
      do
	J=$(echo "16*$NBANDS*$NBANDS*$NKPTS*($NU-1) + 16*$NBANDS*$NKPTS*($M-1) + 16*$NKPTS*($N-1) + 16*($K-1)" | bc -l)
      	echo $K $N $M $NU $(od -A d -t f8 -j $J -N 16 $TEST_MAT_FILE | head -n 1 | tr "," "."| awk '{printf("%15.10f%15.10f    %15.10f\n",$2,$3,sqrt($2**2+$3**2))}') | awk '{print $1,$2,$3,$4,$7}'>> test_file
      done
    done
  done
done



