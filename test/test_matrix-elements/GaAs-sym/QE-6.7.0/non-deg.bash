#!/usr/bin/env bash


# Find non degenerate energies on each k point and save to a file

OUTPUT=gaas.nscf.out

NAT=$(cat $OUTPUT | grep "number of atoms/cell" | tail -n 1 | awk '{print $5}')
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
  echo $i
  l=$(echo $i+1 | bc)
  kpt=$(cat $OUTPUT | grep -A $l "number of k points=" | tail -n 1 | tr ',' ' ' | tr ')' ' ' | awk '{printf("%7.4f%7.4f%7.4f\n",$5,$6,$7)}')
  echo $kpt
  ENERGIES=$(cat $OUTPUT | grep -A $LIN "\\$kpt" | tail -n $LIN | tr '\n' ' ')
  NON_DEG_ENERGIES=$(echo $ENERGIES | tr ' ' '\n' | uniq -u | tr '\n' ' ')
  #
  ib=0
  NON_DEG=""
	for e in $ENERGIES;
	do
    #echo $e
    e=$(echo $e | sed -e 's/-/\\\-/')
		ib=$((ib+1))
    if [ $ib -eq $NBANDS ];
    then
      break
    fi
		echo $NON_DEG_ENERGIES | grep -w -q "$e" && NON_DEG=$(echo "$NON_DEG $ib")
	done
  #echo $NON_DEG
	echo $NON_DEG >> non_deg_file
done


