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
  l_cart=$(echo $i+1 | bc)
  l_crys=$(echo $i+$NKPTS+3 | bc)
  kpt_cart=$(cat $OUTPUT | grep -A $l_cart "number of k points=" | tail -n 1 | tr ',' ' ' | tr ')' ' ' | awk '{printf("%7.4f%7.4f%7.4f\n",$5,$6,$7)}')
  kpt_crys=$(cat $OUTPUT | grep -A $l_crys "number of k points=" | tail -n 1 | tr ',' ' ' | tr ')' ' ' | awk '{printf("%10.4f,%10.4f,%10.4f\n",$5,$6,$7)}')
  echo $kpt_cart
  #
  ENERGIES=$(cat $OUTPUT | grep -A $LIN "\\$kpt_cart" | tail -n $LIN | tr '\n' ' ' | sed -e's/  */ /g' | sed 's/^ *//g' | tr ' ' '\n' | awk '{printf("%10.2f",$1)}' | sed -e's/  */ /g' | sed 's/^ *//g')
  NON_DEG_ENERGIES=$(echo $ENERGIES | tr ' ' '\n' | uniq -u | tr '\n' ' ')
  #
  ib=0
  NON_DEG="k = ($kpt_crys):"
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


