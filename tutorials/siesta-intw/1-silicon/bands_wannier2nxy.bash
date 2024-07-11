#!/usr/bin/env bash

PREFIX=si


BANDS=$(cat ${PREFIX}_band.dat | sed -e's/[[:space:]]*$//' | awk 'BEGIN{b=0} {if ($0=="") {b=b+1}} END{print b}')
POINTS=$(cat ${PREFIX}_band.dat | sed -e's/[[:space:]]*$//' | awk 'BEGIN{k=0} {if ($0=="") {exit} else {k=k+1}} END{print k}')


for ((k = 1 ; k <= $POINTS ; k++)); do
    cat ${PREFIX}_band.dat | sed -e's/[[:space:]]*$//' | awk 'BEGIN{k=0} {if ($0=="") {k=0} else {k=k+1 ; printf("k%i %10.4f\n", k, $2)}}' | grep "k$k " | sed "s/k$k//g" | xargs -n $BANDS | paste <(echo "$k") -
done

