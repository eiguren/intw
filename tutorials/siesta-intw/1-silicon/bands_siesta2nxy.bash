#!/usr/bin/env bash

PREFIX=si

BANDS=$(head -n 3 $PREFIX.bands | tail -n 1 | awk '{print $1}')

POINTS=$(head -n 3 $PREFIX.bands | tail -n 1 | awk '{print $3}')

LINES=$(echo "$BANDS/8" | bc)
L=$(echo "$BANDS/8" | bc -l)
if (( $(echo "$L > $LINES" |bc -l) )); then
  LINES=$(echo "$LINES+1" | bc)
fi
LINES=$(echo "$LINES*$POINTS" | bc)


tail -n $LINES $PREFIX.bands | sed -e ':a;N;$!ba;s/\n    //g' | tr -s ' ' | awk '{$1=$2=$3=""; print $0}' | cat -n
