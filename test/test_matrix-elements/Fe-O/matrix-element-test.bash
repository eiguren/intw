#!/usr/bin/env bash

# Run the calculation

# cd intw
# ../../../../build/src/utilities/ep_melements_haritz_example_proba.x < ep_melements_k.in > ep_melements_k.out
# cd ..


# REFERENCE_MAT_FILE=/home/haritz/bulegoa/Kalkuluak/Probak/Fe-O/k.out_dvtot_all
REFERENCE_MAT_FILE=QE-4.3.2/fe-o-elph.ph.out_dvtot
TEST_MAT_FILE=intw/ep_melements_k.out

MIN_MAT_EL=0.0001
MAX_DIFF_PERCENTAGE=0.15

NUM_BANDS=15
NON_DEGENERATED_BANDS=( 1 2 5 6 9 14 15 )



cat $REFERENCE_MAT_FILE | grep "haritz " | awk '{printf("%2i %2i %2i %15.10f\n",$2,$3,$4,$7)}' > reference_file
cat $TEST_MAT_FILE | grep "haritz " | awk '{printf("%2i %2i %2i %15.10f\n",$2,$3,$4,$7)}' > test_file

echo "iband jband imode                 mat1                 mat2                 diff           percentage"
while IFS= read -r line
do
  #
  iband=$(echo $line | awk '{print $1}')
  iband_=$(echo $line | awk '{print $5}')
  if [[ ! " ${iband} " =~ " ${iband_} " ]]; then
    echo "Error"
    exit 1
  fi
  #
  jband=$(echo $line | awk '{print $2}')
  jband_=$(echo $line | awk '{print $6}')
  if [[ ! " ${jband} " =~ " ${jband_} " ]]; then
    echo "Error"
    exit 1
  fi
  #
  imode=$(echo $line | awk '{print $3}')
  imode_=$(echo $line | awk '{print $7}')
  if [[ ! " ${imode} " =~ " ${imode_} " ]]; then
    echo "Error"
    exit 1
  fi
  #
  if [[ " ${NON_DEGENERATED_BANDS[@]} " =~ " ${iband} " ]]; then
    if [[ " ${NON_DEGENERATED_BANDS[@]} " =~ " ${jband} " ]]; then
      MAT_EL=$(echo $line | awk '{printf("%20.15f", $4)}')
      if (( $(echo "$MAT_EL > $MIN_MAT_EL" | bc -l) )); then
        echo $line | awk '{printf( "%5i %5i %5i %20.15f %20.15f %20.15f %20.5f\n", $1,$2,$3,$4,$8,$4-$8,($4-$8)/$4*100)}'
        DIFF_PERCENTAGE=$(echo $line | awk '{printf( "%20.5f\n", sqrt((($4-$8)/$4*100)**2) )}')
        if (( $(echo "$DIFF_PERCENTAGE > $MAX_DIFF_PERCENTAGE" | bc -l) )); then
          echo "MAX_DIFF_PERCENTAGE exceeded!"
          exit 1
        fi
      fi
    fi
  fi
done <<< $(paste reference_file test_file)

rm reference_file test_file
