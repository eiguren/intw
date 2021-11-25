#!/usr/bin/env bash

# Execute some program
if [ "$#" -eq 2 ]; then
  num1=$1
  num2=$2
else
  echo "Wrong number of arguments"
  exit 1
fi
result=$((num1 + num2))


# Compare the results with previous reference calculations
reference_calculation_result=8
if (( $result == $reference_calculation_result )) ; then
  # Exit with error code 0 if results are correct
  exit 0
else
  # Exit with error code 1 if results are wrong
  exit 1
fi
