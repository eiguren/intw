#!/usr/bin/env python
import sys

if __name__ == "__main__":

    # Execute some program
    num1 = 3
    num2 = 5
    result = num1 + num2

    # Compare the results with previous reference calculations
    reference_calculation_result = 8

    if result == reference_calculation_result:
        # Exit with error code 0 if results are correct
        sys.exit(0)
    else:
        # Exit with error code 1 if results are wrong
        sys.exit(1)
