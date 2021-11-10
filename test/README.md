# Test

This directory contains the test for testing intw.

Remember running ```make``` before executing the tests, since some tests will need to be compiled. Then run ```make test``` or ```ctest``` command to execute all the tests. ```ctest``` gives more options to execute the test, for example with ```ctest -V``` the output of the tests is showed.


## Some theory about testing
There are three main methods to test a software: Unit testing, Integration testing and Functional testing.

Unit testing means testing individual modules of an application in isolation (without any interaction with dependencies) to confirm that the code is doing things right.

Integration testing means checking if different modules are working fine when combined together as a group.

Functional testing means testing a slice of functionality in the system (may interact with dependencies) to confirm that the code is doing the right things.

## intw test suite

intw is prepared to develop tests in different ways:

 1. It is possible to write a script to execute a intw utiliy program, and then compare the obtained results with previous reference calculations (Functional testing).
 2. Another option is to write a fortran function, where it is possible to run individual intw subroutines (Unit testing), or even more complex programs (Integration testing).

#### 1. Scripts

The scripts could be written in different languages, for example bash or python. In the ```test_dummy``` folder there are examples with both languages. This scripts should run a program, check the results and compare them with a previous reference calculation, and exit with error code ```0``` if the result is the expected, or if the result is not correct with error code ```1```. This type of test would be well suited to check the code when major changes are made on it. One example of such a test could be for example to compute the electron-phonon matrix elements of a system and compare it with a previous calculation. In this way we can check if all the subroutines together work as they should be.

#### 2. Fortran functions
In this case, the function could be as simple as reading the calculation parameters, or as complex as computing the electron-phonon matrix elements. But I think that has no sense to write a complex function to compute the electron-phonon matrix elements, because the test should be compiled before runing the test, and for that purpouse I think that it would be better to write a script that runs an existing intw utility.

This type of tests should be used to check very specific parts of the code. For example, one test could check if the symmetry operations are readed and used correctly. Or another test could read the pseudopotentials and check if all the data is readed correctly.

On the ```test_IO``` file there is a very simple test of this type that could be used as an example.

## Test module
The ```module``` directory contains a recopilation of useful functions for testing.

## Adding new tests
Read ```CMakeList.txt``` to learn how to add new tests.
