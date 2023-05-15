# INTW compilation notes

## 05/2023

At this time, the current INTW version is 3.3 hosted in github [here](https://github.com/eiguren/intw.git)

It is compiled with CMake. It requires QE to be compiled with CMake as well. 
Currently supports QE v6.7. The official release contains a bug that does not allow the compilation with CMake, so this particular commit at the [QE github repository](https://github.com/QEF/q-e.git) must be checked out:

626d07a1cff06525d08d931a12672b28717f9932

In practice, to compile this version of QE (in atlas-edr, one also needs to first load `module load intel/2019b` and `module load CMake/3.20.1-GCCcore-10.3.0`):

``
mkdir qe_git
cd qe_git
git init
git remote add origin https://github.com/QEF/q-e.git
git pull origin master
git fetch
git checkout 626d07a1cff06525d08d931a12672b28717f9932
mkdir build
cd build
FC=mpiifort CC=mpiicc cmake ..
make -j 8
``

We also need to compile wannier90 in library mode. The current version of INTW is compatible with the latest wannier90 v3.1.0, which can be downloaded [here](http://www.wannier.org/download/).
To compile, just copy a make.inc, make sure that libraries are properly linked (compare to QE compilation make.inc, for example), and type:

``
make
make lib
``

After this, we are ready to compile INTW. To do that, go to the INTW directory, and type:

``
mkdir build
cd build
FC=mpiifort CC=mpiicc QE_HOME=/path/to/QE W_HOME=/path/to/W cmake ..
make -j 8
``
