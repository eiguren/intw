# INTW project

Hau intw-ren lehenengo bertsioa da.




## Compilation notes

To build INTW CMake (>=3.12), [a supported build tool](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html) like Makefile or Ninja and C and Fortran compilers are required. All of them can be easily installed in most systems.

Additionally, a supported version of [Quantum Espresso](https://www.quantum-espresso.org/) or [SIESTA](https://siesta-project.org/siesta/) (or both) is needed. Currently, the following versions are supported:

| DFT code         | Supported versions                           |
|------------------|----------------------------------------------|
| Quantum Espresso | &bull; 6.6 <br> &bull; 6.7MaX <br>&bull; 6.8 |
| SIESTA           | &bull; 4.1.5                                 |


To build intw to work with both DFT codes:

```bash
mkdir build
cd build
QE_HOME=/path/to/QE SIESTA_HOME=/path/to/SIESTA cmake ..
make
```

If `QE_HOME` and `SIESTA_HOME` are defined as environmental variables CMake will use the value of the environmental variable.

If `QE_HOME` is not specified, CMake will disable the components related to QE automatically. Similarly, if `SIESTA_HOME` is not specified, the components related to SIESTA will be disabled.


### Setting a specific compiler and flags

CMake will find and select automatically the compilers and required libraries. However, it is possible to set a specific set of compilers using CMake's variable `CMAKE_<LANG>_COMPILER`:

```bash
cmake -DCMAKE_Fortran_COMPILER=ifort ..
```

Alternatively, the compiler can be specified using environmental variables:

```bash
export FC=ifort CC=icc
```

or

```bash
FC=ifort cmake ..
```


Similarly, CMake's variable `CMAKE_<LANG>_FLAGS` or `FFLAGS` and `CFLAGS` environmental variables can be used to set a specific set of compilation flags.


### Building intw in Atlas

To build intw in [Atlas](https://dipc.ehu.eus/en/supercomputing-center?set_language=en), one also needs to first load a toolchain and CMake.

In atlas-fdf use for example `module load intel/2019b` and `module load CMake/3.15.3-GCCcore-8.3.0`.

In atlas-edr, use for example `module load intel/2019b`, `module load CMake/3.20.1-GCCcore-10.3.0` and `module load Python/3.7.6-Anaconda3-2020.02`.

In addition, the compilers and the QE, SIESTA or W90 directories need to be specified as follows:

```bash
cmake -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_Compiler=mpiicc -DQE_HOME=/path/to/QE -DW_HOME=/path/to/W90 -DSIESTA_HOME=/path/to/SIESTA  ..
```

## Dependencies


### BLAS/LAPACK (required)

BLAS and LAPACK libraries are found by CMake `find_package` modules. To set a specific BLAS/LAPACK library use `BLA_VENDOR` variable. For example, to enforce CMake to use OpenBLAS use:

```bash
cmake -DBLA_VENDOR=OpenBLAS ..
```

The complete list of allowed BLAS/LAPACK vendors, together with other options to find BLAS and LAPACK with CMake can be found [here](https://cmake.org/cmake/help/latest/module/FindBLAS.html).


### Spglib (required to use SIESTA)

[Spglib](https://spglib.readthedocs.io/) is a C library for finding and handling crystal symmetries, which also has a Fortran interface.

When intw is used in conjunction with SIESTA, Spglib is used to find crystal symmetries and use them to reduce computational costs. In particular, Spglib is a used by the following components of intw:
- `siesta2ph`: Symmetries are used to reduce the required number of atomic displacements to completely determine the interatomic force constants and the induced potentials.
- `siesta2intw`: Symmetries are used to reduce the Brillouin zone sampling to the irreducible q-points.

CMake will attempt to find Spglib by searching for `libsymspg(.so|.a)` in the paths specified in the `LD_LIBRARY_PATH` environmental variable second. If the automatic search fails, `SPGLIB_LIBRARY` variable can be used to set the library:

```bash
cmake -DSPGLIB_LIBRARY=/path/to/libspglib(.so|a) ..
```

Alternatively, if it is supported by the installed Spglib version, `pkg-config` can be used to find the library:

```bash
cmake -DSPGLIB_PREFER_PKGCONFIG=ON ..
```

Remember adjusting the `PKG_CONFIG_PATH` environmental variable if you installed Spglib in a non-standard location.


### triangle and tetgen (optional)

[triangle](https://www.cs.cmu.edu/~quake/triangle.html) is a Two-Dimensional Quality Mesh Generator and Delaunay Triangulator. And [tetgen](https://wias-berlin.de/software/index.jsp?id=TetGen) is a program to generate tetrahedral meshes of any 3D polyhedral domains.

Both, `triangle` and `tetgen` codes, are used by INTW's `triFS.x` utility to obtain a fully symmetric triangulated Fermi surface. In order to `triFS.x` be compiled, ensure that both codes are installed and available in your system's `PATH`. Both codes can be downloaded directly from their official websites or, alternatively, in some distributions they can also be installed from the system's package repositories. E.g.:

```bash
# Ubuntu
sudo apt install triangle-bin tetgen
```

If `triangle` and `tetgen` executables are not found by CMake in the configuring step, `triFS.x` utility will not be compiled, however, the rest of the utilities will still be successfully built.


### OpenMP (optional experimental)

Some parts of intw are parallelized using OpenMP.

:heavy_exclamation_mark: The parallelization with OpenMP is currently untested.


### Wannier90 (optional)

Wannier90 is used only in the test suite. Read [Test suite](#test-suite) section for more information.


### Python 3 (optional)

Python 3 is used only in the test suite. Read [Test suite](#test-suite) section for more information.



## Test suite

Intw's test suite can be executed by using `ctest`, an executable that is available with CMake. Alternatively, when using the Makefile generator, `make test` can be used also.

The test suite has 5 differentiated parts:
- `test_dummy`: It does not contain any actual code test at all, but consists of some simplified tests to serve as examples for developers when creating new tests.
- `test_IO`: To check the IO routines of intw.
- `test_matrix_elements`: A set of tests on different systems to check the electron-phonon matrix element calculation with `ep_melements.x`.
- `test_mmn`: A set of tests on different systems to check the $A_{mn}$ and $M_{mn}$ matrices calculated with `intw2W90.x`.
- `test_w902intw`: A set of tests on different systems to check the interface `w902intw.x`.


### Test dependencies

Some components of the test suite need some specific software to be executed.


#### python3 (required)

Since most tests use `python3` scripts to compare the results, if CMake does not find any `python3` executable, the whole test suite will be disabled.


#### Wannier90 (required for `test_mmn` and `test_w902intw`)

[Wannier90](http://www.wannier.org/) is used in `test_mmn` and `test_w902intw` scripts. Only the `wannier90.x` executable is needed, but the Wannier90 directory has to be specified to CMake:

```bash
W_HOME=/path/to/WANNIER90 QE_HOME=/path/to/QE SIESTA_HOME=/path/to/SIESTA cmake ..
```

If `W_HOME` is not specified, `test_mmn` and `test_w902intw` will be disabled.


#### QE or SIESTA (required for `test_matrix_elements`, `test_mmn` and `test_w902intw`)

`test_matrix_elements`, `test_mmn` and `test_w902intw` tests run QE or SIESTA calculations, therefore, **at least one** of them must be specified to CMake, otherwise `test_matrix_elements`, `test_mmn` and `test_w902intw` will be disabled.
