

# Using intw with SIESTA


### Table of contents
- [Silicon band structure with SIESTA and Wannier](#silicon-band-structure-with-siesta-and-wannier)
   - [Band structure Calculation with SIESTA](#band-structure-calculation-with-siesta)
   - [intw's Wannier90 interface](#intws-wannier90-interface)
   - [SIESTA's Wannier90 interface](#siestas-wannier90-interface)
- [Graphene phonon dispersion](#graphene-phonon-dispersion)
- [Copper electron-phonon matrix elements](#copper-electron-phonon-matrix-elements)


In this tutorial we will use the 4.1.5 version of SIESTA.

SIESTA uses PSF pseudo-potentials. (In the recently released 5.0.0 version PSML pseudo-potential from pseudo-dojo can be used).

The input for the calculation is given in Flexible Data Format (FDF).


## Silicon band structure with SIESTA and Wannier

In this example, first of all, we will do a small introduction to a band structure calculation with SIESTA. And in second place, we will see how to use the interface `siesta2intw.x` and intw's `intw2W90.x` utility to construct the Wannier functions, comparing the Wannier interpolated band structure with SIESTA's band structure. At the end, we will also show how can be used SIESTA's own interface to `wannier90.x`.

All files to run this example can be found in the `1-silicon` directory.

```
$ ls 1-silicon/
intw.in                 s2intw.in    si.win
bands_siesta2nxy.bash   silicon.fdf
bands_wannier2nxy.bash  Si.psf
```

[Back to top :arrow_heading_up:](#using-intw-with-siesta)

### Band structure calculation with SIESTA

First, let's take a look at the input `silicon.fdf` file where the system and the different calculation parameters are specified:

- System name and label:

   ```
   SystemName          bulk silicon
   SystemLabel         si
   ```

- Unit cell information:

   ```
   NumberOfAtoms       2
   NumberOfSpecies     1

   %block ChemicalSpeciesLabel
    1  14  Si
   %endblock ChemicalSpeciesLabel

   LatticeConstant    10.2048933684 Bohr
   %block LatticeVectors
    -0.500  0.000  0.500
     0.000  0.500  0.500
    -0.500  0.500  0.000
   %endblock LatticeVectors

   AtomicCoordinatesFormat  ScaledByLatticeVectors
   %block AtomicCoordinatesAndAtomicSpecies
    -0.125  -0.125  -0.125     1
     0.125   0.125   0.125     1
   %endblock AtomicCoordinatesAndAtomicSpecies
   ```

- Exchange correlation:

   ```
   XC.Functional LDA
   XC.Authors PZ
   ```

- Basis set:

   ```
   PAO.BasisSize          DZ
   ```

- Self-consistent loop:

   ```
   DM.UseSaveDM           true
   DM.MixingWeight        0.3
   DM.NumberPulay         3
   DM.Tolerance           1.d-5
   ```

- Electronic structure parameters:

   ```
   SolutionMethod         diagon
   Diag.Algorithm         expert
   Diag.ParallelOverK     false
   NumberOfEigenStates    10
   OccupationFunction     FD
   ElectronicTemperature  25 meV
   ```

- Brillouin-zone sampling:

   ```
   %block kgrid_Monkhorst_Pack
    8   0   0    0.
    0   8   0    0.
    0   0   8    0.
   %endblock kgrid_Monkhorst_Pack
   ```

- Band block:

   ```
   BandLinesScale ReciprocalLatticeVectors
   %block BandPoints
    0.500000    0.500000    0.500000
    0.475000    0.475000    0.475000
    .
    .
    .
   %endblock BandPoints
   ```

After taking a look at the parameters used in the fdf, we can run the calculation by typing:

```
siesta < silicon.fdf | tee silicon.out
```

SIESTA produces many output files, such as `si.DM` for the density matrix, `si.EIG` for the eigenvalues, `si.FA` for the atomic forces, `Si.ion` and `Si.ion.xml` for the basis orbitals and non-local pseudo-potential, or `si.KP` for the k-points (see SIESTA's documentation for more information). But, for this example, we are interested in the `si.bands` file, which contains the band structure. However, in order to plot the band structure the `si.bands` file needs to be transformed to a format suitable for plotting. Although SIESTA has its own code for that (`gnubands`), in this tutorial we have prepared a script for simplicity:

```
./bands_siesta2nxy.bash > bands_siesta.dat
```

This script transforms the `si.bands` file into a format that can be easily plotted with different programs:

```
xmgrace -nxy bands_siesta.dat
```

or

```
gnuplot -p -e "p for [c=2:*] 'bands_siesta.dat' u 1:c w l lc 'black'"
```


[Back to top :arrow_heading_up:](#using-intw-with-siesta)


### intw's Wannier90 interface

In this second part, we will learn how to interpolate the band structure using intw and `wannier90.x`. In order to constructing the Wannier functions with `wannier90.x` the `si.amn`, `si.mmn` and `si.eig` files are needed. intw's `intw2W90.x` utility can be used used for that, but first of all, we have to transform all the data of the SIESTA calculation to format suitable for intw by running the interface `siesta2intw.x`. So, let's check the `s2intw.in` input file:

```
$ cat s2intw.in
&inputpp
  outdir="./"
  prefix="si"
  nk1=8, nk2=8, nk3=8
  cutoff=40.0
  nbnd_initial=1
  nbnd_final=4
  use_sym=.true.
/
```

:heavy_exclamation_mark:WARNING: `s2intw.in` can't have any other name.

The `outdir` and `prefix` variables specify where will be placed the `si.save.intw` directory and the label used in the siesta calculation, respectively. `nk1`, `nk2` and `nk3` indicate the k-mesh where the Fourier transform of the wave functions will be computed, being `cutoff` the plane-wave cut-off for the wave functions. `nbnd_initial` and `nbnd_final` can be specified to reduce the amount of wave functions to be transformed. And finally, the `use_sym` variable is used to reduce the number of k-points by using the symmetry of the system.

:heavy_exclamation_mark:NOTE: `nk1`, `nk2` and `nk3` don't need to match the k-points of the SIESTA calculation. `siesa2intw.x` will execute a normal self-consistent SIESTA calculation, and then, it will run a non self-consistent calculation to compute the wave functions for the k-points specified by `nk1`, `nk2` and `nk3`.

Now we can run `siesta2intw.x` by typing:

```
siesta2intw.x < silicon.fdf | tee siesta2intw.out
```

This creates the `si.save.intw` directory with all the information about the system in a format readable by intw:

```
$ tree si.save.intw
si.save.intw/
├── 1-KBPP.txt
├── crystal.dat
├── gvectors.dat
├── iGlist.dat
├── kpoints.dat
├── wfc00001.dat
├   .
├   .
├   .
└── wfc00029.dat
```

`1-KBPP.txt` contains the pseudo-potential in its Kleyman-Bylander form. `crystal.dat` stores all the information related to the unit cell, symmetries, etc. `wfc00001.dat` to `wfc00029.dat` contain the Fourier transform of the wave functions for each k-point, where the list of k-points is given in `kpoints.dat`. Finally, `gvectors.dat` contains the global list of G-vectors used in the Fourier transform and `iGlist.dat` has the lists of G-vectors for each k-point.

All the data read by intw is written here, and therefore, now we can already run intw in the same way as if we were using Quantum Espresso and `pw2intw.x`.

To create the `si.amn`, `si.mmn` and `si.eig` files we have to run `wannier90.x` in pre-processing mode first to create the nnkp file:

```
wannier90.x -pp si
```

Since the aim of this tutorial is not to show how to use wannier, we will not analyze the `si.win` file, but for the curious, since we have specified the band range from `nbnd_initial=1` to `nbnd_final=4` in `siesta2intw.x` we will use four s-orbitals, centered in the bonds, as our initial guess.

Now, that we have the `si.save.intw` and the `si.nnkp` file, we can run `intw2W90.x` to generate the `si.amn`, etc. But
first, let's check the input file for intw:

```
$ cat intw.in
&input
    mesh_dir = './'
    prefix = 'si'
    nk1 = 8
    nk2 = 8
    nk3 = 8
    TR_symmetry = .false.
/

&intw2W
    intw2W_fullzone = .false.
    intw2W_method = 'CONVOLUTION'
/

&ph
/
```

In the `input` namelist, `mesh_dir` indicates where is placed the `si.save.intw` directory, and `prefix` the label used in SIESTA. `nk1`, `nk2` and `nk3` indicate the k-mesh, which must be the same used for `siesta2intw.x`, and `TR_symmetry` indicates if time-reversal symmetry can be used to obtain the full k-mesh form the irreducible k-points or not. In the `intw2W` namelist, `intw2W_fullzone` indicates wether the full Brillouin zone is present in the `si.save.intw` directory or not, while `intw2W_method` specifies the method used to compute the amm's. Finally, the `ph` namelist is empty in this example, but it is used to specify information about the phonon structure.

Now, we can run `intw2W90.x` by tying:

```
intw2W90.x < intw.in | tee intw2W90.out
```

Which will create `si.amn`, `si.mmn` and `si.eig` files, and finally we can run `wannier90.x`:

```
wannier90.x si
```

This will construct the Wannier functions, and use them to interpolate the band structure. The interpolated band structure can be found in `si_bands.dat` and can be plotted using the `gnuplot` script `si_bands.gnu`. However, in order to compare the interpolated bands with the original bands calculated with SIESTA, for simplicity, we have prepared an script to write the bands in the same format used previously:

```
./bands_wannier2nxy.bash > bands_wannier_intw.dat
```

And therefore, the band structure can be compared directly by

```
xmgrace -nxy bands_siesta.dat -nxy bands_wannier_intw.dat
```

or by

```
gnuplot -p -e "p for [c=2:*] 'bands_siesta.dat' u 1:c w l lc 'black', for [c=2:*] 'bands_wannier_intw.dat' u 1:c w l lc 'red'"
```


[Back to top :arrow_heading_up:](#using-intw-with-siesta)


### SIESTA's Wannier90 interface


SIESTA has its own interface to `wannier90.x`, so `si.mmn`, `si.amn` and `si.eig` files can be generated adding the following lines to the input fdf file:

```
Siesta2Wannier90.WriteMmn true
Siesta2Wannier90.WriteAmn true
Siesta2Wannier90.WriteEig true
```

:heavy_exclamation_mark:NOTE: Remember executing `wannier90.x` in pre-processing mode to create the `si.nnkp` file before executing SIESTA.


[Back to top :arrow_heading_up:](#using-intw-with-siesta)


## Graphene phonon dispersion

intw contains set of utilities to compute phonon structures using SIESTA and finite differences, `siesta2ph`, which consists of three independent utilities: `siesta2ph.x`, `siesta2fc.x` and `siesta2dv.x`.

The first one, `siesta2ph.x`, reads SIESTA's fdf file, creates a super-cell of the system and computes the minimum set of atomic displacements needed to compute the phonon structure (the irreducible displacements), creating all the fdf input files to run SIESTA.

Then, once SIESTA has been executed for all the irreducible displacements, `siesta2fc.x` can be used to compute the interatomic force constants and dynamical matrices of the system, or, optionally, to compute the phonon dispersion along a path, and `siesta2dv.x` can be used to compute the potential induced by the displacement of the atoms from their equilibrium positions, which is an essential quantity to compute electron-phonon matrix elements.


In this example we will learn how to use this utilities to compute the interatomic force constants, the dynamical matrices and the phonon dispersion of graphene, together with the induced potential.

Folder `2-graphene` has all the files needed to run this example:

```
$ ls 2-graphene
C.psf  graphene.fdf  path.dat  siesta2ph.in
```


First of all, let's check the input file for all the three utilities of `siesta2ph`:

```
$ cat siesta2ph.in
&input
    prefix = "graphene.fdf"
    v0dir = "./441/v0/"
    phdir = "./441/0.01/"
    nr1 = 4
    nr2 = 4
    nr3 = 1
    dx = 0.01
    lpm = .true.
    use_sym = .true.
    irreducible_q = .true.
    kpath_file = "path.dat"
    xsf = .true.
    full_xsf = .false.
    dv_precision = "dp"
    verbose = .true.
/
```

`prefix` specifies the actual fdf file of the system being studied. `v0dir` and `phdir` can be used to properly organize all the files generated by `siesta2ph.x`. `v0dir` and `phdir` specify where will be located the fdf of the super-cell used in the phonon structure calculation, and where will be created the fdf files for the irreducible displacements, respectively. `nr1`, `nr2` and `nr3` indicate the size of the super-cell to be used, while `dx` is the amplitude used for the atomic displacements in Bohr atomic units. `lpm` can be used to displace the atoms also in the negative direction in order to use a more precise second order centered finite difference formula. `use_sym` is used to specify whether symmetries will be used to calculate irreducible displacements or whether all atoms will be displaced in all directions. On the other hand, `irreducible_q` is used to indicate if the dynamical matrices and induced potentials will be computed for the whole q-points commensurate with `nr1`, `nr2` and `nr3`, or just for the irreducible q-points. The `kpath_file` variable designates the name of the file where the k-path for computing the phonon dispersion is stored. While `xsf` indicates if the induced potentials are stored in a 2D XSF file to be plotted by XCrysDen, `full_xsf` determines if the xsf files are stored only for the irreducible displacements, or for all possible atomic displacements. Additionally, `dv_precision` specifies the precision used to compute the induced potentials, and can be used to reduce the amount of disk space used. And finally, `verbose` specifies the level of verbosity of the output.

1. The first step to compute the phonon structure with `siesta2ph` is to create the fdf files for the irreducible displacements using `siesta2ph.x`:

   ```
   siesta2ph.x < siesta2ph.in | tee siesta2ph.out
   ```

   This creates all the directory structure and the fdf files needed:

   ```
   $ tree 441/
   441/
   ├── 0.01
   │   ├── disp-0001
   │   │   └── supercell-graphene.fdf
   │   ├── disp-0002
   │   │   └── supercell-graphene.fdf
   │   ├── disp-0003
   │   │   └── supercell-graphene.fdf
   │   └── disp-0004
   │       └── supercell-graphene.fdf
   └── v0
         └── supercell-graphene.fdf
   ```

   :heavy_exclamation_mark:NOTE: Notice that even if only two irreducible displacements have been found by `siesta2ph.x`, four displacement files have been created because `lpm=.true.` has been specified in the input.

2. Once that all the fdf files have been created for all the irreducible displacements, the next step is to execute a SIESTA calculation for each of them.

   Usually, it is very convenient to run first the calculation of the super-cell with the atoms in their equilibrium positions (the fdf in `v0`), in order to use the density matrix of this calculation as starting point for the displaced atoms calculations:

   ```
   cp C.psf 441/v0/
   cd 441/v0/
   mpirun -n 4 siesta < supercell-graphene.fdf | tee out
   ```

   And after this, we can execute the calculations for all the displacements by using an script included in `siesta2ph` for this purpose:

   ```
   cd ../../
   run-disps.x --dm 441/v0/g.DM -n 4 -v 441/0.01/
   ```

    You can type `run-disps.x --help` for all the information about the command line arguments of the script.

   :heavy_exclamation_mark:NOTE: Remember to include `DM.UseSaveDM true` in the original fdf file in order to take advantage of the close-to-convergence charge density as a starting point for the calculations.

   :heavy_exclamation_mark:NOTE: Since in this calculations the atoms are displaced from their equilibrium positions, it is important to use `MD.Steps 0` in the input fdf to avoid relaxing the system.

3. Now that all SIESTA calculations have been executed, we already can compute the force constants, dynamical matrices and phonon dispersion with `siesta2fc.x`:

   ```
   siesta2fc.x < siesta2ph.in | tee siesta2fc.out
   ```

   This will generate the interatomic force constants file `fc.dat`, the phonon dispersion file `phonons.dat` and the dynamical matrices files `dyn0001.dat` to `dyn0004.dat` for each q-point in the `phdir` directory.

   To plot the phonon dispersion we can use:

   ```
   xmgrace -nxy 441/0.01/phonons.dat
   ```

   or

   ```
   gnuplot -p -e "p for [c=2:*] '441/0.01/phonons.dat' u 1:c w l lc 'black'"
   ```

4. And finally, the induced potential can be computed by `siesta2dv.x`:

   ```
   siesta2dv.x < siesta2ph.in | tee siesta2dv.out
   ```

   This creates the folder `dvscf` inside `phdir`, where induced potentials for each q-point are stored in files `dvscf0001.dat` to `dvscf0004.dat`. Since `xsf=.true.` was specified, the xsf files are also created in `phdir`. This files can be plotted with XCrySDen:

   ```
   xcrysden --xsf 441/0.01/dVR_ia001_id1.xsf
   ```

   `dVR_ia001_id1.xsf` file contains the potential induced in the super-cell by the displacement of the first atom along the first Cartesian direction $x$: $\partial V_{I\alpha} = dV/du_{I}^{\alpha}$ for the atomic index $I=1$ and the Cartesian index $\alpha=x$. In a similar way, `dVq_iq01_ia001_id1.xsf` contains the potential induced in the super-cell, but for the displacement pattern with the periodicity of the first q-point: $\partial V_{I\alpha}(\mathbf{q}) = dV/(u_{1}^{x} e^{i\mathbf{q}\cdot\mathbf{R}})$, being $\mathbf{R}$ a lattice vector. And finally, `dvq_iq01_ia001_id1.xsf` contains the periodic part of the same induced potential: $\partial v_{I\alpha}(\mathbf{q}) = \partial V_{I\alpha}(\mathbf{q}) e^{-i\mathbf{q}\cdot\mathbf{r}}$, which is the same induced potential stored in `dvscf0001.dat`.

   :heavy_exclamation_mark:NOTE: In order to compute the induced potential `SaveTotalPotential T` must be added to de fdf (if you only want to compute the force constants it's not needed).

   :heavy_exclamation_mark:NOTE: At this moment, running the `v0` calculation is mandatory, as some of its output is read by `siesta2dv.x`

[Back to top :arrow_heading_up:](#using-intw-with-siesta)


## Copper electron-phonon matrix elements

In this third example we will show how to calculate the electron-phonon matrix elements for copper using intw's `ep_melements.x` utility. For that, first of all, we will compute the electronic and phonon structures of cooper using SIESTA and `siesta2ph`. In a second step, we will transforming all the calculation data to a format readable by intw with `siesta2intw.x`. And finally, we will execute `ep_melements.x` to compute the electron-phonon matrix elements.

Folder `3-copper` has all the files needed to run the example:

```
$ ls
copper.fdf  Cu.psf  intw.in  s2intw.in  siesta2ph.in
```

1. The first step for a calculation of this type would be computing the electronic structure of the system to be studied:

   ```
   mpirun -n 4 siesta < copper.fdf | tee copper.fdf
   ```

   :heavy_exclamation_mark:NOTE: Before continuing with the calculation of the phonon structure, users should ensure that the parameters used in the fdf yield an electronic structure that is accurate enough for the requirements of their calculations.

   After the electronic structure, we will calculate the dynamical matrices and the induced potentials with `siesta2ph.x`, `siesta2fc.x` and `siesta2dv.x`. Following the procedure that we have seen in the previous example:

   - Create the input files of the irreducible displacements:

     ```
     siesta2ph.x < siesta2ph.in | tee siesta2ph.out
     ```

   - Execute all SIESTA calculations:

     ```
     cp Cu.psf 222/v0/
     cd 222/v0/
     mpirun -n 4 siesta < supercell-copper.fdf | tee out
     cd ../../
     run-disps.x --dm 222/v0/cu.DM -n 4 -v 222/0.02/
     ```

   - Compute the dynamical matrices and induced potentials:

     ```
     siesta2fc.x < siesta2ph.in | tee siesta2fc.out
     siesta2dv.x < siesta2ph.in | tee siesta2dv.out
     ```

   After this steps, we already have computed all the information about the electronic and phonon structures that we will need to calculate the electron-phonon matrix elements.

2. Then, the second step is to run `siesta2intw.x` to store all the data in a format that intw can read. However, in this case, together with the electronic structure, we also need the phonon structure. Therefore, let's check the flags that we have to add in `s2intw.in`:

   ```
   $ cat s2intw.in
   &inputpp
     outdir="./"
     prefix="cu"
     nk1=4, nk2=4, nk3=4
     nbnd_final=10
     cutoff=80.0
     phonons=.true.
     phdir="./222/0.02/"
     nqirr=3
   /
   ```

   `phonons=.true.` indicates that also the phonon structure has to be copied to the `cu.save.intw` directory, `phdir` specifies where are placed the dynamical matrices and induced potentials, and `nqirr` indicates the number of irreducible q-points to be copied.

   So, let's run `siesta2intw.x` by typing:

   ```
   siesta2intw.x < copper.fdf | tee siesta2intw.out
   ```

   This will create the same files that we have seen in the first examples inside `cu.save.intw`, plus `cu.dyn_q*` files for the dynamical matrices, and `cu.dvscf_q*` and `irrq_patterns.dat` files for the induced potentials and the displacement patterns related to them, respectively. Additionally, `qlist.txt` file will also be created, which is needed by intw and contains a list of the irreducible q-points in Cartesian coordinates.

3. Now that we have all the required information in a format that intw can read, we can run intw's `ep_melements.x` utility to compute the electron-phonon matrix elements. But first, let's check `intw.in` file:

   ```
   $ cat intw.in
   &input
     mesh_dir = './'
     prefix = 'cu'
     nk1 = 4
     nk2 = 4
     nk3 = 4
     TR_symmetry = .false.
   /
   &intw2W
     intw2W_fullzone = .false.
     intw2W_method   = 'CONVOLUTION'
   /
   &ph
     nq1 = 2
     nq2 = 2
     nq3 = 2
     nqirr = 3
   /
   ```

   In addition to the parameters that we have already seen in the first example, in this case we also need to specfy the `ph` namelist, where `nq1`, `nq2` and `nq3` indicate the q-mesh, and `nqirr` the number of irreducible q-points.

   So, let's run `ep_melements.x` by typing:

   ```
   ep_melements.x < intw.in | tee ep_melements.out
   ```

   And as you can see, the electron-phonon matrix elements for each q-point have been calculated and saved into files `ep_mat.dat_1` to `ep_mat.dat_8`:

   ```
   $ ls ep_mat.*
   ep_mat.dat_1  ep_mat.dat_2  ep_mat.dat_3  ep_mat.dat_4
   ep_mat.dat_5  ep_mat.dat_6  ep_mat.dat_7  ep_mat.dat_8
   ```

   :heavy_exclamation_mark:NOTE: intw will compute the matrix elements for the entire `nq1, nq2, nq3` q-mesh.

[Back to top :arrow_heading_up:](#using-intw-with-siesta)
