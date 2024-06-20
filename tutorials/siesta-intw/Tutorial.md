

# Using intw with SIESTA


### Table of contents
- [Silicon band structure with SIESTA and Wannier](#silicon-band-structure-with-siesta-and-wannier)
   - [Band structure Calculation with SIESTA](#band-structure-calculation-with-siesta)
   - [SIESTA's Wannier90 interface](#siestas-wannier90-interface)
   - [intw's Wannier90 interface](#intws-wannier90-interface)
- [Graphene phonon dispersion](#graphene-phonon-dispersion)
- [Copper electron-phonon matrix elements](#copper-electron-phonon-matrix-elements)


In this tutorial we will use the 4.1.5 version of SIESTA.

SIESTA uses PSF pseudo-potentials. (In the recently released 5.0.0 version PSML pseudo-potential from pseudo-dojo can be used).

The input for the calculation is given in Flexible Data Format (FDF).


## Silicon band structure with SIESTA and Wannier

In this example, first of all, we will do a small introduction to a band structure calculation with SIESTA. In second place, we will see how is used SIESTA's interface to `wannier90.x`, comparing the wannierized band structure with the original SIESTA calculation. And finally, we will show how to use the interface `siesta2intw.x`, doing the wannierization using intw's `intw2W90.x` utility and, comparing the band structure wannierized with intw and the SIESTA's one.

All files to run this example can be found in the `1-silicon` directory.

```
$ ls
intw.in                 s2intw.in    si.win
bands_siesta2nxy.bash   silicon.fdf
bands_wannier2nxy.bash  Si.psf
```

[Back to top :arrow_heading_up:](#using-intw-with-siesta)

### Band structure calculation with SIESTA

First, lets check the `silicon.fdf` input file:

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
   # PAO.EnergyShift        300 meV
   ```

- Self-consistent loop:

   ```
   DM.UseSaveDM           true
   DM.MixingWeight        0.3
   DM.NumberPulay         3
   DM.Tolerance           1.d-4
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

After this small introduction of some of the fdf parameters we will run the calculation:

```
siesta < silicon.fdf | tee silicon.out
```

SIESTA produces many output files:

```
0_NORMAL_EXIT
BASIS_ENTHALPY
BASIS_HARRIS_ENTHALPY
CLOCK
fdf-40145.log
FORCE_STRESS
INPUT_TMP.40124
MESSAGES
NON_TRIMMED_KP_LIST
OUTVARS.yml
PARALLEL_DIST
si.bands
si.bib
si.BONDS
si.BONDS_FINAL
si.DM
si.EIG
si.FA
Si.ion
Si.ion.xml
si.KP
si.ORB_INDX
si.STRUCT_OUT
si.XV
```

Files that might be important for us are `si.bands` (band structure), `si.DM` (density matrix), `si.EIG` (eigenvalues), `si.FA` (forces), `Si.ion` and `Si.ion.xml` (basis orbitals and non-local PP) and `si.KP` (k-points).

Now that we know what contains each output file, let's plot the band structure. SIESTA has it's own code to transform `si.bands` into a format suitable for plotting, but for this tutorial I have prepared a script to do it in a faster way using `xmgrace`:

```
./bands_siesta2nxy.bash > bands_siesta.dat
xmgrace -nxy bands_siesta.dat
```

[Back to top :arrow_heading_up:](#using-intw-with-siesta)


### SIESTA's Wannier90 interface

```
wannier90.x -pp si
```

To run `wannier90.x`, we just have to add the following lines to the fdf:

```
Siesta2Wannier90.WriteMmn true
Siesta2Wannier90.WriteAmn true
Siesta2Wannier90.WriteEig true
```

Then, run SIESTA again:

```
siesta < silicon.fdf | tee silicon.out
```

And rename `si.eigW`:

```
mv si.eigW si.eig
```

Finally, run `wannier90.x`:

```
wannier90.x si
```

and compare the bands:

```
./bands_wannier2nxy.bash > bands_wannier_siesta.dat
xmgrace -nxy bands_siesta.dat bands_wannier_siesta.dat
```

[Back to top :arrow_heading_up:](#using-intw-with-siesta)

### intw's Wannier90 interface

First of all, we have to run `siesta2intw.x`. For that, let's check the `s2intw.in` file:

```
$ cat s2intw.in
&inputpp
  outdir="./"
  prefix="si"
  nk1=8, nk2=8, nk3=8
  cutoff=40.0
  nbnd_final=4
  use_sym=.true.
/
```

:heavy_exclamation_mark:WARNING: `s2intw.in` can't have any other name.

:heavy_exclamation_mark:NOTE: `nk1`, `nk2` and `nk3` don't need to match the k-points of the SIESTA calculation. `siesa2intw.x` will execute a normal self-consistent SIESTA calculation, and it will run a non self-consistent at the end for the k-mesh specified by `nk1`, `nk2` and `nk3`.

Now we can run `siesta2intw.x`:

```
siesta2intw.x < silicon.fdf | tee siesta2intw.out
```

This creates the `si.save.intw` directory:

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

All the data read by intw is written here, and therefore, now we can already run intw in the same way as if we were using Quantum Espresso.

First, lets check the input file for intw:

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

To create the `si.amn`, `si.mmn` and `si.eig` files:

```
rm si.amn si.mmn si.eig
intw2W90.x < intw.in | tee intw2W90.out
```

And finally run `wannier90.x`:

```
rm si_band.dat
wannier90.x si
```

and compare the bands:

```
./bands_wannier2nxy.bash > bands_wannier_intw.dat
xmgrace bands_wannier_siesta.dat bands_wannier_intw.dat
gvimdiff bands_wannier_siesta.dat bands_wannier_intw.dat
```

[Back to top :arrow_heading_up:](#using-intw-with-siesta)


## Graphene phonon dispersion

:heavy_exclamation_mark:NOTE: It is important to first relax the structure, and then use `MD.Steps 0` for the relaxed structure.

:heavy_exclamation_mark:NOTE: In order to compute the induced potential `SaveTotalPotential T` must be added to de fdf (if you only want to compute the force constants it's not needed).

Folder `2-graphene` has all the files needed to compute the phonon dispersion of graphene:

```
$ ls
C.psf  graphene.fdf  path.dat  siesta2ph.in
```

Let's check the `siesta2ph.in` file:

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


1. The first step is to create the fdf files for the irreducible displacements:

   ```
   siesta2ph.x < siesta2ph.in | tee siesta2ph.out
   ```

   This creates all the directory structure and the fdf files:

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

2. The next step is to run all the supercell calculations.

   Usually, it is very convenient to run first the calculation of the supercell with the atoms in the equilibrium position (`v0`), in order to use the density matrix of this calculation as starting point for the displaced atoms calculations:

   ```
   cp C.psf 441/v0/
   cd 441/v0/
   mpirun -n 4 siesta < supercell-graphene.fdf | tee out
   ```

   And after this, I wrote a bash script to run all the calculations easily:

   ```
   cd ../../
   run-disps.x --dm 441/v0/g.DM -n 4 -v 441/0.01/
   ```

   :heavy_exclamation_mark:NOTE: At this moment, running the `v0` calculation is mandatory, as some of it's output is read by `siesta2dv.x`

3. And now we already can compute the force constants, dynamical matrices and phonon dispersion with `siesta2fc.x`:

   ```
   siesta2fc.x < siesta2ph.in | tee siesta2fc.out
   ```

   This will generate:

   ```
   $ ls 441/0.01
   disp-0001
   disp-0002
   disp-0003
   disp-0004
   dyn0001.dat
   dyn0002.dat
   dyn0003.dat
   dyn0004.dat
   g.FC_1
   g.FC_2
   g.FC_3
   g.FC_4
   g.FC_5
   g.FC_6
   g.FC_7
   g.PH
   ```

   `g.PH` contains the phonon dispersion:

   ```
   xmgrace -nxy 441/0.01/g.PH
   ```

4. And finally, the induced potential can be computed by `siesta2dv.x`:

   ```
   siesta2dv.x < siesta2ph.in | tee siesta2dv.out
   ```

   The xsf files created can be plotted with XCrySDen:

   ```
   xcrysden --xsf 441/0.01/dVR_ia001_id1.xsf
   ```


[Back to top :arrow_heading_up:](#using-intw-with-siesta)


## Copper electron-phonon matrix elements

Folder `3-copper` has all the files needed to compute the electron-phonon matrix elements of copper:

```
$ ls
copper.fdf  Cu.psf  intw.in  s2intw.in  siesta2ph.in
```

1. The first step is to compute the dynamical matrices and the induced potentials with `siesta2ph.x`, `siesta2fc.x` and `siesta2dv.x`:

   ```
   siesta2ph.x < siesta2ph.in | tee siesta2ph.out
   cp Cu.psf 222/v0/
   cd 222/v0/
   mpirun -n 4 siesta < supercell-copper.fdf | tee out
   cd ../../
   run-disps.x --dm 222/v0/cu.DM -n 4 -v 222/0.02/
   siesta2fc.x < siesta2ph.in | tee siesta2fc.out
   siesta2dv.x < siesta2ph.in | tee siesta2dv.out
   ```

2. The second step is to run `siesta2intw.x`. However, in this case, together with the electronic structure, we also need the phonon structure. Therefore, let's check the flags that we have to add in `s2intw.in`:

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

   `phonons=.true.` indicates that also the phonon structure has to be copied to the `save.intw` directory,`phdir` specifies where are placed the dynamical matrices and induced potentials, and `nqirr` indicates the number of irreducible q-points to be copied.

   So, let's run `siesta2intw.x`:

   ```
   siesta2intw.x < copper.fdf | tee siesta2intw.out
   ```

   This will create `cu.dvscf_q*`, `cu.dyn_q*` and `irrq_patterns.dat` files inside `cu.save.intw`, and it will also create automatically `qlist.txt`, which is needed by intw.

3. Now that we have all the required information to compute the electron-phonon matrix elements, we can run intw's `ep_melements.x` utility. But first, let's check `intw.in` file:

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

   And now, let's run `ep_melements.x`:

   ```
   ep_melements.x < intw.in | tee ep_melements.out
   ```

   And check that the matrix elements have been calculated:

   ```
   $ ls ep_mat.*
   ep_mat.dat_1  ep_mat.dat_2  ep_mat.dat_3  ep_mat.dat_4
   ep_mat.dat_5  ep_mat.dat_6  ep_mat.dat_7  ep_mat.dat_8
   ```

   :heavy_exclamation_mark:NOTE: intw will compute the matrix elements for the entire `nq1, nq2, nq3` q-mesh.

[Back to top :arrow_heading_up:](#using-intw-with-siesta)
