### Force field prepapration

* Download `toppar_c36_jul16.tgz` from http://mackerell.umaryland.edu/charmm_ff.shtml
* Extract `top_all_35_ethers.rtf` and `par_all35_ethers.prm` from the package.
* Add the following lines to the respective sections to `par_all35_ethers.prm`.
These are parameters for TIP3P water molecules.
```
OT   HT    450.000     0.9572 ! ALLOW   WAT
                ! FROM TIPS3P GEOM

HT   OT   HT     55.000   104.5200 ! ALLOW WAT
                ! TIP3P GEOMETRY, ADM JR.

HT           0.0    -0.0460    0.2245                 ! TIP3P water

OT     0.000000  -0.152100     1.768200 ! ALLOW   WAT
                !TIP3P OXYGEN PARAMETERS, adm jr., NBFIX obsolete
```

### butane

* Prepare `butane0.pdb`.
The residues name (columns 18-21) is BUTA.
The atom names (columns 14-16) are C1, C2, C3, and C4.
The hydrogen names are H11, H12, H13, H21, H22, ...

* Prepare `butane.pgn`
```
package require psfgen
topology top_all35_ethers.rtf
segment U {pdb butane0.pdb}
coordpdb butane0.pdb U
writepdb butane1.pdb
writepsf butane1.psf
```

* Run
```
vmd butane0.pdb
```
From the menu, choose Extensions -> Tk Console.
In VMD Tk Console, run
```
source butane.pgn
```
This will produce `butane1.pdb` and `butane1.psf`

* Solvation.
In Tk Console, type
```
package require solvate
solvate butane1.psf butane1.pdb -minmax {{-15, -15, -15} {15, 15, 15}} -o butane
```
This creates a box of 30A x 30A x 30A, centered at the origin.

* Quit VMD. Delete `butane.log` and rename `butane.pdb` as `butane2.pdb`.

* Edit the minimization and equilibration script `min.conf`
```
structure           butane.psf
coordinates         butane2.pdb

set outname         butane3

paraTypeCharmm	    on
parameters          par_all35_ethers.prm

cellBasisVector1    30.0  0.0   0.0
cellBasisVector2     0.0 30.0   0.0
cellBasisVector3     0.0  0.0  30.0
cellOrigin           0.0  0.0   0.0

minimize 500
run 10000
```

* Run NAMD
```
namd2 +p2 min.conf
```

* Convert the final structure to PDB file,
```
vmd butane2.pdb
```
In VMD, load coordinates `butane3.restart.corr`.
Then save coordinates of the last frame (First: 1, Last: 1, Stride 1) as `butane.pdb`.
If the ouptut PDB file contains multiple structures,
keep only the last structure and delete others.
This will be used as the starting point for simulations.

* Clean up
```
rm butane3* FFTW*
```
Copy `butane.psf`, `butane.pdb` and `par_all35_ethers.prm` to the running directory.
 
