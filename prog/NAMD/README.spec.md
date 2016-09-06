# NAMD files

This file can be viewed online

https://github.com/cheng2zhang/tempering/blob/master/prog/NAMD/README.spec.md

### spec

This patch contains all modifications to NAMD 2.11 included in `thstat`.
It also includes

 * Reporting quantities of a set of special atoms

 * Note on reconstructing distributions using WHAM

#### Reporting a quantity of a set of special atoms


To compute the alpha-carbon end-to-end distance,
set the following options
```
specAtoms             on
specAtomsFile         end2end.log
specAtomsList         CAY|CA|CAT
specAtomsType         end-to-end distance, radius of gyration
specAtomsFreq         4
```
By default `specAtomsFreq` is set to `dcdFrequency`.

To compute the dihedral of butane,
set the following options
```
specAtoms             on
specAtomsFile         dih.log
specAtomsList         C1,C2,C3,C4
specAtomsType         dihedral
```

The log file specified by `specAtomsFile` contains step and the quantity computed from special atoms.
If adaptive tempering is turned on, the log file also contains the temperature and potential energy.

The format of `specAtomsList` is the following. For example
```
specAtomsList         A1,A2,A3|B|C1,C2
```
This defines three groups.
The first group has three types of atoms: A1, A2 and A3.
The second group has one type of atoms: B.
The third group has two types of atoms: C1 and C2.

Atoms in the first group will be searched first and added to the list of special atoms.
Then atoms in the second group will be searched, then atoms in the third group, ...

For example, with "CAY|CA|CAT", and PDB file lists atoms as
CAY, CA, CA, ..., CAT, CA
the N-terminal cap "CAY" will be searched first,
regular alpha-carbon "CA" will be searched next,
the C-terminal cap "CAT" will be searched last.

The actual code in CollectionMaster.C and CollectionMgr.h.
Particularly, the new routines
`CollectionMaster::receiveSpecPositions()`, `CollectionMaster::enqueueSpecPositions()`,
`CollectionMaster::disposeSpecPositions()`, and `CollectionMgr::submitSpecPositions()`.


#### Note on reconstructing distributions using WHAM

For a constant temperature simulation, we can reconstruct
the distribution of the end-to-end distance or the dihedral angle
by running mkhist.py.
```
./mkhist.py --dx=0.01745329252 dih0.log
```
Then the output histogram is `dih0.his`.
The constant 0.01745329252 is equal to Pi/180, i.e., one degree.

With adaptive tempering, we can reconstruct the distribution
at a certain temperature using `mkhist.py`.
Such a reconstruction can be done either by explicit filtering
trajectory frames with the temperature within a certain range.
For the dihedral, the command is
```
./mkhist.py --dx=0.01745329252 --dT=0.5 -T300 --colT=3 dih.log
```
Then the output histogram is `dih.his`.

Alternatively, it can be done using WHAM (weighted histogram analysis method)
which is more efficient in using data
```
./mkhist.py --dx=0.01745329252 --rst=narrow1.rst -T300 --colT=3 --colE=4 dih.log
```
The third column will be used for temperature, the fourth for the potential energy,
the WHAM weights (average energy) will be read from the NAMD restart file `narrow1.rst`.

For more information about the tool `mkhist.py`, please see the help message displayed from the command
```
./mkhist.py -h
```
