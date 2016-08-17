# NAMD files


### spec

This patch contains all modifications to NAMD 2.11 included in `thstat`.
It also includes an additional change.

 * Reporting a quantity of a set of special atoms


#### Reporting a quantity of a set of special atoms


To compute the alpha-carbon end-to-end distance,
set the following options
```
specAtoms             on
specAtomsList         CAY|CA|CAT
specAtomsType         end-to-end distance
specAtomsFreq         4
```
By default `specAtomsFreq` is set to `dcdFrequency`.

To compute the dihedral of butane,
set the following options
```
specAtoms             on
specAtomsList         CT3,CT2
specAtomsType         dihedral
```

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

For example, with "CAY|CA|CAT", the N-terminal cap "CAY" will be searched first,
regular alpha-carbon "CA" will be searched next,
the C-terminal cap "CAT" will be searched last.

With "CT3,CT2", the two atoms "CT3" and "CT2" are searched with equal priority
as they are within the same group.
Thus the four atoms in butane, CT3, CT2, CT2, CT3 will be searched sequentially.



The actual code in CollectionMaster.C and CollectionMgr.h.
Particularly, the new routines
`CollectionMaster::receiveSpecPositions()`, `CollectionMaster::enqueueSpecPositions()`,
`CollectionMaster::disposeSpecPositions()`, and `CollectionMgr::submitSpecPositions()`.

## Apply patches

http://www.thegeekstuff.com/2014/12/patch-command-examples/

```
make spec.patch
```

To use the patch
```
patch -b -p3 < spec.patch
