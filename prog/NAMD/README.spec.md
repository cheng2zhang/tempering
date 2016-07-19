# NAMD files


### spec

This patch contains all modifications to NAMD 2.11 included in `thstat`.
It also includes two additional changes.

 1. Reporting the alpha-carbon end-to-end distance in every step.
 2. Adding a column (first) of beta to the restart file (due to Justin).


#### alpha-carbon end-to-end distance

The alpha-carbon end-to-end distance is computed.
To use this feature, set `specAtomsOn` and `specAtomsFreq`
```
specAtoms on
specAtomsFreq 4
```
By default `specAtomsFreq` is set to `dcdFrequency`.

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
