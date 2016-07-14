# NAMD files


### spec

This patch contains all modifications to NAMD 2.11 included in `thstat`.
It also includes two additional changes.

 1. Reporting the alpha-carbon end-to-end distance in every step.
 2. Adding a column (first) of beta to the restart file (due to Justin).


#### alpha-carbon end-to-end distance

The alpha-carbon end-to-end distance is computed in CollectionMaster.C and CollectionMgr.h.
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
